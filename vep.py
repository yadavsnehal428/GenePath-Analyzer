"""Functional annotation with the Ensembl VEP REST API.

For each HGVS description VEP returns the consequence on every overlapping
transcript plus "colocated variants" (known variants at the same position).
From that we keep, for the RefSeq transcript the user chose:
  - SIFT and PolyPhen-2 predictions (missense variants only)
  - VEP impact class (HIGH / MODERATE / LOW / MODIFIER)
  - ClinVar clinical significance for the same allele, via Ensembl
  - gnomAD allele frequency (exomes and genomes)
  - dbSNP rsID

Docs: https://rest.ensembl.org/documentation/info/vep_hgvs_post
"""
import time

import requests

SERVER = "https://rest.ensembl.org"
PARAMS = {"refseq": 1, "hgvs": 1, "canonical": 1, "mane": 1,
          "af_gnomade": 1, "af_gnomadg": 1, "variant_class": 1}
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}


MINIMAL_PARAMS = {"refseq": 1}


class VEPError(Exception):
    pass


TRANSIENT = {500, 502, 503, 504}  # Ensembl overloaded: worth retrying


def _post(hgvs: list[str], params: dict, retries: int = 5) -> list[dict]:
    """POST to VEP, retrying with back-off when Ensembl is busy (429/5xx)."""
    wait = 1.0
    for attempt in range(retries):
        try:
            r = requests.post(f"{SERVER}/vep/human/hgvs", headers=HEADERS,
                              params=params, json={"hgvs_notations": hgvs}, timeout=60)
        except requests.RequestException as e:
            if attempt == retries - 1:
                raise VEPError(f"network error: {e}")
            time.sleep(wait); wait *= 2
            continue
        if r.status_code == 429:
            time.sleep(float(r.headers.get("Retry-After", wait))); wait *= 2
            continue
        if r.status_code in TRANSIENT and attempt < retries - 1:
            time.sleep(wait); wait *= 2
            continue
        if not r.ok:
            try:
                msg = r.json().get("error", r.text)
            except ValueError:
                msg = r.text
            raise VEPError(f"HTTP {r.status_code}: {str(msg)[:300]}")
        return r.json()
    raise VEPError("Ensembl is busy (rate limited or unavailable); try again in a minute")


def query_vep(hgvs_list: list[str]) -> tuple[dict[str, dict], dict[str, str]]:
    """Return ({hgvs: VEP record}, {hgvs: error message}).

    All variants go in one request first (fewer requests = less chance of
    hitting a busy server). If that fails, each variant is retried on its
    own, first with the full option set and then with only refseq=1."""
    records, errors = {}, {}
    try:
        for rec in _post(hgvs_list, PARAMS):
            records[rec.get("input")] = rec
    except VEPError:
        pass
    for h in hgvs_list:
        if h in records:
            continue
        last_err = ""
        for params in (PARAMS, MINIMAL_PARAMS):
            try:
                recs = _post([h], params)
                if recs:
                    records[h] = recs[0]
                    break
                last_err = "empty response"
            except VEPError as e:
                last_err = str(e)
        else:
            errors[h] = last_err
    return records, errors


def _pick_transcript(rec: dict, accession: str) -> dict | None:
    """The consequence for our RefSeq transcript (any version), else MANE/canonical."""
    base = accession.split(".")[0]
    tcs = rec.get("transcript_consequences", [])
    for tc in tcs:
        if tc.get("transcript_id", "").split(".")[0] == base:
            return tc
    for key in ("mane_select", "canonical"):
        for tc in tcs:
            if tc.get(key):
                return tc
    return tcs[0] if tcs else None


def _alt_allele(rec: dict) -> str:
    allele_string = rec.get("allele_string", "")
    return allele_string.split("/")[-1] if "/" in allele_string else ""


def _clinvar(colocated: list[dict], alt: str) -> str:
    """Allele-specific ClinVar significance. clin_sig_allele looks like
    'T:pathogenic;T:likely_pathogenic;C:benign'."""
    sigs = []
    for cv in colocated:
        for item in (cv.get("clin_sig_allele") or "").split(";"):
            allele, _, sig = item.partition(":")
            if sig and (allele == alt or not alt):
                sigs.append(sig)
    if not sigs:  # fall back to position-level significance
        for cv in colocated:
            sigs.extend(cv.get("clin_sig", []))
    # most severe first, so "pathogenic" is not hidden behind "uncertain significance"
    order = ["pathogenic", "likely_pathogenic", "risk_factor", "drug_response",
             "uncertain_significance", "conflicting_interpretations_of_pathogenicity",
             "likely_benign", "benign"]
    rank = lambda x: order.index(x) if x in order else len(order)
    uniq = sorted(dict.fromkeys(sigs), key=rank)
    return ", ".join(s.replace("_", " ") for s in uniq)


def _gnomad(colocated: list[dict], alt: str) -> float | None:
    best = None
    for cv in colocated:
        freqs = (cv.get("frequencies") or {}).get(alt, {})
        for key in ("gnomade", "gnomadg"):
            if key in freqs:
                best = max(best or 0.0, float(freqs[key]))
    return best


def summarise(rec: dict | None, accession: str, error: str = "") -> dict:
    """Flatten one VEP record into table columns."""
    empty = {"Impact": "", "SIFT": "", "PolyPhen": "", "ClinVar": "",
             "gnomAD AF": None, "dbSNP": "", "VEP transcript": ""}
    if not rec:
        return {**empty, "VEP note": error or "not annotated"}
    tc = _pick_transcript(rec, accession) or {}
    colocated = rec.get("colocated_variants", [])
    alt = _alt_allele(rec)
    fmt = lambda pred, score: f"{pred.replace('_', ' ')} ({score:.2f})" if pred else ""
    rsids = [cv["id"] for cv in colocated if str(cv.get("id", "")).startswith("rs")]
    return {
        "Impact": tc.get("impact", ""),
        "SIFT": fmt(tc.get("sift_prediction", ""), tc.get("sift_score", 0)),
        "PolyPhen": fmt(tc.get("polyphen_prediction", ""), tc.get("polyphen_score", 0)),
        "ClinVar": _clinvar(colocated, alt),
        "gnomAD AF": _gnomad(colocated, alt),
        "dbSNP": ", ".join(rsids),
        "VEP transcript": tc.get("transcript_id", ""),
        "VEP note": "",
    }


def annotate_rows(rows: list[dict], accession: str) -> list[dict]:
    """Add VEP columns to the rows produced by pipeline.annotate()."""
    keys = [r["HGVS (transcript)"] for r in rows]
    records, errors = query_vep(keys)
    return [{**r, **summarise(records.get(k), accession, errors.get(k, ""))}
            for r, k in zip(rows, keys)]
