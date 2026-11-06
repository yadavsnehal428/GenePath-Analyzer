"""Parsing tests for genepath.vep using a hand-made record shaped like a VEP
REST response (field names as in the Ensembl documentation). The values are
synthetic test data, not real annotations."""
from genepath.vep import summarise

RECORD = {
    "input": "NM_000000.1:c.100C>T",
    "allele_string": "C/T",
    "transcript_consequences": [
        {"transcript_id": "ENST00000000001", "canonical": 1, "impact": "MODIFIER"},
        {"transcript_id": "NM_000000.2", "impact": "MODERATE",
         "consequence_terms": ["missense_variant"],
         "sift_prediction": "deleterious", "sift_score": 0.01,
         "polyphen_prediction": "probably_damaging", "polyphen_score": 0.998},
    ],
    "colocated_variants": [
        {"id": "rs0000001", "clin_sig_allele": "T:pathogenic;G:benign",
         "clin_sig": ["pathogenic", "benign"],
         "frequencies": {"T": {"gnomade": 0.00012, "gnomadg": 0.0002}}},
        {"id": "COSV00000001"},
    ],
}


def test_picks_our_refseq_transcript_even_if_version_differs():
    s = summarise(RECORD, "NM_000000.1")
    assert s["VEP transcript"] == "NM_000000.2"
    assert s["Impact"] == "MODERATE"


def test_predictions_formatted():
    s = summarise(RECORD, "NM_000000.1")
    assert s["SIFT"] == "deleterious (0.01)"
    assert s["PolyPhen"] == "probably damaging (1.00)"


def test_clinvar_is_allele_specific():
    s = summarise(RECORD, "NM_000000.1")
    assert s["ClinVar"] == "pathogenic"  # G:benign belongs to another allele


def test_gnomad_takes_highest_population_af_and_rsid():
    s = summarise(RECORD, "NM_000000.1")
    assert s["gnomAD AF"] == 0.0002
    assert s["dbSNP"] == "rs0000001"


def test_missing_record():
    s = summarise(None, "NM_000000.1")
    assert s["VEP note"] == "not annotated"


def test_retries_on_503(monkeypatch):
    """A busy server (503) is retried instead of failing straight away."""
    import genepath.vep as vep

    class Resp:
        def __init__(self, code, body):
            self.status_code, self._body, self.headers, self.text = code, body, {}, str(body)
            self.ok = code == 200

        def json(self):
            return self._body

    calls = iter([Resp(503, {}), Resp(503, {}), Resp(200, [RECORD])])
    monkeypatch.setattr(vep.requests, "post", lambda *a, **k: next(calls))
    monkeypatch.setattr(vep.time, "sleep", lambda s: None)
    records, errors = vep.query_vep([RECORD["input"]])
    assert RECORD["input"] in records and not errors


def test_clinvar_most_severe_first():
    rec = {**RECORD, "colocated_variants": [
        {"id": "rs1", "clin_sig_allele": "T:uncertain_significance;T:pathogenic"}]}
    assert summarise(rec, "NM_000000.1")["ClinVar"] == "pathogenic, uncertain significance"
