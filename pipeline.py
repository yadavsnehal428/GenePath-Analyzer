"""Glue: reference + query -> annotated variant table."""
from .align import align
from .hgvs import c_hgvs, protein_effect
from .reference import Transcript
from .variants import call_variants


def clean_dna(text: str) -> str:
    lines = [l for l in text.strip().splitlines() if not l.startswith(">")]
    seq = "".join(lines).upper().replace(" ", "").replace("U", "T")
    bad = set(seq) - set("ACGT")
    if bad:
        raise ValueError(f"Sequence contains invalid characters: {''.join(sorted(bad))}")
    return seq


def annotate(tx: Transcript, query: str) -> list[dict]:
    aln_ref, aln_q = align(tx.sequence, query)
    rows = []
    for v in call_variants(aln_ref, aln_q, tx.sequence):
        c = c_hgvs(v, tx.cds_start, tx.cds_end)
        p, consequence = protein_effect(v, tx.sequence, tx.cds_start, tx.cds_end)
        rows.append({
            "HGVS (transcript)": f"{tx.accession}:{c}",
            "c.": c,
            "p.": p,
            "Consequence": consequence,
            "Type": {"sub": "substitution", "del": "deletion", "ins": "insertion",
                     "dup": "duplication", "delins": "deletion-insertion"}[v.kind],
            "Transcript position": v.start if v.start == v.end else f"{v.start}-{v.end}",
        })
    return rows


def reference_protein(tx: Transcript) -> str:
    from .hgvs import translate_until_stop
    return translate_until_stop(tx.sequence[tx.cds_start - 1:tx.cds_end])[0]
