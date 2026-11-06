"""Offline tests on a small synthetic transcript (no internet needed).

Transcript layout (1-based):
  1-5    5' UTR      GGCAC
  6-35   CDS         ATG AAA CGT TTA GAA AAG CTG GGA CCC TGA
                      M   K   R   L   E   K   L   G   P  *
  36-44  3' UTR      CCATTTGGG
"""
import pytest

from align import _align_python
from hgvs import c_hgvs, protein_effect
from reference import Transcript, accession_from_header
from variants import call_variants
import pipeline

UTR5, CDS, UTR3 = "GGCAC", "ATGAAACGTTTAGAAAAGCTGGGACCCTGA", "CCATTTGGG"
REF = UTR5 + CDS + UTR3
TX = Transcript("NM_TEST.1", "TEST", REF, 6, 35, "NP_TEST.1")


@pytest.fixture(autouse=True)
def pure_python_aligner(monkeypatch):
    monkeypatch.setattr(pipeline, "align", _align_python)


def mutate(pos, ref, alt):
    """Replace ref (starting at 1-based pos) with alt."""
    assert REF[pos - 1:pos - 1 + len(ref)] == ref
    return REF[:pos - 1] + alt + REF[pos - 1 + len(ref):]


def run(query):
    return pipeline.annotate(TX, query)


def test_no_variants():
    assert run(REF) == []


def test_missense():
    # codon 3 CGT (Arg) -> GGT (Gly): c.7C>G
    rows = run(mutate(12, "C", "G"))
    assert rows[0]["c."] == "c.7C>G"
    assert rows[0]["p."] == "p.(Arg3Gly)"
    assert rows[0]["Consequence"] == "missense_variant"


def test_synonymous():
    # codon 3 CGT -> CGC, still Arg
    rows = run(mutate(14, "T", "C"))
    assert rows[0]["c."] == "c.9T>C"
    assert rows[0]["p."] == "p.(Arg3=)"


def test_nonsense():
    # codon 7 CTG? use codon 6 AAG (Lys) -> TAG (stop): c.16A>T
    rows = run(mutate(21, "A", "T"))
    assert rows[0]["c."] == "c.16A>T"
    assert rows[0]["p."] == "p.(Lys6Ter)"
    assert rows[0]["Consequence"] == "stop_gained"


def test_frameshift_deletion_is_shifted_3prime():
    # delete one A from the AAA run (codon 2): HGVS puts it at the last A
    rows = run(mutate(9, "A", ""))
    assert rows[0]["c."] == "c.6del"
    assert rows[0]["Consequence"] == "frameshift_variant"
    # ATG AAC GTT TAG -> Met Asn Val Ter: Lys2 changes first, stop 3 codons on
    assert rows[0]["p."] == "p.(Lys2AsnfsTer3)"


def test_inframe_deletion():
    # delete codon 4 TTA (Leu)
    rows = run(mutate(15, "TTA", ""))
    assert rows[0]["c."] == "c.10_12del"
    assert rows[0]["p."] == "p.(Leu4del)"
    assert rows[0]["Consequence"] == "inframe_deletion"


def test_duplication():
    # insert an extra A in the AAA run -> c.6dup
    rows = run(REF[:11] + "A" + REF[11:])
    assert rows[0]["c."] == "c.6dup"
    assert rows[0]["Type"] == "duplication"


def test_utr_variants():
    rows = run(mutate(2, "G", "T"))
    assert rows[0]["c."] == "c.-4G>T"
    assert rows[0]["Consequence"] == "5_prime_UTR_variant"
    rows = run(mutate(37, "C", "G"))
    assert rows[0]["c."] == "c.*2C>G"
    assert rows[0]["Consequence"] == "3_prime_UTR_variant"


def test_start_lost():
    rows = run(mutate(6, "A", "C"))
    assert rows[0]["Consequence"] == "start_lost"


def test_stop_lost():
    # stop codon TGA (33-35) -> TGG (Trp)
    rows = run(mutate(35, "A", "G"))
    assert rows[0]["Consequence"] == "stop_lost"
    assert rows[0]["p."].startswith("p.(Ter10Trpext")


def test_partial_query_covers_part_of_transcript():
    # query = CDS only, with one missense; UTRs missing must not be called as deletions
    q = mutate(12, "C", "G")[5:35]
    rows = run(q)
    assert [r["c."] for r in rows] == ["c.7C>G"]


def test_accession_from_header():
    h = ">NM_007194.4 Homo sapiens checkpoint kinase 2 (CHEK2), transcript variant 1, mRNA"
    assert accession_from_header(h) == "NM_007194.4"


def test_invalid_characters_rejected():
    with pytest.raises(ValueError):
        pipeline.clean_dna("ATGXXT")
