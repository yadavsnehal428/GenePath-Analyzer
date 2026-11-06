"""Fetch a RefSeq transcript and its coding region (CDS) from NCBI."""
import re
from dataclasses import dataclass
from io import StringIO

import requests

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ACCESSION_RE = re.compile(r"\b([NX][MR]_\d+(?:\.\d+)?)\b")


@dataclass
class Transcript:
    accession: str
    gene: str
    sequence: str
    cds_start: int   # 1-based first base of the start codon
    cds_end: int     # 1-based last base of the stop codon
    protein_id: str


def accession_from_header(header: str) -> str | None:
    """Find a RefSeq transcript accession such as NM_007194.4 in a FASTA header."""
    m = ACCESSION_RE.search(header)
    return m.group(1) if m else None


def parse_genbank(text: str) -> Transcript:
    from Bio import SeqIO
    rec = SeqIO.read(StringIO(text), "genbank")
    cds = next((f for f in rec.features if f.type == "CDS"), None)
    if cds is None:
        raise ValueError(f"{rec.id} has no CDS (is it a non-coding transcript?)")
    gene = cds.qualifiers.get("gene", ["Unknown"])[0]
    prot = cds.qualifiers.get("protein_id", [""])[0]
    return Transcript(rec.id, gene, str(rec.seq).upper(),
                      int(cds.location.start) + 1, int(cds.location.end), prot)


def fetch_transcript(accession: str, email: str = "") -> Transcript:
    params = {"db": "nuccore", "id": accession, "rettype": "gb", "retmode": "text"}
    if email:
        params["email"] = email
    r = requests.get(EUTILS, params=params, timeout=30)
    r.raise_for_status()
    if not r.text.startswith("LOCUS"):
        raise ValueError(f"NCBI returned no GenBank record for {accession}")
    return parse_genbank(r.text)


def find_transcript_by_blast(sequence: str) -> str | None:
    """Fallback when no accession is given: BLAST against human RefSeq mRNAs
    and return the accession of the best hit. Slow (often 30 s or more)."""
    from Bio.Blast import NCBIWWW, NCBIXML
    handle = NCBIWWW.qblast("blastn", "refseq_rna", sequence,
                            entrez_query="Homo sapiens[organism] AND biomol_mrna[prop]",
                            hitlist_size=5)
    rec = NCBIXML.read(handle)
    for aln in rec.alignments:
        acc = accession_from_header(aln.title)
        if acc and acc.startswith("NM_"):
            return acc
    return None
