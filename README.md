# GenePath Analyzer

GenePath Analyzer finds the differences between a DNA/cDNA sequence and a human RefSeq transcript and reports each one the way clinical genomics does: in **HGVS notation**, with its predicted effect on the protein.

Version 2.0 is a rebuild of GenePath Analyzer v1.4. Compared with v1.4 it:

- aligns against the **specific RefSeq transcript** (e.g. `NM_007194.4`) instead of the top BLAST hit in `nt`
- uses **affine-gap alignment**, so a substitution is no longer misreported as a deletion plus an insertion
- reports positions in **transcript and coding (c.) coordinates** rather than alignment columns
- **normalises indels** to the most 3′ position and recognises **duplications**, following HGVS rules
- translates the coding region to give the **protein change (p.)** and a **Sequence Ontology consequence** (missense, stop gained, frameshift, in-frame deletion, start/stop lost, UTR …)
- replaces the old count-based "harmfulness" score with real predictions from **Ensembl VEP**: SIFT, PolyPhen-2, ClinVar significance and gnomAD frequency

## Example

`examples/CHEK2_query.fasta` is a CHEK2 cDNA sequence whose header names its transcript (`NM_007194.4`). Load it in the app to see the variant table.

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Tests

The core logic (alignment, variant calling, HGVS, protein effect) and the VEP response parsing are tested offline:

```bash
python -m pytest -q
```

## How it works

1. **Reference**: the RefSeq transcript and its CDS coordinates are fetched from NCBI E-utilities (GenBank format).
2. **Alignment**: Biopython `PairwiseAligner`, global mode, match +2, mismatch −3, gap open −5, gap extend −2, free end gaps for the query (so partial sequences work).
3. **Variant calling**: consecutive differing columns are grouped; indels are shifted 3′; insertions that repeat the preceding bases become duplications.
4. **Annotation**: each variant is applied to the reference transcript, the coding region is translated, and the first changed amino acid determines the p. description and consequence.
5. **Functional annotation** (optional): each HGVS description is sent to the [Ensembl VEP REST API](https://rest.ensembl.org) with RefSeq transcripts enabled. GenePath Analyzer keeps the result for the chosen transcript: VEP impact, SIFT and PolyPhen-2 (missense only), allele-specific ClinVar significance, gnomAD exome/genome allele frequency and the dbSNP rsID.

## Version history

- **v2.0** Transcript-specific alignment, HGVS c./p., 3′ normalisation, Sequence Ontology consequences, Ensembl VEP annotation, offline tests, validation against Mutalyzer
- **v1.4** BLAST-based mutation listing with KEGG pathway lookup

## Validation

See [validation.md](validation.md). HGVS descriptions agree with Mutalyzer 3 for all variants tested so far.

## Limitations

Each variant is annotated independently against the reference (as VEP does), so the combined effect of several variants in the same codon is not modelled. Only human RefSeq transcripts are supported. VEP annotation needs an internet connection and depends on Ensembl's data release; SIFT/PolyPhen are predictions, not evidence of pathogenicity. For research and education only, not for clinical use.

## Author

Snehal Yadav, M.Sc. Bioinformatics, DES Pune University. Originally developed under the mentorship of Dr. Kushagra Kashyap.
