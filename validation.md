# Validation

GenePath Analyzer v2's HGVS descriptions are checked against **Mutalyzer 3** (https://mutalyzer.nl), the reference HGVS name checker maintained by Leiden University Medical Center.

## Round 1: CHEK2 test sequence (6 Nov 2026)

Input: `CHEK2_query.fasta`, aligned to RefSeq transcript `NM_007194.4` (CHEK2; coding region 59–1690, 543 amino acids, protein `NP_009125.1`).

| # | GenePath Analyzer c. | GenePath Analyzer p. | Mutalyzer c. (normalised) | Mutalyzer p. | Consequence | Match |
|---|---|---|---|---|---|---|
| 1 | `c.327_329delinsT` | `p.(Asn110Ter)` | `c.327_329delinsT` | `p.(Asn110Ter)` | stop_gained | ✅ |
| 2 | `c.676_677del` | `p.(Leu226TrpfsTer18)` | `c.676_677del` | `p.(Leu226TrpfsTer18)` | frameshift_variant | ✅ |
| 3 | `c.1079A>G` | `p.(Glu360Gly)` | `c.1079A>G` | `p.(Glu360Gly)` | missense_variant | ✅ |

**Result: 3/3 variants agree on both c. and p. descriptions.**

## Round 1b: known pathogenic variant, VEP annotation (6 Nov 2026)

Input: `CHEK2_1100delC.fasta`, a CHEK2 fragment carrying the CHEK2 1100delC founder variant associated with breast cancer risk.

| Check | Expected | GenePath Analyzer | Match |
|---|---|---|---|
| Variant called | c.1100del (transcript position 1158) | c.1100del, position 1158 | ✅ |
| Consequence / impact | frameshift, HIGH | frameshift_variant, HIGH | ✅ |
| dbSNP | rs555607708 | rs555607708 | ✅ |
| ClinVar | pathogenic (with some uncertain-significance submissions) | pathogenic, uncertain significance | ✅ |
| gnomAD AF | present, rare (founder variant) | 2.2 × 10⁻³ | ✅ |

This confirms the Ensembl VEP integration retrieves the correct transcript, rsID, allele-specific ClinVar significance and gnomAD frequency for a real variant.

## Offline unit tests

`test_core.py` and `test_vep.py` check the core logic on a synthetic transcript: missense, synonymous, stop gained, frameshift (with 3′ shifting), in-frame deletion, duplication, 5′ and 3′ UTR, start lost, stop lost, partial query sequences, and invalid input. VEP response parsing (transcript choice, SIFT/PolyPhen formatting, allele-specific ClinVar, gnomAD, retry on a busy server) is tested with a mock response. Run with `python -m pytest -q`.

## Planned

Round 2 will test 20–30 ClinVar variants across several genes, chosen to cover the harder cases: duplications, deletions within repeats (3′ shifting), insertions, UTR and start/stop variants.
