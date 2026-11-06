"""HGVS c. (coding DNA) and p. (protein) descriptions and variant consequences."""
from .variants import Variant, apply_variant

CODON = {}
_bases = "TCAG"
_aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_bases):
    for j, b in enumerate(_bases):
        for k, c in enumerate(_bases):
            CODON[a + b + c] = _aas[16 * i + 4 * j + k]

THREE = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
         "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
         "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
         "Y": "Tyr", "V": "Val", "*": "Ter", "X": "Xaa"}


def translate_until_stop(seq: str) -> tuple[str, bool]:
    """Translate from the first base; stop at the first stop codon.
    Returns (protein without '*', whether a stop codon was reached)."""
    prot = []
    for i in range(0, len(seq) - 2, 3):
        aa = CODON.get(seq[i:i + 3], "X")
        if aa == "*":
            return "".join(prot), True
        prot.append(aa)
    return "".join(prot), False


def c_pos(pos: int, cds_start: int, cds_end: int) -> str:
    """Transcript position -> c. coordinate (cds_end = last base of stop codon)."""
    if pos < cds_start:
        return f"-{cds_start - pos}"
    if pos > cds_end:
        return f"*{pos - cds_end}"
    return str(pos - cds_start + 1)


def c_hgvs(v: Variant, cds_start: int, cds_end: int) -> str:
    c = lambda p: c_pos(p, cds_start, cds_end)
    rng = c(v.start) if v.start == v.end else f"{c(v.start)}_{c(v.end)}"
    if v.kind == "sub":
        return f"c.{c(v.start)}{v.ref}>{v.alt}"
    if v.kind == "del":
        return f"c.{rng}del"
    if v.kind == "dup":
        return f"c.{rng}dup"
    if v.kind == "ins":
        return f"c.{c(v.start)}_{c(v.end)}ins{v.alt}"
    return f"c.{rng}delins{v.alt}"


def protein_effect(v: Variant, ref_tx: str, cds_start: int, cds_end: int) -> tuple[str, str]:
    """Return (p. description, Sequence Ontology consequence term)."""
    # Region checks first
    last = v.end if v.kind != "ins" else v.start
    first = v.start if v.kind != "ins" else v.end
    if last < cds_start:
        return "p.(=)", "5_prime_UTR_variant"
    if first > cds_end:
        return "p.(=)", "3_prime_UTR_variant"

    ref_prot, _ = translate_until_stop(ref_tx[cds_start - 1:cds_end])
    mut_tx = apply_variant(ref_tx, v)
    mut_prot, has_stop = translate_until_stop(mut_tx[cds_start - 1:])

    if mut_tx[cds_start - 1:cds_start + 2] != "ATG":
        return "p.(Met1?)", "start_lost"

    R, M = ref_prot + "*", mut_prot + ("*" if has_stop else "")
    i = next((k for k in range(min(len(R), len(M))) if R[k] != M[k]), None)
    if i is None and len(R) == len(M):
        codon_no = (max(v.start, cds_start) - cds_start) // 3 + 1
        aa = R[codon_no - 1] if codon_no <= len(ref_prot) else "*"
        return f"p.({THREE[aa]}{codon_no}=)", "synonymous_variant"
    if i is None:
        i = min(len(R), len(M))

    ref_aa = THREE[R[i]] if i < len(R) else "?"
    alt_aa = THREE[M[i]] if i < len(M) else "?"
    net = len(v.alt) - len(v.ref) if v.kind != "dup" else len(v.alt)

    if net % 3 != 0:  # frameshift
        if alt_aa == "Ter":
            return f"p.({ref_aa}{i + 1}Ter)", "stop_gained"
        n_ter = (len(mut_prot) - i + 1) if has_stop else None
        tail = f"fsTer{n_ter}" if n_ter else "fsTer?"
        return f"p.({ref_aa}{i + 1}{alt_aa}{tail})", "frameshift_variant"

    if R[i] == "*" and i == len(ref_prot):  # the stop codon was lost
        ext = f"extTer{len(mut_prot) - i}" if has_stop else "extTer?"
        return f"p.(Ter{i + 1}{alt_aa}{ext})", "stop_lost"
    if alt_aa == "Ter":
        return f"p.({ref_aa}{i + 1}Ter)", "stop_gained"

    # In-frame change: trim the common suffix to find the altered stretch
    j = 0
    while j < min(len(R), len(M)) - i and R[len(R) - 1 - j] == M[len(M) - 1 - j]:
        j += 1
    r_part, m_part = R[i:len(R) - j], M[i:len(M) - j]
    three = lambda s: "".join(THREE[a] for a in s)
    if len(r_part) == 1 and len(m_part) == 1:
        return f"p.({ref_aa}{i + 1}{alt_aa})", "missense_variant"
    if not m_part:
        loc = f"{THREE[R[i]]}{i + 1}" if len(r_part) == 1 else f"{THREE[R[i]]}{i + 1}_{THREE[r_part[-1]]}{i + len(r_part)}"
        return f"p.({loc}del)", "inframe_deletion"
    if not r_part:
        return (f"p.({THREE[R[i - 1]]}{i}_{THREE[R[i]]}{i + 1}ins{three(m_part)})",
                "inframe_insertion")
    loc = f"{THREE[R[i]]}{i + 1}" if len(r_part) == 1 else f"{THREE[R[i]]}{i + 1}_{THREE[r_part[-1]]}{i + len(r_part)}"
    return f"p.({loc}delins{three(m_part)})", "protein_altering_variant"
