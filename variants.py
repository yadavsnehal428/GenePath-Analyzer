"""Call variants from an alignment and normalise them the HGVS way.

All positions are 1-based coordinates on the reference transcript.
"""
from dataclasses import dataclass


@dataclass
class Variant:
    kind: str    # "sub", "del", "ins", "dup", "delins"
    start: int   # first affected transcript position (for ins: base before the insertion)
    end: int     # last affected position (for ins: base after the insertion)
    ref: str     # reference bases removed ("" for ins/dup)
    alt: str     # bases added ("" for del)


def call_variants(aln_ref: str, aln_query: str, ref_seq: str) -> list[Variant]:
    """Walk the alignment columns and group consecutive differences."""
    # Only call variants inside the region the query actually covers.
    first = next(i for i, c in enumerate(aln_query) if c != "-")
    last = len(aln_query) - 1 - next(i for i, c in enumerate(reversed(aln_query)) if c != "-")

    variants, block_ref, block_alt, block_start = [], "", "", None
    pos = sum(1 for c in aln_ref[:first] if c != "-")  # ref bases before coverage

    def flush():
        nonlocal block_ref, block_alt, block_start
        if block_start is not None:
            variants.append(_make_variant(block_start, block_ref, block_alt, ref_seq))
        block_ref, block_alt, block_start = "", "", None

    for col in range(first, last + 1):
        r, q = aln_ref[col], aln_query[col]
        if r != "-":
            pos += 1
        if r == q:
            flush()
            continue
        if block_start is None:
            # deletion/substitution: this ref base; insertion: the base before it
            block_start = pos
        elif not block_ref and r != "-":
            block_start = pos  # block began with an insertion, now has ref bases
        block_ref += "" if r == "-" else r
        block_alt += "" if q == "-" else q
    flush()
    # Adjacent substitutions are kept as separate SNVs (most common case);
    # mixed gap/mismatch blocks become delins.
    out = []
    for v in variants:
        if v.kind == "delins" and len(v.ref) == len(v.alt):
            for k, (a, b) in enumerate(zip(v.ref, v.alt)):
                if a != b:
                    out.append(Variant("sub", v.start + k, v.start + k, a, b))
        else:
            out.append(v)
    return out


def _make_variant(start: int, ref: str, alt: str, ref_seq: str) -> Variant:
    if ref and not alt:
        return normalise(Variant("del", start, start + len(ref) - 1, ref, ""), ref_seq)
    if alt and not ref:
        return normalise(Variant("ins", start, start + 1, "", alt), ref_seq)
    if len(ref) == 1 and len(alt) == 1:
        return Variant("sub", start, start, ref, alt)
    return Variant("delins", start, start + len(ref) - 1, ref, alt)


def normalise(v: Variant, ref_seq: str) -> Variant:
    """Shift deletions/insertions as far 3' as possible (HGVS rule) and turn
    insertions that copy the preceding bases into duplications."""
    n = len(ref_seq)
    if v.kind == "del":
        s, e, d = v.start, v.end, v.ref
        while e < n and ref_seq[s - 1] == ref_seq[e]:  # ref_seq is 0-based
            d = d[1:] + ref_seq[e]
            s, e = s + 1, e + 1
        return Variant("del", s, e, d, "")
    if v.kind == "ins":
        p, ins = v.start, v.alt
        while p < n and ref_seq[p] == ins[0]:
            ins = ins[1:] + ins[0]
            p += 1
        k = len(ins)
        if p - k >= 0 and ref_seq[p - k:p] == ins:
            return Variant("dup", p - k + 1, p, "", ins)
        return Variant("ins", p, p + 1, "", ins)
    return v


def apply_variant(ref_seq: str, v: Variant) -> str:
    """Return the transcript sequence carrying this single variant."""
    if v.kind in ("sub", "del", "delins"):
        return ref_seq[: v.start - 1] + v.alt + ref_seq[v.end:]
    if v.kind == "ins":
        return ref_seq[: v.start] + v.alt + ref_seq[v.start:]
    if v.kind == "dup":
        return ref_seq[: v.end] + v.alt + ref_seq[v.end:]
    raise ValueError(v.kind)
