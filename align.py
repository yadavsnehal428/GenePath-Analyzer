"""Pairwise alignment of a query sequence against a reference transcript.

Uses Biopython's PairwiseAligner (fast, C implementation) with affine gap
penalties. End gaps in the query are free, so a query that covers only part
of the transcript is aligned without penalty.

Why not pairwise2.globalxx (used in GenePath v1.4)? globalxx scores matches
only and gives mismatches and gaps zero penalty, so a single substitution is
often reported as one deletion plus one insertion. Affine scoring fixes that.
"""

MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND = 2, -3, -5, -2


def align(ref: str, query: str) -> tuple[str, str]:
    """Return (aligned_ref, aligned_query) with '-' for gaps."""
    try:
        from Bio.Align import PairwiseAligner
    except ImportError:  # pragma: no cover - Biopython is in requirements.txt
        return _align_python(ref, query)

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = MATCH
    aligner.mismatch_score = MISMATCH
    aligner.open_gap_score = GAP_OPEN
    aligner.extend_gap_score = GAP_EXTEND
    aligner.query_end_gap_score = 0  # query may be a fragment of the transcript
    best = aligner.align(ref, query)[0]
    return str(best[0]), str(best[1])


def _align_python(ref: str, query: str) -> tuple[str, str]:
    """Pure-Python Gotoh alignment with the same scoring. Slow (O(n*m));
    used only when Biopython is unavailable, e.g. in quick tests."""
    NEG = float("-inf")
    n, m = len(ref), len(query)
    # M: match/mismatch, X: gap in query (deletion), Y: gap in ref (insertion)
    M = [[NEG] * (m + 1) for _ in range(n + 1)]
    X = [[NEG] * (m + 1) for _ in range(n + 1)]
    Y = [[NEG] * (m + 1) for _ in range(n + 1)]
    M[0][0] = 0
    for i in range(1, n + 1):  # leading gaps in query are free
        X[i][0] = 0
    for j in range(1, m + 1):
        Y[0][j] = GAP_OPEN + (j - 1) * GAP_EXTEND
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = MATCH if ref[i - 1] == query[j - 1] else MISMATCH
            M[i][j] = max(M[i-1][j-1], X[i-1][j-1], Y[i-1][j-1]) + s
            free = j == m  # trailing gaps in query are free
            op = 0 if free else GAP_OPEN
            ex = 0 if free else GAP_EXTEND
            X[i][j] = max(M[i-1][j] + op, X[i-1][j] + ex, Y[i-1][j] + op)
            Y[i][j] = max(M[i][j-1] + GAP_OPEN, Y[i][j-1] + GAP_EXTEND, X[i][j-1] + GAP_OPEN)
    # traceback
    i, j = n, m
    state = max((M[i][j], "M"), (X[i][j], "X"), (Y[i][j], "Y"))[1]
    a, b = [], []
    while i > 0 or j > 0:
        if state == "M":
            s = MATCH if ref[i-1] == query[j-1] else MISMATCH
            prev = M[i][j] - s
            a.append(ref[i-1]); b.append(query[j-1]); i -= 1; j -= 1
            state = "M" if M[i][j] == prev else ("X" if X[i][j] == prev else "Y")
        elif state == "X":
            free = j == m or j == 0
            op = 0 if free else GAP_OPEN
            ex = 0 if free else GAP_EXTEND
            cur = X[i][j]
            a.append(ref[i-1]); b.append("-"); i -= 1
            if i == 0 and j == 0:
                break
            state = "X" if X[i][j] + ex == cur else ("M" if M[i][j] + op == cur else "Y")
        else:
            cur = Y[i][j]
            a.append("-"); b.append(query[j-1]); j -= 1
            state = "Y" if Y[i][j] + GAP_EXTEND == cur else ("M" if M[i][j] + GAP_OPEN == cur else "X")
    return "".join(reversed(a)), "".join(reversed(b))
