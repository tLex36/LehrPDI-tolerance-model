"""Binomial density smoothing over a protein, with area-based edge correction and optional
substitution-matrix weighting.

    from convolve import convolve_v4
    y = convolve_v4(["W596T", "R601C", 610], protein_len=800, window_size=5)
    y = convolve_v4(variants, protein_len=800, matrix="grantham")   # weight by severity

Edge handling: the fraction of the kernel's AREA that falls outside the protein is measured
directly, and the surviving contributions are divided by the retained area.  A uniformly
covered sequence therefore returns a flat line, so terminal residues are weighted equally
instead of being penalised by the false zeros beyond the ends.
"""
import re
from functools import lru_cache
from scipy.stats import binom

# Grantham's distance (composition, polarity, molecular volume); 0-215, higher = more radical
_GRANTHAM_ORDER = ['S','R','L','P','T','A','V','G','I','F','Y','C','H','Q','N','K','D','E','M','W']
_GRANTHAM_UPPER = [
 [0,110,145,74,58,99,124,56,142,155,144,112,89,68,46,121,65,80,135,177],
 [0,0,102,103,71,112,96,125,97,97,77,180,29,43,86,26,96,54,91,101],
 [0,0,0,98,92,96,32,138,5,22,36,198,99,113,153,107,172,138,15,61],
 [0,0,0,0,38,27,68,42,95,114,110,169,77,76,91,103,108,93,87,147],
 [0,0,0,0,0,58,69,59,89,103,92,149,47,42,65,78,85,65,81,128],
 [0,0,0,0,0,0,64,60,94,113,112,195,86,91,111,106,126,107,84,148],
 [0,0,0,0,0,0,0,109,29,50,55,192,84,96,133,97,152,121,21,88],
 [0,0,0,0,0,0,0,0,135,153,147,159,98,87,80,127,94,98,127,184],
 [0,0,0,0,0,0,0,0,0,21,33,198,94,109,149,102,168,134,10,61],
 [0,0,0,0,0,0,0,0,0,0,22,205,100,116,158,102,177,140,28,40],
 [0,0,0,0,0,0,0,0,0,0,0,194,83,99,143,85,160,122,36,37],
 [0,0,0,0,0,0,0,0,0,0,0,0,174,154,139,202,154,170,196,215],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,24,68,32,81,40,87,115],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,46,53,61,29,101,130],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,94,23,42,142,174],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,101,56,95,110],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,45,160,181],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,126,152],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,67],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]

@lru_cache(maxsize=None)
def _load_blosum(name):
    """Load a Biopython substitution matrix once, plus its (min, max) for rescaling."""
    from Bio.Align import substitution_matrices
    M = substitution_matrices.load(name)
    return M, float(min(M.values())), float(max(M.values()))

@lru_cache(maxsize=None)
def grantham_matrix() -> dict:
    """{ref: {alt: distance}} -- symmetric, 0 on the diagonal."""
    g = {r: {} for r in _GRANTHAM_ORDER}
    for i, r in enumerate(_GRANTHAM_ORDER):
        for j, rr in enumerate(_GRANTHAM_ORDER):
            d = _GRANTHAM_UPPER[i][j] or _GRANTHAM_UPPER[j][i]
            g[r][rr] = 0 if r == rr else d
    return g

_VARIANT_RE = re.compile(r"^(?:p\.)?([A-Z])(\d+)([A-Z*=])$")

def parse_substitution(v):
    """'W596T' / 'p.W596T' -> ('W', 596, 'T').  An int -> (None, int, None)."""
    if isinstance(v, (int, float)):
        return None, int(v), None
    m = _VARIANT_RE.match(str(v).strip())
    if not m:
        return None, None, None
    return m.group(1), int(m.group(2)), m.group(3)

def substitution_weight(ref, alt, matrix="grantham", normalize=True):
    """Severity weight for one substitution. Higher = more radical.

    matrix: 'grantham'  -> Grantham distance (0-215), /215 when normalize
            'blosum62'  -> (max - score) rescaled so dissimilar substitutions score high
            dict        -> {ref: {alt: value}} used directly
            callable    -> f(ref, alt) -> value
            None        -> every variant weighs 1
    """
    if matrix is None or ref is None or alt is None:
        return 1.0
    if callable(matrix):
        return float(matrix(ref, alt))
    if isinstance(matrix, dict):
        return float(matrix.get(ref, {}).get(alt, 0.0))
    key = str(matrix).lower()
    if key in ("grantham", "granthams"):
        d = grantham_matrix().get(ref, {}).get(alt)
        if d is None: return 1.0
        return d / 215.0 if normalize else float(d)
    if key.startswith("blosum"):
        M, lo, hi = _load_blosum(key.upper())
        try:
            s = float(M[ref][alt])
        except KeyError:
            return 1.0
        sev = hi - s                                   # dissimilar -> large
        return sev / (hi - lo) if normalize else sev
    raise ValueError(f"unknown matrix: {matrix!r}")

def _resolve_per_residue(value, adjustment_list, threshold, protein_len, name):
    """Expand a setting into one value per residue.

    adjustment_list is None   -> the single value is used everywhere
    length-2 [low, high]      -> low where adjustment < threshold, high where >= threshold
    length-protein_len list   -> used as given, one entry per residue
    """
    if adjustment_list is None:
        return [value] * protein_len
    if len(adjustment_list) != protein_len:
        raise ValueError(f"adjustment_list has {len(adjustment_list)} entries, "
                         f"expected protein_len={protein_len}")
    if isinstance(value, (list, tuple)):
        if len(value) == protein_len:
            return list(value)
        if len(value) == 2:
            low, high = value
            return [low if a < threshold else high for a in adjustment_list]
        raise ValueError(f"{name} must be a 2-item [below, at/above] pair or one value per "
                         f"residue when adjustment_list is given; got {len(value)} items")
    return [value] * protein_len            # a single value still applies everywhere


def convolve_v4(variants, protein_len, window_size=5, matrix=None, start=1,
                edge="renormalize", normalize_weights=True, return_weights=False,
                adjustment_list=None, adjustment_threshold=0.5,
                matrix_from="variant", aggregate="sum"):
    """Smoothed variant density along a protein.

    variants:    substitution strings ('W596T', 'p.R601C') and/or integer positions.
                 Strings are required when `matrix` is set, since the weight depends on the
                 reference and substituted residues.
    protein_len: length of the region being smoothed
    window_size: w, the binomial kernel spans 2w+1 residues.  With `adjustment_list`, pass
                 [w_below, w_at_or_above] (e.g. [31, 5]) or one w per residue.
    matrix:      None, 'grantham', 'blosum62', a {ref:{alt:v}} dict, or a callable.  With
                 `adjustment_list`, pass [m_below, m_at_or_above] (e.g. ['blosum62',
                 'grantham']) or one matrix per residue.
    adjustment_list: per-residue values (e.g. AIUPred3 disorder scores) that select which
                 matrix and window apply at each residue, split at `adjustment_threshold`.
    matrix_from: 'variant' (default) weights each variant by the matrix chosen at its OWN
                 residue -- a substitution's severity does not depend on where you measure
                 from.  'output' instead re-weights every variant per output position.
    edge:        'renormalize' (default) divides by the kernel area that landed inside the
                 protein; 'none' leaves the raw convolution, which sags at the termini.
    aggregate:   how several variants at ONE residue combine --
                 'sum'    (default) add their weights: measures variant load
                 'max'    keep the most severe: measures worst-case impact
                 'mean'   average them: measures typical severity
                 'binary' 1 if any variant is present, else 0 -- what convolve_v3 used
                          (see assemble_indexed_residue_list), so pass this to reproduce
                          earlier figures
    """
    windows   = _resolve_per_residue(window_size, adjustment_list, adjustment_threshold,
                                     protein_len, "window_size")
    matrices  = _resolve_per_residue(matrix, adjustment_list, adjustment_threshold,
                                     protein_len, "matrix")

    kernels = {}
    def kernel(w):
        if w not in kernels:
            n = w * 2
            kernels[w] = [binom.pmf(r, n, 0.5) for r in range(n + 1)]
        return kernels[w]

    parsed = []
    for v in variants:
        ref, pos, alt = parse_substitution(v)
        if pos is None:
            continue
        i = pos - start
        if 0 <= i < protein_len:
            parsed.append((i, ref, alt))

    def combine(per_residue):
        """per_residue: {index: [weight, ...]} -> a single value per index."""
        w = [0.0] * protein_len
        for i, vals in per_residue.items():
            if not vals: continue
            if   aggregate == "sum":    w[i] = sum(vals)
            elif aggregate == "max":    w[i] = max(vals)
            elif aggregate == "mean":   w[i] = sum(vals) / len(vals)
            elif aggregate == "binary": w[i] = 1.0
            else: raise ValueError(f"unknown aggregate: {aggregate!r}")
        return w

    # weights keyed on each variant's own residue (the usual case)
    by_residue = {}
    for i, ref, alt in parsed:
        by_residue.setdefault(i, []).append(
            substitution_weight(ref, alt, matrices[i], normalize_weights))
    weights = combine(by_residue)

    out = []
    for i in range(protein_len):
        w = windows[i]
        dist = kernel(w)
        if matrix_from == "output":          # re-weight everything by this residue's matrix
            tmp = {}
            for k, ref, alt in parsed:
                tmp.setdefault(k, []).append(
                    substitution_weight(ref, alt, matrices[i], normalize_weights))
            local = combine(tmp)
        else:
            local = weights
        raw = retained = 0.0
        for j in range(len(dist)):
            k = i - w + j
            if 0 <= k < protein_len:
                raw += local[k] * dist[j]
                retained += dist[j]          # kernel area inside the protein = 1 - f
        out.append(raw / retained if (edge == "renormalize" and retained) else raw)
    return (out, weights) if return_weights else out
