"""Per-residue tolerance scores (gnomAD density - ClinVar density) built with convolve_v4,
and the lookup that maps them back onto individual variants.

    from tolerance import build_tolerance_columns, score_variants
    psd_seq_df = build_tolerance_columns(psd_seq_df, train_df, name="v4_train")
    val_df     = score_variants(val_df, psd_seq_df, "tolerance_v4_train")

Pass only the TRAINING variants to build_tolerance_columns, so the held-out variants never
contribute to the density they are later scored against.
"""
import numpy as np
import pandas as pd
from convolve import convolve_v4

# the recipe from misc_figures (NLGN3 convolve_v4 cell): disorder < 0.5 -> BLOSUM62, w=31;
# disorder >= 0.5 -> Grantham, w=5; the most severe substitution at a residue wins
V4_PARAMS = dict(window_size=[31, 5], matrix=["blosum62", "grantham"],
                 adjustment_threshold=0.5, aggregate="max")

def normalize(y):
    """Min-max to [0, 1]; a flat track (e.g. a gene with no ClinVar variants) -> all zeros."""
    y = np.asarray(y, float)
    mn, mx = y.min(), y.max()
    return np.zeros_like(y) if mx == mn else (y - mn) / (mx - mn)

def tolerance_track(variants_gnomad, variants_clinvar, disorder, **params):
    """One gene -> (tolerance, gnomad_y, clinvar_y), each one value per residue (0-based)."""
    kw = {**V4_PARAMS, **params}
    disorder = list(disorder)
    L = len(disorder)
    gnomad_y  = np.asarray(convolve_v4(variants_gnomad,  L, adjustment_list=disorder, **kw))
    clinvar_y = np.asarray(convolve_v4(variants_clinvar, L, adjustment_list=disorder, **kw))
    return normalize(gnomad_y) - normalize(clinvar_y), gnomad_y, clinvar_y

def build_tolerance_columns(seq_df, variant_df, name="v4", gnomad_min_count=1,
                            verbose=True, **params):
    """Add tolerance_<name>, gnomad_y_<name>, clinvar_y_<name> columns (one array per gene)
    to a copy of seq_df, using only the variants in variant_df.

    gnomad_min_count: drop gnomAD variants with total_count below this before convolving
                      (matches the count_threshold used in the distribution figure).
    params:           override any of V4_PARAMS, e.g. aggregate="sum".
    """
    by_gene = dict(tuple(variant_df.groupby("gene")))
    empty = variant_df.iloc[:0]

    tol, gnom, clin = [], [], []
    for n, row in enumerate(seq_df.itertuples(), 1):
        d = by_gene.get(row.gene, empty)
        is_gnomad = (d["database"] == "GnomAD") & (d["total_count"] >= gnomad_min_count)
        is_clinvar = d["database"] == "ClinVar"
        t, g, c = tolerance_track(d.loc[is_gnomad, "protein_consequence"].tolist(),
                                  d.loc[is_clinvar, "protein_consequence"].tolist(),
                                  row.disorder_scores, **params)
        tol.append(t); gnom.append(g); clin.append(c)
        if verbose and n % 100 == 0:
            print(f"{n}/{len(seq_df)} genes")

    out = seq_df.copy()
    out[f"tolerance_{name}"] = tol
    out[f"gnomad_y_{name}"]  = gnom
    out[f"clinvar_y_{name}"] = clin
    return out

def score_variants(variant_df, seq_df, column):
    """Look up seq_df[column][gene][position-1] for every variant -> new column on a copy.
    Variants whose gene is missing or whose position is out of range get NaN."""
    tracks = dict(zip(seq_df["gene"], seq_df[column]))
    scores = np.full(len(variant_df), np.nan)
    pos = variant_df["position"].to_numpy()
    for gene, idx in variant_df.groupby("gene").indices.items():   # idx = row numbers
        track = tracks.get(gene)
        if track is None:
            continue
        p = pos[idx] - 1
        ok = (p >= 0) & (p < len(track))
        scores[idx[ok]] = np.asarray(track)[p[ok]]
    out = variant_df.copy()
    out[column] = scores
    return out
