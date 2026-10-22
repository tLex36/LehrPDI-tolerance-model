"""Significance of the ClinVar vs gnomAD separation in a per-variant score.

Variants are NOT independent: neighbours share a smoothed density, and a few genes contribute
thousands of rows.  The tests below therefore go from naive (variant = unit) to honest
(gene = unit, and finally the whole pipeline re-run under a null):

    variant_level_tests   Welch t / Mann-Whitney / KS on all variants  -- effect sizes are fine,
                          p-values are anti-conservative (n is inflated)
    gene_bootstrap        95% CIs for the mean difference and AUROC, resampling whole genes
    gene_level_test       one paired difference per gene, Wilcoxon signed-rank + sign test
    permutation_test      shuffle ClinVar/gnomAD labels within each gene in TRAINING, rebuild
                          the tolerance tracks, re-score validation -> null distribution

    from significance import significance_report, permutation_test
    report = significance_report(val_df, "tolerance_v4_train")
    perm   = permutation_test(psd_seq_df, train_df, val_df, n_perm=100)

Sign convention everywhere: diff = mean(gnomAD) - mean(ClinVar), so a positive diff means
gnomAD variants sit in more tolerant positions; AUROC = P(gnomAD score > ClinVar score),
so 0.5 is chance and > 0.5 is the expected direction.
"""
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind, ks_2samp, wilcoxon, binomtest
from sklearn.metrics import roc_auc_score


def _arrays(df, score_col):
    d = df.dropna(subset=[score_col])
    is_gnomad = (d["database"] == "GnomAD").to_numpy()
    return d, d[score_col].to_numpy(float), is_gnomad


def variant_level_tests(df, score_col):
    """Every variant treated as an independent observation."""
    _, s, g = _arrays(df, score_col)
    gnom, clin = s[g], s[~g]
    t, p_t = ttest_ind(gnom, clin, equal_var=False)
    u, p_u = mannwhitneyu(gnom, clin, alternative="two-sided")
    ks, p_ks = ks_2samp(gnom, clin)
    pooled_sd = np.sqrt((gnom.var(ddof=1) + clin.var(ddof=1)) / 2)
    return dict(n_gnomad=len(gnom), n_clinvar=len(clin),
                mean_gnomad=gnom.mean(), mean_clinvar=clin.mean(),
                diff=gnom.mean() - clin.mean(), cohen_d=(gnom.mean() - clin.mean()) / pooled_sd,
                auroc=u / (len(gnom) * len(clin)),
                welch_t=t, p_welch=p_t, p_mannwhitney=p_u, ks_D=ks, p_ks=p_ks)


def gene_bootstrap(df, score_col, n_boot=2000, n_boot_auroc=500, seed=0):
    """Percentile CIs from resampling GENES with replacement (cluster bootstrap).

    Each replicate draws as many genes as there are, with replacement, and keeps every variant
    of each drawn gene, so within-gene correlation is carried into the spread.
    p_boot: two-sided bootstrap p for diff = 0 (2 x the fraction of replicates on the far side
    of 0); with n_boot replicates the smallest reportable value is 1/n_boot.
    """
    rng = np.random.default_rng(seed)
    d, s, g = _arrays(df, score_col)
    codes, genes = pd.factorize(d["gene"])
    G = len(genes)

    # per-gene sums / counts -> a replicate is just a weighted sum over genes
    sum_g = np.bincount(codes, weights=np.where(g, s, 0), minlength=G)
    n_g   = np.bincount(codes, weights=g.astype(float), minlength=G)
    sum_c = np.bincount(codes, weights=np.where(~g, s, 0), minlength=G)
    n_c   = np.bincount(codes, weights=(~g).astype(float), minlength=G)

    W = np.stack([np.bincount(rng.integers(G, size=G), minlength=G) for _ in range(n_boot)])
    diffs = (W @ sum_g) / (W @ n_g) - (W @ sum_c) / (W @ n_c)

    rows_by_gene = [np.flatnonzero(codes == k) for k in range(G)]
    aucs = []
    for w in W[:n_boot_auroc]:
        idx = np.concatenate([np.repeat(rows_by_gene[k], w[k]) for k in np.flatnonzero(w)])
        if g[idx].all() or (~g[idx]).all():
            continue
        aucs.append(roc_auc_score(g[idx], s[idx]))

    obs = s[g].mean() - s[~g].mean()
    p_boot = min(1.0, 2 * min((diffs <= 0).mean(), (diffs >= 0).mean()))
    return dict(diff=obs, diff_ci=tuple(np.percentile(diffs, [2.5, 97.5])),
                auroc=roc_auc_score(g, s), auroc_ci=tuple(np.percentile(aucs, [2.5, 97.5])),
                p_boot=max(p_boot, 1 / n_boot), n_genes=G, n_boot=n_boot)


def gene_level_test(df, score_col, min_per_group=5):
    """One number per gene: mean(gnomAD) - mean(ClinVar), using genes with at least
    `min_per_group` variants in each database.  Genes are the unit of replication.

    Wilcoxon signed-rank: are the per-gene differences centred above 0?
    Sign test:            is the fraction of genes with diff > 0 above 1/2?
    """
    d, _, _ = _arrays(df, score_col)
    m = d.groupby(["gene", "database"])[score_col].agg(["mean", "size"]).unstack("database")
    ok = (m[("size", "GnomAD")] >= min_per_group) & (m[("size", "ClinVar")] >= min_per_group)
    per_gene = (m.loc[ok, ("mean", "GnomAD")] - m.loc[ok, ("mean", "ClinVar")]).rename("diff")
    w, p_w = wilcoxon(per_gene, alternative="two-sided")
    k = int((per_gene > 0).sum())
    return dict(n_genes=len(per_gene), median_gene_diff=per_gene.median(),
                frac_genes_positive=k / len(per_gene),
                wilcoxon_W=w, p_wilcoxon=p_w,
                p_sign=binomtest(k, len(per_gene), 0.5).pvalue,
                per_gene=per_gene.sort_values())


def significance_report(df, score_col, n_boot=2000, n_boot_auroc=500, min_per_group=5, seed=0):
    """Run the three cheap analyses and print a summary.  Returns the three dicts."""
    v = variant_level_tests(df, score_col)
    b = gene_bootstrap(df, score_col, n_boot, n_boot_auroc, seed)
    gl = gene_level_test(df, score_col, min_per_group)
    print(f"{score_col}:  gnomAD n={v['n_gnomad']:,}  ClinVar n={v['n_clinvar']:,}  "
          f"genes={b['n_genes']}")
    print(f"  mean gnomAD {v['mean_gnomad']:.3f}  ClinVar {v['mean_clinvar']:.3f}  "
          f"diff {v['diff']:+.3f}  Cohen's d {v['cohen_d']:+.2f}")
    print(f"  [variant-level]  Welch p={v['p_welch']:.1e}  Mann-Whitney p={v['p_mannwhitney']:.1e}"
          f"  KS D={v['ks_D']:.2f}  (anti-conservative)")
    print(f"  [gene bootstrap] diff 95% CI {b['diff_ci'][0]:+.3f} to {b['diff_ci'][1]:+.3f}  "
          f"AUROC {b['auroc']:.3f} (95% CI {b['auroc_ci'][0]:.3f}-{b['auroc_ci'][1]:.3f})  "
          f"p_boot<={b['p_boot']:.1e}")
    print(f"  [gene-level]     {gl['n_genes']} genes, median diff {gl['median_gene_diff']:+.3f}, "
          f"{gl['frac_genes_positive']:.0%} positive  Wilcoxon p={gl['p_wilcoxon']:.1e}  "
          f"sign p={gl['p_sign']:.1e}")
    return dict(variant=v, bootstrap=b, gene_level=gl)


def permutation_test(seq_df, train_df, val_df, n_perm=100, seed=0, verbose=True, **build_params):
    """Full-pipeline null: is the held-out separation larger than tolerance tracks built from
    label-shuffled training data would produce?

    For each permutation the ClinVar/gnomAD labels of the TRAINING variants are shuffled within
    each gene (positions and each gene's ClinVar:gnomAD ratio are kept; only the position <->
    database link is broken), the tracks are rebuilt, and the validation variants -- with their
    real labels -- are re-scored.  ~15-20 s per permutation.

    p = (1 + #null diffs >= observed) / (1 + n_perm), one-sided (gnomAD more tolerant);
    the smallest reportable p is 1/(n_perm+1).

    build_params are passed to build_tolerance_columns and must match the observed build.
    Count filters must be applied to train_df BEFORE calling this: shuffled labels would
    otherwise put ClinVar's NaN total_count onto 'GnomAD' rows, so the count filter is
    switched off inside.
    """
    from tolerance import build_tolerance_columns, score_variants

    def val_stats(train):
        tracks = build_tolerance_columns(seq_df, train, name="perm", gnomad_min_count=-np.inf,
                                         verbose=False, **build_params)
        _, s, g = _arrays(score_variants(val_df, tracks, "tolerance_perm"), "tolerance_perm")
        return s[g].mean() - s[~g].mean(), roc_auc_score(g, s)

    train = train_df.copy()
    train["total_count"] = train["total_count"].fillna(np.inf)   # see docstring
    obs_diff, obs_auc = val_stats(train)

    rng = np.random.default_rng(seed)
    labels = train["database"].to_numpy()
    groups = list(train.groupby("gene").indices.values())        # positional row indices
    null_diff, null_auc = [], []
    for i in range(n_perm):
        shuffled = labels.copy()
        for idx in groups:
            shuffled[idx] = rng.permutation(labels[idx])
        train["database"] = shuffled
        dd, aa = val_stats(train)
        null_diff.append(dd); null_auc.append(aa)
        if verbose and (i + 1) % 10 == 0:
            print(f"{i + 1}/{n_perm} permutations")

    null_diff, null_auc = np.array(null_diff), np.array(null_auc)
    return dict(obs_diff=obs_diff, obs_auroc=obs_auc,
                null_diff_mean=null_diff.mean(), null_diff_95=np.percentile(null_diff, 95),
                null_auroc_mean=null_auc.mean(), null_auroc_95=np.percentile(null_auc, 95),
                p_diff=(1 + (null_diff >= obs_diff).sum()) / (1 + n_perm),
                p_auroc=(1 + (null_auc >= obs_auc).sum()) / (1 + n_perm),
                null_diff=null_diff, null_auroc=null_auc)
