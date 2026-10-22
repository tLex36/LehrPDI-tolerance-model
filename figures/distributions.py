"""ClinVar vs gnomAD score distributions with Mann-Whitney / Welch / KS effect sizes.

    from distributions import plot_score_distribution, compare_train_val
    plot_score_distribution(val_df, "tolerance_v4_train", title="Held-out variants")
    compare_train_val(train_df, val_df, "tolerance_v4_train")
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, ttest_ind, ks_2samp

CLIN_C, GNOM_C = '#0072B2', '#E69F00'   # colorblind-safe blue / orange

def split_scores(df, score_col, count_threshold=1, disorder_range=None, genes=None):
    """-> (clinvar_scores, gnomad_scores) as NaN-free arrays.

    count_threshold: minimum gnomAD total_count
    disorder_range:  (lo, hi) keeps lo <= disorder_scores <= hi; None keeps everything
    genes:           iterable of gene symbols to keep; None keeps everything
    """
    keep = np.ones(len(df), bool)
    if disorder_range is not None:
        lo, hi = disorder_range
        keep &= df['disorder_scores'].between(lo, hi).to_numpy()
    if genes is not None:
        keep &= df['gene'].isin(genes).to_numpy()
    d = df[keep]
    clin = d.loc[d['database'] == 'ClinVar', score_col].to_numpy(float)
    gnom = d.loc[(d['database'] == 'GnomAD') & (d['total_count'] >= count_threshold),
                 score_col].to_numpy(float)
    return clin[~np.isnan(clin)], gnom[~np.isnan(gnom)]

def separation_stats(clin, gnom):
    """AUROC = P(ClinVar score > gnomAD score); < 0.5 means ClinVar scores lower."""
    u, p_u = mannwhitneyu(clin, gnom, alternative='two-sided')
    _, p_t = ttest_ind(clin, gnom, equal_var=False)
    d_ks, p_ks = ks_2samp(clin, gnom)
    return dict(auroc=u / (len(clin) * len(gnom)),
                cohend=(clin.mean() - gnom.mean()) / np.sqrt((clin.var() + gnom.var()) / 2),
                ks_d=d_ks, p_mw=p_u, p_t=p_t, p_ks=p_ks,
                clin_mean=clin.mean(), gnom_mean=gnom.mean(),
                n_clin=len(clin), n_gnom=len(gnom))

def plot_score_distribution(df, score_col, title=None, xlabel=None, ax=None, bins=50,
                            count_threshold=1, disorder_range=None, genes=None):
    """Overlaid ClinVar/gnomAD density histograms with mean lines and a stats box.
    Returns the separation_stats dict."""
    clin, gnom = split_scores(df, score_col, count_threshold, disorder_range, genes)
    s = separation_stats(clin, gnom)

    if ax is None:
        _, ax = plt.subplots(figsize=(9, 5))
    if np.isscalar(bins):
        bins = np.linspace(min(clin.min(), gnom.min()), max(clin.max(), gnom.max()), bins)
    ax.hist(gnom, bins=bins, density=True, alpha=0.55, color=GNOM_C, label=f'gnomAD (n={len(gnom):,})')
    ax.hist(clin, bins=bins, density=True, alpha=0.55, color=CLIN_C, label=f'ClinVar (n={len(clin):,})')

    gm, cm = s['gnom_mean'], s['clin_mean']
    ax.axvline(gm, color=GNOM_C, ls='--', lw=1.5)
    ax.axvline(cm, color=CLIN_C, ls='--', lw=1.5)
    tf = ax.get_xaxis_transform()   # x = data coords, y = axes fraction
    offset = 0.02 * (bins[-1] - bins[0])
    ax.text(gm + offset, 0.6, f'gnomAD {gm:.3f}', color=GNOM_C, transform=tf,
            ha='left', va='top', fontsize=9, fontweight='bold')
    ax.text(cm + offset, 0.7, f'ClinVar {cm:.3f}', color=CLIN_C, transform=tf,
            ha='left', va='top', fontsize=9, fontweight='bold')

    lines = [
        f"{'Test':<13}{'effect':>13}",
        f"{'Mann-Whitney':<13}" + f"AUROC={s['auroc']:.2f}".rjust(13),
        f"{'Welch t':<13}"      + f"d={s['cohend']:+.2f}".rjust(13),
        f"{'KS':<13}"           + f"D={s['ks_d']:.2f}".rjust(13),
    ]
    ax.text(0.03, 0.97, "\n".join(lines), transform=ax.transAxes, ha='left', va='top',
            fontsize=8, family='monospace',
            bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.85))

    ax.set_xlabel(xlabel or score_col); ax.set_ylabel('Density')
    ax.set_title(title or f'{score_col} distribution: ClinVar vs gnomAD')
    ax.legend(frameon=False)
    ax.spines[['top', 'right']].set_visible(False)
    return s

def compare_train_val(train_df, val_df, score_col, xlabel=None, suptitle=None, **filters):
    """Training (in-sample) and held-out distributions side by side on shared bins/axes.
    filters: count_threshold, disorder_range, genes -- applied to both panels."""
    both = np.concatenate([np.concatenate(split_scores(d, score_col, **filters))
                           for d in (train_df, val_df)])
    bins = np.linspace(both.min(), both.max(), 50)

    fig, axes = plt.subplots(1, 2, figsize=(15, 5), sharex=True, sharey=True)
    stats = {}
    for ax, (name, d) in zip(axes, [('Training (in-sample)', train_df),
                                    ('Validation (held-out)', val_df)]):
        stats[name] = plot_score_distribution(d, score_col, title=name, xlabel=xlabel,
                                              ax=ax, bins=bins, **filters)
    axes[1].set_ylabel('')
    fig.suptitle(suptitle or f'{score_col}: ClinVar vs gnomAD')
    fig.tight_layout()
    return fig, stats
