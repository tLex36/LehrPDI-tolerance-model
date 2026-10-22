"""Volcano plots with gene labelling and KEGG pathway colouring, plus an enrichment figure.

    from volcano import volcano, enrichment_plot, load_sheet
    df = load_sheet(path, "P6")
    volcano(df, label_genes=["Nlgn3", "Nrxn1"], color_pathways=["Glutamatergic synapse"])
"""
import re as _re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

UP_COLOR, DOWN_COLOR, NS_COLOR = "#c0392b", "#2c7fb8", "#c8c8c8"
PATHWAY_PALETTE = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                   "#a65628", "#f781bf", "#17becf", "#bcbd22", "#7f7f7f"]



# common MS contaminants: skin/hair keratins, the IP antibody itself, blood and the
# digestion enzyme. Matched case-insensitively against the START of the gene symbol.
CONTAMINANT_PATTERNS = (r"Krt", r"Ig[hkl]", r"Alb$", r"Try", r"Prss1", r"Hb[ab]")


def drop_contaminants(d, patterns=CONTAMINANT_PATTERNS, gene_col="gene", report=False):
    """Remove likely contaminant proteins by gene-symbol prefix.

    Keratins come from skin and hair during handling; Ig genes in an IP are the capture
    antibody, not an interactor. Both routinely carry very strong p-values and will
    dominate a volcano plot's labels if left in.
    """
    pat = r"^(" + "|".join(patterns) + r")"
    m = d[gene_col].astype(str).str.match(pat, case=False, na=False)
    if report:
        print(f"dropped {int(m.sum())} of {len(d)}: " + ", ".join(sorted(d[gene_col][m])[:12])
              + (" ..." if m.sum() > 12 else ""))
    return d[~m].reset_index(drop=True)


# ---------------------------------------------------------------- column discovery
def _find_col(df, patterns, name):
    for pat in patterns:
        for c in df.columns:
            if _re.search(pat, str(c), _re.I):
                return c
    raise KeyError(f"no {name} column found in {list(df.columns)[:12]}")

def gene_col_of(df):  return _find_col(df, [r"^gene symbol$", r"gene.?symbol", r"^gene$", r"^symbol$"], "gene")
def p_col_of(df):     return _find_col(df, [r"p[-_ ]?value"], "p-value")
def ratio_col_of(df): return _find_col(df, [r"abundance ratio(?!.*p-?value)", r"log2.?fold", r"fold.?change", r"^ratio$"], "ratio")


def load_sheet(path, sheet, header=0) -> pd.DataFrame:
    """Read one Excel sheet, skipping Excel lock files."""
    import os
    if os.path.basename(path).startswith("~$"):
        raise ValueError(f"{path} is an Excel lock file, not data")
    return pd.read_excel(path, sheet_name=sheet, header=header)


def prepare(df, gene_col=None, ratio_col=None, p_col=None, log2=True,
            drop_na=True, contaminants=False, require_p=True,
            missing_p=None) -> pd.DataFrame:
    """Normalise any of these sheets to gene / log2fc / pval / neglog10p columns.

    require_p: drop rows with a fold change but no p-value. Proteome Discoverer leaves the
        p blank when a protein is absent from the denominator condition -- in a pulldown
        those are the MOST enriched proteins, so dropping them removes the best hits.
        Set False to keep them (pval stays NaN), or give `missing_p` a value to fill with
        so they can be plotted; 'min' uses the smallest p actually observed.
    """
    g = gene_col or gene_col_of(df)
    r = ratio_col or ratio_col_of(df)
    p = p_col or p_col_of(df)
    out = pd.DataFrame({"gene": df[g].astype(str).str.strip()})
    ratio = pd.to_numeric(df[r], errors="coerce")
    out["log2fc"] = np.log2(ratio.where(ratio > 0)) if log2 else ratio
    out["pval"] = pd.to_numeric(df[p], errors="coerce")
    out["neglog10p"] = -np.log10(out.pval.where(out.pval > 0))
    out = out.replace([np.inf, -np.inf], np.nan)
    out = out[~out.gene.str.lower().isin(["nan", "none", ""])]      # blank Gene Symbol cells
    out["p_missing"] = out.pval.isna()
    if missing_p is not None:
        fill = float(out.pval.min()) if missing_p == "min" else float(missing_p)
        out.loc[out.p_missing, "pval"] = fill
        out["neglog10p"] = -np.log10(out.pval.where(out.pval > 0))
    if drop_na:
        out = out.dropna(subset=["log2fc"])
        if require_p:
            out = out.dropna(subset=["neglog10p"])
    if contaminants is False:
        out = drop_contaminants(out)
    elif contaminants not in (True, None):          # a custom pattern tuple
        out = drop_contaminants(out, patterns=contaminants)
    return out.reset_index(drop=True)


# ---------------------------------------------------------------- gene selection
def _pathway_gene_set(spec, organism):
    from kegg import pathway_members
    return pathway_members(spec, organism=organism)


def select_genes(d, genes=(), pathway=None, top=0, top_by="p", organism="mmu",
                 p_threshold=0.05, fc_threshold=0.0) -> set:
    """Union of: an explicit gene list, a KEGG pathway's members, and the top N hits."""
    chosen = {str(x).strip() for x in genes}
    if pathway:
        specs = [pathway] if isinstance(pathway, str) else list(pathway)
        for s in specs:
            chosen |= _pathway_gene_set(s, organism)
    if top:
        sig = d[(d.pval < p_threshold) & (d.log2fc.abs() >= fc_threshold)]
        key = "neglog10p" if top_by == "p" else "log2fc"
        if top_by == "fc":
            sig = sig.reindex(sig.log2fc.abs().sort_values(ascending=False).index)
        else:
            sig = sig.sort_values(key, ascending=False)
        chosen |= set(sig.gene.head(top))
    return {g for g in chosen if g in set(d.gene)}


# ---------------------------------------------------------------- the plot
def volcano(df, gene_col=None, ratio_col=None, p_col=None, log2=True,
            p_threshold=0.05, fc_threshold=1.0,
            label_genes=(), label_pathway=None, label_top=0, label_top_by="p",
            color_pathways=None, pathway_palette=PATHWAY_PALETTE, organism="mmu",
            up_color=UP_COLOR, down_color=DOWN_COLOR, ns_color=NS_COLOR,
            point_size=12, label_size=8, adjust=True, legend=True,
            title=None, xlabel=None, ylabel=None, figsize=(7, 6), ax=None,
            return_data=False):
    """Volcano plot of log2 fold-change vs -log10 p.

    label_genes:    explicit gene symbols to label
    label_pathway:  KEGG pathway id or name substring -- labels every member present
    label_top:      also label the N most significant hits (label_top_by='p' or 'fc')
    color_pathways: list of pathway ids/names -- members are coloured by pathway instead
                    of by up/down.  A gene in several listed pathways takes the first.
    Returns the axes, or (axes, DataFrame) with return_data=True.
    """
    d = prepare(df, gene_col, ratio_col, p_col, log2=log2)

    sig = (d.pval < p_threshold) & (d.log2fc.abs() >= fc_threshold)
    d["direction"] = np.where(~sig, "ns", np.where(d.log2fc > 0, "up", "down"))

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # ---- colour groups
    d["group"] = d.direction
    pathway_sets = {}
    if color_pathways:
        specs = [color_pathways] if isinstance(color_pathways, str) else list(color_pathways)
        from kegg import pathway_genes
        allp = pathway_genes(organism)
        for s in specs:
            label = allp[s]["name"] if s in allp else s
            pathway_sets[label] = _pathway_gene_set(s, organism)
        assigned = set()
        for label, members in pathway_sets.items():
            m = d.gene.isin(members) & ~d.gene.isin(assigned)
            d.loc[m, "group"] = label
            assigned |= set(d.gene[m])

    base = {"ns": ns_color, "up": up_color, "down": down_color}
    colors = dict(base)
    for i, label in enumerate(pathway_sets):
        colors[label] = pathway_palette[i % len(pathway_palette)]

    # ---- draw: grey first, then significant, then pathway members on top
    order = ["ns", "down", "up"] + list(pathway_sets)
    for g in order:
        sub = d[d.group == g]
        if sub.empty:
            continue
        is_path = g in pathway_sets
        ax.scatter(sub.log2fc, sub.neglog10p, s=point_size * (2.2 if is_path else 1),
                   c=colors[g], alpha=0.9 if g != "ns" else 0.45,
                   edgecolors="black" if is_path else "none",
                   linewidths=0.4 if is_path else 0,
                   zorder=3 if is_path else (2 if g != "ns" else 1),
                   label={"ns": "not significant", "up": "up", "down": "down"}.get(g, g))

    # ---- threshold guides
    ax.axhline(-np.log10(p_threshold), color="#888", lw=0.8, ls="--", zorder=0)
    if fc_threshold:
        for s in (-fc_threshold, fc_threshold):
            ax.axvline(s, color="#888", lw=0.8, ls="--", zorder=0)

    # ---- labels
    chosen = select_genes(d, label_genes, label_pathway, label_top, label_top_by,
                          organism, p_threshold, fc_threshold)
    texts = []
    for _, r in d[d.gene.isin(chosen)].iterrows():
        texts.append(ax.text(r.log2fc, r.neglog10p, r.gene, fontsize=label_size, zorder=5))
    if texts and adjust:
        try:
            from adjustText import adjust_text
            adjust_text(texts, ax=ax,
                        arrowprops=dict(arrowstyle="-", color="#666", lw=0.5))
        except ImportError:
            for t in texts:
                t.set_position((t.get_position()[0] + 0.08, t.get_position()[1] + 0.08))

    ax.set_xlabel(xlabel or "log$_2$ fold change")
    ax.set_ylabel(ylabel or "-log$_{10}$ p-value")
    if title:
        ax.set_title(title, fontsize=11, loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    if legend:
        ax.legend(frameon=False, fontsize=8, loc="best", markerscale=1.4)
    return (ax, d) if return_data else ax


# ---------------------------------------------------------------- enrichment figure
def enrichment_plot(enr, top=15, alpha=0.05, by="q_value", min_show=5, ax=None,
                    figsize=(8, 5), cmap="viridis_r", size_scale=14, title=None,
                    annotate=True):
    """Dot plot of KEGG over-representation results from kegg.enrich().

    x = fold enrichment, dot size = overlapping genes, colour = significance.
    alpha/min_show: pathways passing `alpha` are shown; if fewer than `min_show` pass,
    the top `min_show` by p are shown anyway with a dashed line marking the threshold,
    so a null result reads as a null result instead of an empty panel.
    """
    e = enr.copy()
    if by not in e.columns:
        by = "p_value"
    n_pass = int((e[by] <= alpha).sum())
    e = (e[e[by] <= alpha] if n_pass >= min_show else e.sort_values(by).head(max(min_show, n_pass)))
    e = e.sort_values(by).head(top).iloc[::-1]
    if e.empty:
        raise ValueError("no pathways to plot")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    y = np.arange(len(e))
    sig = e[by] <= alpha
    sc = ax.scatter(e.fold_enrichment, y, s=e.n_overlap * size_scale,
                    c=-np.log10(e[by].clip(lower=1e-300)), cmap=cmap,
                    edgecolors=np.where(sig, "black", "#bbbbbb"),
                    linewidths=np.where(sig, 0.9, 0.5), zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{n}" + ("" if s_ else "  (n.s.)")
                        for n, s_ in zip(e.name, sig)], fontsize=9)
    for tick, s_ in zip(ax.get_yticklabels(), sig):
        tick.set_color("#222" if s_ else "#999")
    ax.set_xlabel("fold enrichment")
    ax.grid(axis="x", color="#e8e8e8", zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.axvline(1, color="#888", lw=0.8, ls="--", zorder=1)
    ax.margins(x=0.18, y=0.08)
    if annotate:
        for yi, (_, r) in zip(y, e.iterrows()):
            ax.text(r.fold_enrichment, yi, f"  {r.n_overlap}/{r.n_pathway}",
                    va="center", fontsize=7.5, color="#444")
    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label(f"-log$_{{10}}$ {by.replace('_', ' ')}", fontsize=9)

    sizes = sorted({int(e.n_overlap.min()), int(e.n_overlap.median()), int(e.n_overlap.max())})
    handles = [plt.scatter([], [], s=n * size_scale, c="#bbb", edgecolors="black",
                           linewidths=0.5) for n in sizes]
    ax.legend(handles, [str(n) for n in sizes], title="genes", frameon=False, fontsize=8,
              title_fontsize=8, loc="upper left", bbox_to_anchor=(1.18, 1.0),
              labelspacing=1.4, borderpad=1.0)
    if title:
        ax.set_title(title, fontsize=11, loc="left",
                     pad=18 if n_pass < min_show else 6)
    if n_pass < min_show:
        ax.text(0.0, 1.015, f"{n_pass} pathway(s) at {by} \u2264 {alpha}; "
                            f"showing top {len(e)} by p",
                transform=ax.transAxes, fontsize=8, color="#999")
    fig.tight_layout()
    return ax


# ---------------------------------------------------------------- two-group contrast
def abundance_columns(df, group) -> list:
    """Per-replicate abundance columns for one sample group.

    Matches Proteome Discoverer's layout, e.g. 'Abundance: F5: Sample, P6, 3'.
    """
    pat = _re.compile(r"^abundance:.*\b" + _re.escape(str(group)) + r"\b", _re.I)
    return [c for c in df.columns if pat.search(str(c))]


def contrast(df, group_a, group_b, gene_col=None, min_valid=2, log2=True,
             equal_var=False, method="fdr_bh") -> pd.DataFrame:
    """Test group_a against group_b from per-replicate abundances.

    Use this instead of subtracting two ratio-vs-control sheets: those p-values each test
    their own ratio against 1, not against each other, and they share the same control
    replicates, so differencing them double-counts the control noise.

    Returns gene / log2fc (a over b) / pval / qval / n_a / n_b, ready to hand to volcano()
    with ratio_col='log2fc', p_col='pval', log2=False.
    """
    from scipy import stats

    g = gene_col or gene_col_of(df)
    ca, cb = abundance_columns(df, group_a), abundance_columns(df, group_b)
    if not ca or not cb:
        raise KeyError(f"found {len(ca)} columns for {group_a!r} and {len(cb)} for {group_b!r}; "
                       f"available: {[c for c in df.columns if str(c).lower().startswith('abundance:')][:6]}")

    A = df[ca].apply(pd.to_numeric, errors="coerce")
    B = df[cb].apply(pd.to_numeric, errors="coerce")
    if log2:                                   # abundances are intensities -> log before testing
        A, B = np.log2(A.where(A > 0)), np.log2(B.where(B > 0))

    keep = (A.notna().sum(axis=1) >= min_valid) & (B.notna().sum(axis=1) >= min_valid)
    A, B = A[keep], B[keep]
    t, p = stats.ttest_ind(A, B, axis=1, equal_var=equal_var, nan_policy="omit")

    out = pd.DataFrame({
        "gene": df.loc[keep, g].astype(str).str.strip().values,
        "log2fc": (A.mean(axis=1) - B.mean(axis=1)).values,
        "pval": np.asarray(p, dtype=float),
        "n_a": A.notna().sum(axis=1).values,
        "n_b": B.notna().sum(axis=1).values,
    }).dropna(subset=["log2fc", "pval"])
    out = out[~out.gene.str.lower().isin(["nan", "none", ""])]

    try:
        from statsmodels.stats.multitest import multipletests
        out["qval"] = multipletests(out.pval, method=method)[1]
    except ImportError:
        o = out.sort_values("pval").reset_index(drop=True)
        out = o.assign(qval=(o.pval * len(o) / (o.index + 1)).cummin()[::-1].cummin()[::-1].clip(upper=1.0))

    out.attrs["group_a"], out.attrs["group_b"] = str(group_a), str(group_b)
    out.attrs["n_rep"] = (len(ca), len(cb))
    return out.reset_index(drop=True)


def ratio_scatter(df_a, df_b, label_a="A", label_b="B", gene_col=None, label_genes=(),
                  label_pathway=None, color_pathways=None, organism="mmu",
                  p_threshold=0.05, label_size=8, point_size=12, figsize=(6.2, 6),
                  ax=None, title=None):
    """Each condition's log2 ratio-vs-control on its own axis.

    The honest summary-statistic comparison when you only have per-sheet ratios: it shows
    where the two conditions agree and disagree without inventing a p-value for the
    difference. Points off the diagonal are the genes that changed between conditions.
    """
    a, b = prepare(df_a, gene_col=gene_col), prepare(df_b, gene_col=gene_col)
    m = a.merge(b, on="gene", suffixes=("_a", "_b"))
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    sig = (m.pval_a < p_threshold) | (m.pval_b < p_threshold)
    ax.scatter(m.log2fc_a[~sig], m.log2fc_b[~sig], s=point_size, c=NS_COLOR, alpha=0.45, zorder=1)
    ax.scatter(m.log2fc_a[sig], m.log2fc_b[sig], s=point_size, c="#555", alpha=0.75, zorder=2)

    if color_pathways:
        from kegg import pathway_genes
        allp = pathway_genes(organism)
        specs = [color_pathways] if isinstance(color_pathways, str) else list(color_pathways)
        assigned = set()
        for i, s in enumerate(specs):
            name = allp[s]["name"] if s in allp else s
            members = _pathway_gene_set(s, organism) - assigned
            sub = m[m.gene.isin(members)]
            assigned |= set(sub.gene)
            ax.scatter(sub.log2fc_a, sub.log2fc_b, s=point_size * 2.2,
                       c=PATHWAY_PALETTE[i % len(PATHWAY_PALETTE)],
                       edgecolors="black", linewidths=0.4, zorder=3, label=name)
        ax.legend(frameon=False, fontsize=8)

    lim = [min(m.log2fc_a.min(), m.log2fc_b.min()), max(m.log2fc_a.max(), m.log2fc_b.max())]
    ax.plot(lim, lim, color="#888", lw=0.8, ls="--", zorder=0)
    ax.axhline(0, color="#ddd", lw=0.8, zorder=0); ax.axvline(0, color="#ddd", lw=0.8, zorder=0)

    chosen = select_genes(m.rename(columns={"log2fc_a": "log2fc", "pval_a": "pval"}),
                          label_genes, label_pathway, 0, "p", organism, p_threshold, 0.0)
    texts = [ax.text(r.log2fc_a, r.log2fc_b, r.gene, fontsize=label_size)
             for _, r in m[m.gene.isin(chosen)].iterrows()]
    if texts:
        try:
            from adjustText import adjust_text
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#666", lw=0.5))
        except ImportError:
            pass

    ax.set_xlabel(f"log$_2$ {label_a} / control")
    ax.set_ylabel(f"log$_2$ {label_b} / control")
    ax.set_aspect("equal", adjustable="datalim")
    ax.spines[["top", "right"]].set_visible(False)
    if title:
        ax.set_title(title, fontsize=11, loc="left")
    return ax


# ---------------------------------------------------------------- pathway direction
def pathway_activity(d, organism="mmu", min_genes=5, max_genes=500,
                     method="fdr_bh") -> pd.DataFrame:
    """Per-pathway summary of fold-change direction for one dataset.

    d: a prepared frame (gene / log2fc / pval), e.g. from prepare() or contrast().
    Each pathway's member genes are compared against all other detected genes with a
    Mann-Whitney U test, so a pathway is "up" when its genes are shifted higher as a
    group -- not merely when a few members happen to be significant.
    """
    from scipy.stats import mannwhitneyu
    from kegg import pathway_genes

    vals = d.set_index("gene").log2fc.groupby(level=0).mean()     # one value per gene
    paths = pathway_genes(organism, min_size=min_genes, max_size=max_genes)
    rows = []
    for pid, info in paths.items():
        members = [g for g in info["genes"] if g in vals.index]
        if len(members) < min_genes:
            continue
        inside = vals.loc[members]
        outside = vals.drop(index=members, errors="ignore")
        if len(outside) < min_genes:
            continue
        try:
            u, p = mannwhitneyu(inside, outside, alternative="two-sided")
        except ValueError:
            continue
        rows.append({"pathway": pid, "name": info["name"], "n_genes": len(members),
                     "median_log2fc": float(inside.median()),
                     "mean_log2fc": float(inside.mean()),
                     "median_shift": float(inside.median() - outside.median()),
                     "p_value": float(p), "genes": ",".join(sorted(members))})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    try:
        from statsmodels.stats.multitest import multipletests
        out["q_value"] = multipletests(out.p_value, method=method)[1]
    except ImportError:
        o = out.sort_values("p_value").reset_index(drop=True)
        out = o.assign(q_value=(o.p_value * len(o) / (o.index + 1)).cummin()[::-1].cummin()[::-1].clip(upper=1.0))
    out.attrs["background_median"] = float(vals.median())
    return out.sort_values("median_shift", ascending=False).reset_index(drop=True)


def pathway_activity_plot(d, act=None, organism="mmu", top=15, direction="up",
                          alpha=0.05, min_genes=5, show_points=True, ax=None,
                          figsize=(9, 6), title=None, cmap="RdBu_r", vlim=None):
    """Boxplot of member-gene log2 fold-change per KEGG pathway, ranked by shift.

    direction: 'up' (most positively shifted), 'down', or 'both' (top/2 of each).
    Boxes are coloured by their median shift; pathways passing `alpha` are labelled in
    black with their q, the rest greyed, so an unconvincing result stays visibly so.
    """
    from matplotlib import cm, colors as mcolors

    act = pathway_activity(d, organism=organism, min_genes=min_genes) if act is None else act
    if act.empty:
        raise ValueError("no pathways with enough detected genes")

    if direction == "up":
        sel = act.sort_values("median_shift", ascending=False).head(top)
    elif direction == "down":
        sel = act.sort_values("median_shift").head(top)
    else:
        k = max(1, top // 2)
        sel = pd.concat([act.sort_values("median_shift", ascending=False).head(k),
                         act.sort_values("median_shift").head(k)]).drop_duplicates("pathway")
    sel = sel.sort_values("median_shift")

    vals = d.set_index("gene").log2fc.groupby(level=0).mean()
    data = [vals.loc[[g for g in row.genes.split(",") if g in vals.index]].values
            for row in sel.itertuples()]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    lim = vlim or float(sel["median_shift"].abs().max()) or 1.0
    norm = mcolors.TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=lim)
    cmo = cm.get_cmap(cmap) if hasattr(cm, "get_cmap") else plt.get_cmap(cmap)

    y = np.arange(len(sel)) + 1
    bp = ax.boxplot(data, vert=False, widths=0.62, patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", lw=1.4),
                    whiskerprops=dict(color="#888", lw=0.9),
                    capprops=dict(color="#888", lw=0.9))
    for patch, shift in zip(bp["boxes"], sel["median_shift"]):
        patch.set_facecolor(cmo(norm(shift)))
        patch.set_edgecolor("#444")
        patch.set_linewidth(0.7)
        patch.set_alpha(0.95)
    if show_points:
        for yi, vv in zip(y, data):
            ax.scatter(vv, np.full(len(vv), yi) + np.random.uniform(-0.16, 0.16, len(vv)),
                       s=6, color="#333", alpha=0.35, zorder=3, linewidths=0)

    sig = (sel.q_value <= alpha).values
    ax.set_yticks(y)
    ax.set_yticklabels([f"{n}  ({int(k)})" for n, k in zip(sel.name, sel.n_genes)], fontsize=9)
    for tick, s_ in zip(ax.get_yticklabels(), sig):
        tick.set_color("#222" if s_ else "#999")

    ax.axvline(0, color="#888", lw=0.9, ls="--", zorder=0)
    bgm = d.log2fc.median()
    ax.axvline(bgm, color="#2ca02c", lw=0.9, ls=":", zorder=0,
               label=f"all genes (median {bgm:.2f})")
    x0 = ax.get_xlim()[1]
    for yi, q, s_ in zip(y, sel.q_value, sig):
        ax.text(x0, yi, f"  q={q:.1e}" if s_ else "  n.s.", va="center", fontsize=7.5,
                color="#222" if s_ else "#aaa")
    ax.margins(x=0.14)
    ax.set_xlabel("log$_2$ fold change of pathway member genes")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="upper center",
              bbox_to_anchor=(0.5, -0.10), ncol=2)
    sm = cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, pad=0.10)
    cb.set_label("median shift vs all genes", fontsize=9)
    if title:
        ax.set_title(title, fontsize=11, loc="left")
    fig.tight_layout()
    return ax


def enrich_from(d, p_threshold=0.05, fc_threshold=1.0, direction="up", organism="mmu",
                p_col="pval", **kw):
    """KEGG over-representation for one prepared frame, with the background set correctly.

    d:          a prepared frame (prepare() or contrast())
    direction:  'up', 'down' or 'both' -- which side of the volcano forms the query
    The background is every gene in `d`, i.e. everything detected and quantified, which is
    the right universe; passing no background at all would test against the whole genome
    and inflate every p-value.
    """
    from kegg import enrich
    p = d[p_col] if p_col in d.columns else d["pval"]
    sel = p < p_threshold
    if direction == "up":
        sel &= d.log2fc >= fc_threshold
    elif direction == "down":
        sel &= d.log2fc <= -fc_threshold
    else:
        sel &= d.log2fc.abs() >= fc_threshold
    out = enrich(d.gene[sel].tolist(), background=d.gene.tolist(), organism=organism, **kw)
    out.attrs["direction"] = direction
    out.attrs["n_selected"] = int(sel.sum())
    out.attrs["n_total"] = len(d)
    return out


# ---------------------------------------------------------------- timepoint series
def ratio_series(df, denom="Ctrl", order=None, strip=r"[-_ ]?(NL3|NLGN3|HA)$",
                 **prep) -> dict:
    """{sample: prepared frame} for every 'Abundance Ratio: (X) / (denom)' column.

    One Proteome Discoverer sheet usually holds several conditions scored against a shared
    control; this pulls each out as its own gene/log2fc/pval frame so they can be compared
    as a series.

    denom: the control's name, or None to accept any denominator -- needed when each
           sample has its own control, as in '(CTX-NL3) / (CTX-Ctrl)'.
    strip: regex removed from the numerator to name the series, so 'CTX-NL3' -> 'CTX'.
    Extra keyword arguments go to prepare() -- require_p=False keeps proteins that have a
    ratio but no p-value (absent from the denominator condition).
    """
    pat = _re.compile(r"abundance ratio:\s*\((?P<num>[^)]+)\)\s*/\s*\((?P<den>[^)]+)\)", _re.I)
    out = {}
    for c in df.columns:
        m = pat.match(str(c).strip())
        if not m:
            continue
        den = m.group("den").strip()
        if denom is not None and den.lower() != str(denom).lower():
            continue
        num = m.group("num").strip()
        p = next((pc for pc in df.columns
                  if "p-value" in str(pc).lower() and f"({num})" in str(pc)
                  and f"({den})" in str(pc)), None)
        if p is None:
            continue
        label = _re.sub(strip, "", num, flags=_re.I).strip() or num
        out[label] = prepare(df, ratio_col=c, p_col=p, **prep)
    if order:
        out = {k: out[k] for k in order if k in out}
    return out


def _numeric_x(labels):
    """'P6' -> 6, falling back to evenly spaced positions when labels aren't numeric."""
    xs = []
    for s in labels:
        m = _re.search(r"(\d+(?:\.\d+)?)", str(s))
        xs.append(float(m.group(1)) if m else None)
    return xs if all(x is not None for x in xs) else list(range(len(labels)))


MODULE_COLORS = ["#d62728", "#ff7f0e", "#1f77b4", "#9467bd", "#8c564b", "#2ca02c",
                 "#e377c2", "#17becf"]


def module_trajectory(series, modules, ax=None, figsize=(7.5, 4.6), x=None, iqr=True,
                      show_all=True, colors=None, span=None, span_label=None,
                      lw=2.2, ms=5, legend=True, title=None, ylabel=None, xlabel=None,
                      return_data=False):
    """Median enrichment of each gene module across a series of samples.

    series:  {sample label: prepared frame}, e.g. from ratio_series()
    modules: {module name: [gene symbols]}
    span:    (lo, hi) in x units to shade -- e.g. (10, 14) for myelination onset
    Returns the axes, or (axes, DataFrame of medians) with return_data=True.
    """
    labels = list(series)
    xs = x if x is not None else _numeric_x(labels)
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    cols = colors or {m: MODULE_COLORS[i % len(MODULE_COLORS)]
                      for i, m in enumerate(modules)}

    if span:
        ax.axvspan(span[0], span[1], color="#ffe9a8", alpha=0.5, zorder=0)
        if span_label:
            ax.text(sum(span) / 2, ax.get_ylim()[0], span_label, ha="center", va="bottom",
                    fontsize=8, color="#8a6d00", zorder=1)
    if show_all:
        ax.plot(xs, [float(series[t].log2fc.median()) for t in labels], color="0.5",
                lw=1.4, ls=":", marker="o", ms=4, label="all proteins", zorder=1)

    rows = []
    for mod, genes in modules.items():
        med, lo, hi, n = [], [], [], []
        for t in labels:
            v = series[t].set_index("gene").log2fc
            g = [x_ for x_ in genes if x_ in v.index]
            n.append(len(g))
            med.append(float(v[g].median()) if g else np.nan)
            lo.append(float(v[g].quantile(0.25)) if g else np.nan)
            hi.append(float(v[g].quantile(0.75)) if g else np.nan)
        c = cols[mod] if isinstance(cols, dict) else cols[list(modules).index(mod)]
        if iqr:
            ax.fill_between(xs, lo, hi, color=c, alpha=0.11, lw=0, zorder=2)
        ax.plot(xs, med, color=c, lw=lw, marker="o", ms=ms, zorder=3,
                label=f"{mod} (n={max(n)})")
        rows.append(dict(module=mod, n=max(n), **dict(zip(labels, [round(v, 2) for v in med])),
                         delta=round(med[-1] - med[0], 2)))

    ax.axhline(0, color="black", lw=0.8, zorder=1)
    ax.set_xticks(xs); ax.set_xticklabels(labels)
    ax.set_xlabel(xlabel or "sample")
    ax.set_ylabel(ylabel or "median log$_2$ enrichment over control")
    if title:
        ax.set_title(title, fontsize=11, loc="left")
    if legend:
        ax.legend(frameon=False, fontsize=8.5, loc="best")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    tab = pd.DataFrame(rows)
    return (ax, tab) if return_data else ax


def gene_trajectory(series, genes, ax=None, figsize=(6.5, 4.6), x=None, colors=None,
                    label_last=True, lw=1.2, alpha=0.6, ms=2.8, span=None,
                    title=None, ylabel=None, xlabel=None, label_size=8.5, pad=0.14):
    """One line per gene across the series.

    genes: a flat list, or {group name: [genes]} to colour by group (colours from
           `colors`, a {group: colour} dict).
    """
    labels = list(series)
    xs = x if x is not None else _numeric_x(labels)
    groups = genes if isinstance(genes, dict) else {"": list(genes)}
    cols = colors or {g: MODULE_COLORS[i % len(MODULE_COLORS)]
                      for i, g in enumerate(groups)}
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    if span:
        ax.axvspan(span[0], span[1], color="#ffe9a8", alpha=0.5, zorder=0)

    ends = []
    for grp, gl in groups.items():
        c = cols[grp] if isinstance(cols, dict) else cols
        for g in gl:
            y = [float(series[t].set_index("gene").log2fc[g])
                 if g in set(series[t].gene) else np.nan for t in labels]
            if np.isnan(y).all():
                continue
            ax.plot(xs, y, color=c, lw=lw, alpha=alpha, marker="o", ms=ms, zorder=2)
            ends.append((y[-1], g, c))
    if label_last:                       # stack end labels so they do not overlap
        ends.sort(reverse=True)
        span_y = (max(e[0] for e in ends) - min(e[0] for e in ends)) or 1
        gap, last = span_y * pad, None
        for y, g, c in ends:
            yy = y if last is None else min(y, last - gap)
            ax.text(max(xs) + (max(xs) - min(xs)) * 0.02, yy, g, fontsize=label_size,
                    va="center", color=c)
            last = yy
        ax.set_xlim(min(xs) - (max(xs) - min(xs)) * 0.03,
                    max(xs) + (max(xs) - min(xs)) * 0.18)
    ax.axhline(0, color="black", lw=0.8, zorder=1)
    ax.set_xticks(xs); ax.set_xticklabels(labels)
    ax.set_xlabel(xlabel or "sample")
    ax.set_ylabel(ylabel or "log$_2$ enrichment over control")
    if title:
        ax.set_title(title, fontsize=11, loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return ax
