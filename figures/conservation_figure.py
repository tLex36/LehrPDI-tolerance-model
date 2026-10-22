"""Tree + pairwise-identity heatmap figure.

    tree = build_tree(keep=set(mat.index)); tree.ladderize(reverse=True)
    fig = tree_heatmap_figure(tree, res["identity_matrix"], title="NLGN1 ICD")
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from species_tree import draw_species_tree, label_nodes

# sequential single-hue ramp (light -> dark)
BLUE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
SEQ_BLUE = LinearSegmentedColormap.from_list("seq_blue", BLUE_RAMP)
# diverging red -> neutral -> blue (low = red, high = blue)
RED_BLUE = LinearSegmentedColormap.from_list("red_blue", ["#b2182b", "#ef8a62", "#f7f7f7", "#67a9cf", "#2166ac"])

def _ink(rgba):
    """Black or white text depending on cell luminance."""
    r, g, b = rgba[:3]
    return "#1a1a1a" if 0.299 * r + 0.587 * g + 0.114 * b > 0.55 else "white"

def plot_identity_heatmap(mat, order=None, ax=None, cmap=RED_BLUE, vmin=None, vmax=None,
                          triangle="lower", diagonal=True, row_labels=True, col_labels=True,
                          cbar=True, cbar_label="% identity", aspect="equal", y0=0, fontsize=8):
    """mat: square DataFrame of % identity.
    order:    row/col order (e.g. tree tip order).
    vmin/vmax: colour range; default = min/max of the off-diagonal values.
    triangle: "lower", "upper" or None (full matrix).
    y0:       y-coordinate of the first row centre (Phylo.draw puts tips at 1..n; pass 1 to align)."""
    if order is not None:
        order = [s for s in order if s in mat.index]
        mat = mat.loc[order, order]
    vals = mat.to_numpy(dtype=float).copy()
    n = len(mat)
    off = ~np.eye(n, dtype=bool)
    if vmin is None: vmin = np.nanmin(vals[off])
    if vmax is None: vmax = np.nanmax(vals[off])
    # mask the unwanted triangle / diagonal
    if triangle == "lower":   vals[np.triu_indices(n, 1)] = np.nan
    elif triangle == "upper": vals[np.tril_indices(n, -1)] = np.nan
    if not diagonal:          vals[np.diag_indices(n)] = np.nan

    if ax is None:
        fig, ax = plt.subplots(figsize=(0.5 * n + 1.5, 0.5 * n + 1))
    cmap = plt.get_cmap(cmap).copy(); cmap.set_bad(alpha=0)
    ys = np.arange(y0, y0 + n)
    im = ax.imshow(vals, cmap=cmap, vmin=vmin, vmax=vmax, aspect=aspect,
                   extent=(-0.5, n - 0.5, y0 + n - 0.5, y0 - 0.5))
    # 2px surface gap between cells
    ax.set_xticks(np.arange(-0.5, n, 1), minor=True); ax.set_yticks(ys - 0.5, minor=True)
    ax.grid(which="minor", color="white", linewidth=2); ax.tick_params(which="minor", length=0)
    # integer % in every drawn cell
    norm = (vals - vmin) / (vmax - vmin) if vmax > vmin else np.zeros_like(vals)
    for i in range(n):
        for j in range(n):
            if np.isnan(vals[i, j]): continue
            ax.text(j, ys[i], f"{vals[i, j]:.0f}", ha="center", va="center",
                    fontsize=fontsize, color=_ink(cmap(norm[i, j])))
    ax.set_xticks(range(n)); ax.set_xticklabels(mat.columns if col_labels else [], rotation=90, fontsize=9)
    ax.set_yticks(ys);       ax.set_yticklabels(mat.index if row_labels else [], fontsize=9)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(length=0)
    if cbar:
        cb = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        cb.set_label(cbar_label); cb.outline.set_visible(False)
    return ax

def tree_heatmap_figure(tree, mat, title=None, figsize=None, cell=0.45, labels_on="tree",
                        clade_labels=None, **kw):
    """Species tree on the left, identity heatmap on the right, rows aligned to tree tips.
    labels_on: "tree" (names next to tips) or "heatmap" (names on the heatmap y-axis).
    Extra keyword args go to plot_identity_heatmap (cmap, triangle, vmin, ...)."""
    order = [t.name for t in tree.get_terminals()]
    n = len(order)
    figsize = figsize or (cell * n + 5, cell * n + 1.5)
    fig, (ax_tree, ax_hm) = plt.subplots(1, 2, figsize=figsize,
                                         gridspec_kw={"width_ratios": [2.2, cell * n / 1.6 + 1],
                                                      "wspace": 0.05 if labels_on == "tree" else 0.45})
    draw_species_tree(tree, ax=ax_tree, title="", labels=(labels_on == "tree"))
    if clade_labels:                                   # list of clade names, or True for all named
        label_nodes(tree, ax_tree, names=None if clade_labels is True else clade_labels)
    plot_identity_heatmap(mat, order=order, ax=ax_hm, aspect="auto", y0=1,
                          row_labels=(labels_on == "heatmap"), **kw)
    for ax in (ax_tree, ax_hm):                        # identical y-coords: tips and rows at 1..n
        ax.set_ylim(n + 0.5, 0.5)
    if labels_on == "tree":                            # leave room for tip names before the heatmap
        ax_tree.set_xlim(ax_tree.get_xlim()[0], ax_tree.get_xlim()[1] * 1.6)
    if title: fig.suptitle(title, y=0.99)
    return fig

def plot_matrix_heatmap(df, ax=None, cmap="BuPu", vmin=None, vmax=None, row_labels=True,
                        cbar=True, cbar_label="% identity", aspect="equal", y0=0, fontsize=8):
    """Rectangular heatmap (rows x cols DataFrame), integer labels in cells, NaN cells blank."""
    vals = df.to_numpy(dtype=float)
    n, m = vals.shape
    if vmin is None: vmin = np.nanmin(vals)
    if vmax is None: vmax = np.nanmax(vals)
    if ax is None:
        fig, ax = plt.subplots(figsize=(0.6 * m + 2, 0.45 * n + 1))
    cmap = plt.get_cmap(cmap).copy(); cmap.set_bad(alpha=0)
    ys = np.arange(y0, y0 + n)
    im = ax.imshow(vals, cmap=cmap, vmin=vmin, vmax=vmax, aspect=aspect,
                   extent=(-0.5, m - 0.5, y0 + n - 0.5, y0 - 0.5))
    ax.set_xticks(np.arange(-0.5, m, 1), minor=True); ax.set_yticks(ys - 0.5, minor=True)
    ax.grid(which="minor", color="white", linewidth=2); ax.tick_params(which="minor", length=0)
    norm = (vals - vmin) / (vmax - vmin) if vmax > vmin else np.zeros_like(vals)
    for i in range(n):
        for j in range(m):
            if np.isnan(vals[i, j]): continue
            ax.text(j, ys[i], f"{vals[i, j]:.0f}", ha="center", va="center",
                    fontsize=fontsize, color=_ink(cmap(norm[i, j])))
    ax.set_xticks(range(m)); ax.set_xticklabels(df.columns, fontsize=9)
    ax.set_yticks(ys);       ax.set_yticklabels(df.index if row_labels else [], fontsize=9)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(length=0)
    if cbar:
        cb = plt.colorbar(im, ax=ax, fraction=0.06, pad=0.04)
        cb.set_label(cbar_label); cb.outline.set_visible(False)
    return ax

def tree_gene_figure(tree, df, title=None, figsize=None, cell=0.45, labels_on="tree",
                     clade_labels=None, **kw):
    """Species tree on the left, species x gene heatmap on the right, rows aligned to tips.
    df: rows=species, cols=genes (e.g. from conservation.species_vs_reference)."""
    order = [t.name for t in tree.get_terminals()]
    df = df.reindex(order)                              # tree order; species missing from df -> blank row
    n, m = df.shape
    figsize = figsize or (0.6 * m + 6, cell * n + 1.5)
    fig, (ax_tree, ax_hm) = plt.subplots(1, 2, figsize=figsize,
                                         gridspec_kw={"width_ratios": [2.2, 0.6 * m + 0.8],
                                                      "wspace": 0.05 if labels_on == "tree" else 0.45})
    draw_species_tree(tree, ax=ax_tree, title="", labels=(labels_on == "tree"))
    if clade_labels:                                   # list of clade names, or True for all named
        label_nodes(tree, ax_tree, names=None if clade_labels is True else clade_labels)
    plot_matrix_heatmap(df, ax=ax_hm, aspect="auto", y0=1, row_labels=(labels_on == "heatmap"), **kw)
    for ax in (ax_tree, ax_hm):
        ax.set_ylim(n + 0.5, 0.5)
    if labels_on == "tree":
        ax_tree.set_xlim(ax_tree.get_xlim()[0], ax_tree.get_xlim()[1] * 1.6)
    if title: fig.suptitle(title, y=0.99)
    return fig

# ---------------------------------------------------------------- two domains in one heatmap
from matplotlib.patches import Polygon

def _prep(dfs, order=None):
    """dfs: {domain: DataFrame rows=species, cols=genes}. Align all to common rows/cols."""
    names = list(dfs)
    genes = list(dfs[names[0]].columns)
    rows = order if order is not None else list(dfs[names[0]].index)
    dfs = {k: v.reindex(index=rows, columns=genes) for k, v in dfs.items()}
    return names, genes, rows, dfs

def _range(dfs, vmin, vmax):
    allv = np.concatenate([d.to_numpy(dtype=float).ravel() for d in dfs.values()])
    return (np.nanmin(allv) if vmin is None else vmin, np.nanmax(allv) if vmax is None else vmax)

def plot_grouped_heatmap(dfs, order=None, ax=None, cmap="BuPu", vmin=None, vmax=None,
                         row_labels=True, cbar=True, cbar_label="% identity", aspect="equal",
                         y0=0, fontsize=8, group_gap=0.35):
    """Option A: one column per (gene, domain), grouped by gene with the gene name above.
    dfs: {"ECD": df, "ICD": df} — column order within a group follows dict order."""
    names, genes, rows, dfs = _prep(dfs, order)
    vmin, vmax = _range(dfs, vmin, vmax)
    n, k, m = len(rows), len(names), len(genes)
    if ax is None:
        fig, ax = plt.subplots(figsize=(0.6 * k * m + 2, 0.45 * n + 1.5))
    cmap = plt.get_cmap(cmap).copy(); cmap.set_bad(alpha=0)
    ys = np.arange(y0, y0 + n)
    # x-position of each column, with an extra gap between gene groups
    xs = [g * (k + group_gap) + d for g in range(m) for d in range(k)]
    norm = lambda v: (v - vmin) / (vmax - vmin) if vmax > vmin else 0
    for gi, gene in enumerate(genes):
        for di, dom in enumerate(names):
            x = xs[gi * k + di]
            for ri, sp in enumerate(rows):
                v = dfs[dom].loc[sp, gene]
                if np.isnan(v): continue
                c = cmap(norm(v))
                ax.add_patch(plt.Rectangle((x - 0.5, ys[ri] - 0.5), 1, 1, facecolor=c, edgecolor="white", lw=2))
                ax.text(x, ys[ri], f"{v:.0f}", ha="center", va="center", fontsize=fontsize, color=_ink(c))
        gx = np.mean(xs[gi * k:(gi + 1) * k])
        ax.text(gx, y0 - 0.9, gene, ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_xlim(-0.6, xs[-1] + 0.6); ax.set_ylim(y0 + n - 0.5, y0 - 1.6)
    ax.set_aspect(aspect)
    ax.set_xticks(xs); ax.set_xticklabels(names * m, fontsize=8)
    ax.set_yticks(ys); ax.set_yticklabels(rows if row_labels else [], fontsize=9)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(length=0)
    if cbar:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin, vmax))
        cb = plt.colorbar(sm, ax=ax, fraction=0.05, pad=0.04)
        cb.set_label(cbar_label); cb.outline.set_visible(False)
    return ax

def plot_split_heatmap(dfs, order=None, ax=None, cmap="BuPu", vmin=None, vmax=None,
                       row_labels=True, cbar=True, cbar_label="% identity", aspect="equal",
                       y0=0, fontsize=7):
    """Option B: each cell split diagonally — first domain in dfs lower-left, second upper-right."""
    names, genes, rows, dfs = _prep(dfs, order)
    assert len(names) == 2, "split cells take exactly two domains"
    lo, hi = names
    vmin, vmax = _range(dfs, vmin, vmax)
    n, m = len(rows), len(genes)
    if ax is None:
        fig, ax = plt.subplots(figsize=(0.9 * m + 2, 0.55 * n + 1.5))
    cmap = plt.get_cmap(cmap).copy()
    ys = np.arange(y0, y0 + n)
    norm = lambda v: (v - vmin) / (vmax - vmin) if vmax > vmin else 0
    for gi, gene in enumerate(genes):
        for ri, sp in enumerate(rows):
            x, y = gi, ys[ri]
            L, R, T, B = x - 0.5, x + 0.5, y - 0.5, y + 0.5      # y axis is inverted: T < B
            for dom, tri, tx, ty in ((lo, [(L, T), (L, B), (R, B)], x - 0.18, y + 0.2),
                                     (hi, [(L, T), (R, T), (R, B)], x + 0.18, y - 0.2)):
                v = dfs[dom].loc[sp, gene]
                if np.isnan(v): continue
                c = cmap(norm(v))
                ax.add_patch(Polygon(tri, closed=True, facecolor=c, edgecolor="white", lw=1.5))
                ax.text(tx, ty, f"{v:.0f}", ha="center", va="center", fontsize=fontsize, color=_ink(c))
    ax.set_xlim(-0.6, m - 0.4); ax.set_ylim(y0 + n - 0.5, y0 - 0.5)
    ax.set_aspect(aspect)
    ax.set_xticks(range(m)); ax.set_xticklabels(genes, fontsize=10)
    ax.set_yticks(ys); ax.set_yticklabels(rows if row_labels else [], fontsize=9)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(length=0)
    # legend key: one split cell showing which triangle is which
    kx, ky = m - 0.5 + 0.9, y0 - 0.5 + 0.5
    ax.add_patch(Polygon([(kx-0.4, ky-0.4), (kx-0.4, ky+0.4), (kx+0.4, ky+0.4)], facecolor="#ddd", edgecolor="white", clip_on=False))
    ax.add_patch(Polygon([(kx-0.4, ky-0.4), (kx+0.4, ky-0.4), (kx+0.4, ky+0.4)], facecolor="#999", edgecolor="white", clip_on=False))
    ax.text(kx - 0.15, ky + 0.17, lo, fontsize=6, ha="center", va="center", clip_on=False)
    ax.text(kx + 0.15, ky - 0.17, hi, fontsize=6, ha="center", va="center", color="white", clip_on=False)
    if cbar:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin, vmax))
        cb = plt.colorbar(sm, ax=ax, fraction=0.05, pad=0.12)
        cb.set_label(cbar_label); cb.outline.set_visible(False)
    return ax

def tree_domains_figure(tree, dfs, style="grouped", title=None, figsize=None, cell=0.45,
                        labels_on="tree", clade_labels=None, **kw):
    """Species tree + combined-domain heatmap. dfs: {"ECD": df, "ICD": df}.
    style: "grouped" (two columns per gene) or "split" (diagonal cells)."""
    order = [t.name for t in tree.get_terminals()]
    n = len(order)
    names, genes, _, _ = _prep(dfs, order)
    width = (0.6 * len(names) * len(genes) + 1) if style == "grouped" else (0.9 * len(genes) + 1.5)
    figsize = figsize or (width + 5.5, cell * n + 1.8)
    fig, (ax_tree, ax_hm) = plt.subplots(1, 2, figsize=figsize,
                                         gridspec_kw={"width_ratios": [2.2, width],
                                                      "wspace": 0.05 if labels_on == "tree" else 0.45})
    draw_species_tree(tree, ax=ax_tree, title="", labels=(labels_on == "tree"))
    if clade_labels:
        label_nodes(tree, ax_tree, names=None if clade_labels is True else clade_labels)
    fn = plot_grouped_heatmap if style == "grouped" else plot_split_heatmap
    fn(dfs, order=order, ax=ax_hm, aspect="auto", y0=1, row_labels=(labels_on == "heatmap"), **kw)
    lo = 1 - 1.6 if style == "grouped" else 0.5
    for ax in (ax_tree, ax_hm):
        ax.set_ylim(n + 0.5, lo)
    if labels_on == "tree":
        ax_tree.set_xlim(ax_tree.get_xlim()[0], ax_tree.get_xlim()[1] * 1.6)
    if title: fig.suptitle(title, y=0.99)
    return fig

def plot_two_triangle_heatmap(lower, upper, lower_label, upper_label, ax=None, cmap="BuPu",
                              vmin=None, vmax=None, cbar_label="% identity", label_fontsize=10, **kw):
    """One square heatmap from two symmetric matrices: `lower` fills the bottom-left triangle,
    `upper` the top-right; diagonal left blank. Colour scale is shared."""
    names = list(lower.index)
    upper = upper.reindex(index=names, columns=names)
    combo = lower.copy().astype(float)
    n = len(names)
    iu = np.triu_indices(n, 1); il = np.tril_indices(n, -1)
    vals = combo.to_numpy(); uv = upper.to_numpy(dtype=float)
    vals[iu] = uv[iu]; vals[np.diag_indices(n)] = np.nan
    combo = pd.DataFrame(vals, index=names, columns=names)
    ax = plot_identity_heatmap(combo, ax=ax, cmap=cmap, vmin=vmin, vmax=vmax, triangle=None,
                               diagonal=False, cbar_label=cbar_label, **kw)
    # labels in the blank diagonal: top-left cell heads the `upper` rows, bottom-right ends the `lower` rows
    y0 = kw.get("y0", 0)
    ax.text(0, y0, upper_label, ha="center", va="center", fontsize=label_fontsize, fontweight="bold", color="#444")
    ax.text(n - 1, y0 + n - 1, lower_label, ha="center", va="center", fontsize=label_fontsize, fontweight="bold", color="#444")
    return ax

def plot_combo_heatmap(full, ecd, icd, secondary=None, ax=None, cmap="BuPu", vmin=None, vmax=None,
                       cbar_label="% identity", fontsize=9, split_fontsize=7, key=True):
    """5x5-style paralog chart. Lower-left triangle: `full` (one number per cell; if `secondary`
    is given, its value is printed smaller underneath, e.g. similarity under identity).
    Upper-right triangle: each cell split diagonally, `ecd` lower-left / `icd` upper-right.
    All three share one colour scale."""
    names = list(full.index)
    n = len(names)
    ecd = ecd.reindex(index=names, columns=names); icd = icd.reindex(index=names, columns=names)
    allv = np.concatenate([m.to_numpy(dtype=float)[~np.eye(n, dtype=bool)] for m in (full, ecd, icd)])
    if vmin is None: vmin = np.nanmin(allv)
    if vmax is None: vmax = np.nanmax(allv)
    cmap = plt.get_cmap(cmap)
    norm = lambda v: (v - vmin) / (vmax - vmin) if vmax > vmin else 0
    if ax is None:
        fig, ax = plt.subplots(figsize=(0.75 * n + 2, 0.75 * n + 1))
    for i in range(n):
        for j in range(n):
            x, y = j, i
            L, R, T, B = x - 0.5, x + 0.5, y - 0.5, y + 0.5
            if j < i:                                              # lower-left: full-length
                v = float(full.iloc[i, j])
                if np.isnan(v): continue
                c = cmap(norm(v))
                ax.add_patch(plt.Rectangle((L, T), 1, 1, facecolor=c, edgecolor="white", lw=2))
                if secondary is not None:
                    v2 = float(secondary.iloc[i, j])
                    ax.text(x, y - 0.12, f"{v:.0f}", ha="center", va="center", fontsize=fontsize, color=_ink(c), fontweight="bold")
                    ax.text(x, y + 0.22, f"{v2:.0f}", ha="center", va="center", fontsize=fontsize - 2, color=_ink(c))
                else:
                    ax.text(x, y, f"{v:.0f}", ha="center", va="center", fontsize=fontsize, color=_ink(c))
            elif j > i:                                            # upper-right: split ECD / ICD
                for m, tri, tx, ty in ((ecd, [(L, T), (L, B), (R, B)], x - 0.18, y + 0.2),
                                       (icd, [(L, T), (R, T), (R, B)], x + 0.18, y - 0.2)):
                    v = float(m.iloc[i, j])
                    if np.isnan(v): continue
                    c = cmap(norm(v))
                    ax.add_patch(Polygon(tri, closed=True, facecolor=c, edgecolor="white", lw=1.5))
                    ax.text(tx, ty, f"{v:.0f}", ha="center", va="center", fontsize=split_fontsize, color=_ink(c))
    ax.set_xlim(-0.6, n - 0.4); ax.set_ylim(n - 0.4, -0.6); ax.set_aspect("equal")
    ax.set_xticks(range(n)); ax.set_xticklabels(names, rotation=90, fontsize=9)
    ax.set_yticks(range(n)); ax.set_yticklabels(names, fontsize=9)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(length=0)
    if key:                                                        # key drawn in the diagonal cells
        ax.text(0, 0, "full", ha="center", va="center", fontsize=8, color="#444", fontweight="bold")
        kx, ky = n - 1, n - 1
        ax.add_patch(Polygon([(kx-0.5, ky-0.5), (kx-0.5, ky+0.5), (kx+0.5, ky+0.5)], facecolor="#ddd", edgecolor="white"))
        ax.add_patch(Polygon([(kx-0.5, ky-0.5), (kx+0.5, ky-0.5), (kx+0.5, ky+0.5)], facecolor="#aaa", edgecolor="white"))
        ax.text(kx - 0.18, ky + 0.2, "ECD", fontsize=7, ha="center", va="center", color="#333")
        ax.text(kx + 0.18, ky - 0.2, "ICD", fontsize=7, ha="center", va="center", color="white")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin, vmax))
    cb = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(cbar_label); cb.outline.set_visible(False)
    return ax
