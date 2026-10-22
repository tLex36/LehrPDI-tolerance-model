"""Poster figures as editable SVG (text kept as text, Arial), plus PNG previews.

    cd figures/notebooks && python ../results/tolerance_model/poster/poster_figs.py

Rebuilds what each figure needs from the notebook's own cells and saved results (no kernel needed):
§1 setup + the §14 reload cell (held-out §6g variants, tables), the 14h setup cell (the variant table with its
42 features), the cached 14h ablation predictions, and the official §16 test predictions.
"""
import os, sys, json, time
import numpy as np, pandas as pd

NB_DIR = "/Users/awlehr/LehrPDI/figures/notebooks"
POSTER = "/Users/awlehr/LehrPDI/figures/results/tolerance_model/poster"
SVG, PNG = os.path.join(POSTER, "svg"), os.path.join(POSTER, "png_preview")
os.makedirs(SVG, exist_ok=True); os.makedirs(PNG, exist_ok=True)
os.chdir(NB_DIR)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

C = json.load(open("05_tolerance_score_modelling.ipynb"))["cells"]
src = lambda prefix: next("".join(c["source"]) for c in C
                          if c["cell_type"] == "code" and "".join(c["source"]).startswith(prefix))
t0 = time.time()
exec(src("import sys, os, random"), globals())                      # §1 setup
exec(src("# OPTIONAL, after a kernel restart"), globals())         # §14 reload: §2-§6g state + saved tables
import torch
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score, roc_curve
probe_train = torch.load(os.path.join(DATA_DIR, "train.pt"), weights_only=True)
probe_val = torch.load(os.path.join(DATA_DIR, "val.pt"), weights_only=True)
AA20 = "ACDEFGHIKLMNPQRSTVWY"                                       # as in 14g
ANN_COLS = [c for c in residue_df.columns if c.startswith("ann_") and c != "ann_disorder_curated"] + ["n_annotations"]
exec(src("# 14h setup:"), globals())                                # lab: 15,339 variants x 42 features
plt.close("all")
print(f"state rebuilt in {time.time() - t0:.0f} s: {len(lab):,} labelled variants", flush=True)

# ------------------------------------------------------------------ poster style
plt.rcParams.update({"font.family": "Arial", "font.size": 15, "axes.titlesize": 17, "axes.labelsize": 15,
                     "xtick.labelsize": 13.5, "ytick.labelsize": 13.5, "legend.fontsize": 13.5, "svg.fonttype": "none",
                     "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#8a8984",
                     "axes.labelcolor": "#2b2b2a", "xtick.color": "#52514e", "ytick.color": "#52514e",
                     "figure.facecolor": "white", "axes.titlelocation": "left", "axes.titleweight": "bold"})
BLUE, RED, VIOLET, AQUA, YELLOW, ORANGE, GREEN = "#2a78d6", "#e34948", "#6250d6", "#1baf7a", "#eda100", "#eb6834", "#008300"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e5e4e0"
TINT = {BLUE: "#e3eefb", RED: "#fbe6e6", VIOLET: "#ebe8fa", AQUA: "#dff4ec", YELLOW: "#fcf1d6", ORANGE: "#fde9e0",
        GREEN: "#e2f1e2", "grey": "#f0efec"}
GENE, VAR_POS = "GLRA1", 300
official = os.path.join(OUT, "final_test")
pred = pd.read_parquet(os.path.join(official, "predictions.parquet")); pred["benign"] = pred.clinsig_group.eq("benign")
result = json.load(open(os.path.join(official, "results.json")))
PRIMARY = "logistic: all features"


def save(fig, name):
    fig.savefig(os.path.join(SVG, name + ".svg"), bbox_inches="tight")
    fig.savefig(os.path.join(PNG, name + ".png"), dpi=80, bbox_inches="tight")
    plt.close(fig); print("wrote", name, flush=True)


def gene_bootstrap(stat, genes, n=1000, seed=0):
    rows_of = pd.Series(np.arange(len(genes))).groupby(np.asarray(genes)).indices
    keys = np.array(list(rows_of)); rng = np.random.default_rng(seed); out = []
    for _ in range(n):
        v = stat(np.concatenate([rows_of[k] for k in rng.choice(keys, len(keys))]))
        if np.isfinite(v):
            out.append(v)
    return np.percentile(out, [2.5, 97.5])


# ================================================================== I1 / M2: scores by variant group (held-out, §6g)
v = val_scored.sort_values("database").drop_duplicates("id").copy()
v["group"] = np.where(v.database == "GnomAD", "gnomAD only", "ClinVar " + v.clinsig_group.astype(str))
GROUPS = [("gnomAD only", "gnomAD\n(healthy\npopulation)", BLUE), ("ClinVar benign", "ClinVar\nbenign", AQUA),
          ("ClinVar vus", "ClinVar\nVUS", YELLOW), ("ClinVar pathogenic", "ClinVar\npathogenic", RED)]
group_stats = {}


def by_group(col, ylabel, title, name, clip=None):
    fig, ax = plt.subplots(figsize=(9, 6.2))
    data = [v.loc[v.group == g, col].dropna() for g, _, _ in GROUPS]
    bp = ax.boxplot(data, showfliers=False, widths=0.55, patch_artist=True, medianprops=dict(color=INK, lw=2))
    for patch, (_, _, c) in zip(bp["boxes"], GROUPS):
        patch.set(facecolor=TINT[c], edgecolor=c, lw=2)
    for el in ["whiskers", "caps"]:
        for k, line in enumerate(bp[el]):
            line.set(color=GROUPS[k // 2][2], lw=1.6)
    lo_, hi_ = ax.get_ylim()
    for k, (d, (g, _, c)) in enumerate(zip(data, GROUPS), 1):
        ax.plot(k, d.mean(), "D", color=c, ms=10, mec="white", mew=1.5, zorder=4)
        ax.text(k + 0.31, d.mean(), f"mean\n{d.mean():+.2f}", va="center", fontsize=12.5, color=INK2)
        ax.text(k, lo_ - (hi_ - lo_) * 0.04, f"n = {len(d):,}", ha="center", va="top", fontsize=12.5, color=INK2)
    ax.set_ylim(lo_ - (hi_ - lo_) * 0.14, hi_)
    ax.set_xticks(range(1, len(GROUPS) + 1), [lab_ for _, lab_, _ in GROUPS])
    ax.set_ylabel(ylabel); ax.set_title(title)
    ax.set_xlim(0.45, len(GROUPS) + 0.85)
    gn, path, vus = data[0], data[3], data[2]
    stats = {"means": {g: float(d.mean()) for (g, _, _), d in zip(GROUPS, data)},
             "medians": {g: float(d.median()) for (g, _, _), d in zip(GROUPS, data)},
             "n": {g: int(len(d)) for (g, _, _), d in zip(GROUPS, data)},
             "auroc_gnomad_vs_pathogenic": roc_auc_score(np.r_[np.ones(len(gn)), np.zeros(len(path))], np.r_[gn, path]),
             "auroc_gnomad_vs_vus": roc_auc_score(np.r_[np.ones(len(gn)), np.zeros(len(vus))], np.r_[gn, vus]),
             "auroc_vus_vs_pathogenic": roc_auc_score(np.r_[np.ones(len(vus)), np.zeros(len(path))], np.r_[vus, path])}
    group_stats[col] = stats
    fig.tight_layout(); save(fig, name)


by_group("esm_llr", "ESM-1v LLR of the substitution\n(higher = looks harmless)",
         "ESM-1v separates healthy-population from pathogenic variants", "I1_esm_by_variant_group")
by_group("tolerance_v4_train", "tolerance at the variant's residue\n(higher = more tolerant)",
         "Tolerance, built without these variants, separates them too", "M2_tolerance_by_variant_group")

# ================================================================== M1: building the tolerance score (GLRA1)
from scipy.stats import binom  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
tol = pd.read_parquet(os.path.join(OUT, "tolerance_v4.parquet")).set_index("gene").loc[GENE]
L = len(tol.sequence); x = np.arange(1, L + 1)
gv = tolerance_variants[tolerance_variants.gene == GENE]
gn_cnt = np.bincount(gv.position[gv.database == "GnomAD"] - 1, minlength=L)
cv_cnt = np.bincount(gv.position[gv.database == "ClinVar"] - 1, minlength=L)
norm = lambda a: (np.asarray(a) - np.min(a)) / (np.ptp(a) or 1)
dis = np.asarray(tol.disorder_scores); t = np.asarray(tol.tolerance_v4)
fig, axes = plt.subplots(3, 1, figsize=(13, 9.5), sharex=True, gridspec_kw={"height_ratios": [1, 1, 1.1]})
ax = axes[0]
ax.vlines(x[gn_cnt > 0], 0, gn_cnt[gn_cnt > 0], color=BLUE, lw=1.4)
ax.vlines(x[cv_cnt > 0], 0, -cv_cnt[cv_cnt > 0], color=RED, lw=1.4)
ax.axhline(0, color=MUTED, lw=0.9); ax.set_ylim(-cv_cnt.max() - 1, gn_cnt.max() + 4.5)
ax.set_ylabel("variants\nper residue")
ax.legend(handles=[Line2D([], [], color=BLUE, lw=3, label=f"gnomAD (up): {int((gv.database == 'GnomAD').sum())} variants, healthy population"),
                   Line2D([], [], color=RED, lw=3, label=f"ClinVar (down): {int((gv.database == 'ClinVar').sum())} variants, any classification")],
          frameon=False, loc="upper left", fontsize=13)
ax.set_title(f"1. Observed missense variants ({GENE}, {L} residues)")
ax = axes[1]
ax.plot(x, norm(tol.gnomad_y_v4), color=BLUE, lw=2.2, label="gnomAD density")
ax.plot(x, norm(tol.clinvar_y_v4), color=RED, lw=2.2, label="ClinVar density")
ax.set_ylabel("smoothed density\n(scaled 0–1)"); ax.legend(frameon=False, ncol=2, loc="upper left", fontsize=13)
ax.set_ylim(0, 1.25)
ax.set_title("2. Weighted by substitution severity, smoothed along the protein")
ax = axes[2]
ax.fill_between(x, 0, 1, where=dis >= 0.5, transform=ax.get_xaxis_transform(), color=TINT[YELLOW], lw=0, step="mid")
ax.plot(x, t, color=INK, lw=1.8)
ax.fill_between(x, t, 0, where=t > 0, color=BLUE, alpha=0.3, interpolate=True)
ax.fill_between(x, t, 0, where=t < 0, color=RED, alpha=0.3, interpolate=True)
ax.axhline(0, color=MUTED, lw=0.9)
ax.annotate(f"p.A{VAR_POS}P\n(pathogenic)", (VAR_POS, t[VAR_POS - 1]), (VAR_POS + 22, -0.55), fontsize=13,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=1), color=INK2)
ax.text(4, 0.43, "tolerated", color=BLUE, fontsize=13, weight="bold", va="top")
ax.text(60, -0.5, "intolerant", color=RED, fontsize=13, weight="bold")
ax.set_ylim(-0.72, 0.5); ax.set_ylabel("tolerance\n(gnomAD − ClinVar)"); ax.set_xlabel(f"{GENE} residue")
ax.set_title("3. Tolerance = gnomAD density − ClinVar density   (yellow: disordered regions)")
fig.tight_layout(); save(fig, "M1_tolerance_score_build")

# ================================================================== M3: structural and annotation tracks (GLRA1)
rg = residue_df[residue_df.gene == GENE].sort_values("pos").reset_index(drop=True)
pos = rg.pos.to_numpy()
tracks = [("ann_signal", "signal peptide", MUTED), ("ann_domain", "domain (InterPro)", VIOLET),
          ("ann_transmembrane", "transmembrane", ORANGE), ("ann_ligand_site", "ligand site", GREEN),
          ("ann_disulfide", "disulfide", AQUA), ("ann_ptm", "PTM site", YELLOW),
          ("ann_motif", "motif", RED), ("ann_motif_predicted", "predicted motif", RED)]
tracks = [tr for tr in tracks if tr[0] in rg and rg[tr[0]].fillna(False).astype(bool).any()]
fig, axes = plt.subplots(5, 1, figsize=(13, 10.5), sharex=True,
                         gridspec_kw={"height_ratios": [1.1, 0.9, 1.1, 0.42, 0.32 * len(tracks) + 0.3]})
ax = axes[0]
tc = rg.tolerance_v4_centered.to_numpy()
ax.plot(pos, tc, color=INK, lw=1.6)
ax.fill_between(pos, tc, 0, where=tc > 0, color=BLUE, alpha=0.3, interpolate=True)
ax.fill_between(pos, tc, 0, where=tc < 0, color=RED, alpha=0.3, interpolate=True)
ax.axhline(0, color=MUTED, lw=0.9); ax.set_ylabel("tolerance\n(centred)")
ax.set_title(f"{GENE}: tolerance and the structural features at every residue")
ax = axes[1]
d_ = rg.disorder.to_numpy()
ax.fill_between(pos, d_, color=YELLOW, alpha=0.55, lw=0)
ax.plot(pos, rg.binding.to_numpy(), color=ORANGE, lw=1.4)
ax.axhline(0.5, color=MUTED, lw=0.9, ls="--")
ax.text(pos[-1] + 6, 0.62, "disordered", ha="left", va="center", fontsize=12, color=INK2, clip_on=False)
ax.text(pos[-1] + 6, 0.38, "structured", ha="left", va="center", fontsize=12, color=INK2, clip_on=False)
ax.set_ylim(0, 1.05); ax.set_ylabel("AIUPred\nscore")
ax.legend(handles=[Line2D([], [], color=YELLOW, lw=6, alpha=0.6, label="disorder"),
                   Line2D([], [], color=ORANGE, lw=2, label="binding")], frameon=False, ncol=2, loc="upper center",
          bbox_to_anchor=(0.42, 1.02), fontsize=12.5)
ax = axes[2]
rsa = rg.rsa.to_numpy(dtype=float)
cls = np.where(rsa < 0.2, AQUA, np.where(rsa <= 0.4, "#8fcfb6", "#cfe9de"))
ok_ = np.isfinite(rsa)
ax.bar(pos[ok_], rsa[ok_], width=1.0, color=cls[ok_], lw=0)
ax.axhline(0.2, color=MUTED, lw=0.9, ls="--"); ax.axhline(0.4, color=MUTED, lw=0.9, ls=":")
ax.set_ylim(0, 1.05); ax.set_ylabel("RSA")
ax.text(pos[-1] + 6, 0.7, "surface\n(> 0.4)", ha="left", va="center", fontsize=12, color=INK2, clip_on=False)
ax.text(pos[-1] + 6, 0.1, "core\n(< 0.2)", ha="left", va="center", fontsize=12, color=INK2, clip_on=False)
ax = axes[3]
iface = rg.iface_any.fillna(False).astype(bool).to_numpy()
ax.bar(pos[iface], 1, width=1.0, color=VIOLET, lw=0)
ax.set_yticks([]); ax.set_ylim(0, 1); ax.set_ylabel("interface", rotation=0, ha="right", va="center")
ax.spines["left"].set_visible(False)
ax = axes[4]
for k, (col, labl, colr) in enumerate(tracks):
    m = rg[col].fillna(False).astype(bool).to_numpy()
    ax.bar(pos[m], 0.7, bottom=len(tracks) - 1 - k + 0.15, width=1.0, color=colr, lw=0)
ax.set_yticks([len(tracks) - 1 - k + 0.5 for k in range(len(tracks))], [labl for _, labl, _ in tracks], fontsize=12.5)
ax.set_ylim(0, len(tracks)); ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
ax.set_xlabel(f"{GENE} residue")
for a in axes:
    a.axvline(VAR_POS, color=INK2, lw=1, ls=":")
axes[0].text(VAR_POS + 3, axes[0].get_ylim()[1] * 0.9, f"p.A{VAR_POS}P", fontsize=12, color=INK2, va="top")
fig.tight_layout(); save(fig, "M3_structural_feature_tracks")

# ================================================================== M4: one row per residue -> one row per variant
rw = residue_df[(residue_df.gene == GENE) & residue_df.pos.between(VAR_POS - 3, VAR_POS + 3)].set_index("pos")
cols6 = [("tolerance_v4_centered", "tolerance\n(centred)", BLUE, "{:+.2f}"), ("esm_mean_llr", "ESM mean\nLLR", VIOLET, "{:.1f}"),
         ("esm_entropy", "ESM\nentropy", VIOLET, "{:.2f}"), ("rsa", "RSA", AQUA, "{:.2f}"), ("plddt", "pLDDT", AQUA, "{:.0f}"),
         ("iface_any", "inter-\nface", AQUA, "{}"), ("ann_transmembrane", "trans-\nmembrane", YELLOW, "{}"),
         ("disorder", "disorder", YELLOW, "{:.2f}"), ("msa_conservation", "conser-\nvation", YELLOW, "{:.1f}")]
fig, ax = plt.subplots(figsize=(15, 7.4)); ax.set_xlim(-0.01, 1.01); ax.set_ylim(-0.01, 1.01); ax.axis("off")
ax.text(0, 0.995, f"One row per residue: every source lined up on the same numbering ({GENE} residues {VAR_POS - 3}–{VAR_POS + 3})",
        fontsize=16, weight="bold", va="top")
heads = [("pos", "residue", "grey"), ("aa", "amino\nacid", "grey")] + [(c, h, colr) for c, h, colr, _ in cols6]
x0, cw, rh, top = 0.0, 0.0885, 0.06, 0.9
for j, (_, h, colr) in enumerate(heads):
    ax.add_patch(plt.Rectangle((x0 + j * cw, top - rh * 1.45), cw - 0.004, rh * 1.45, color=TINT.get(colr, TINT["grey"])))
    ax.text(x0 + j * cw + cw / 2, top - rh * 0.72, h, ha="center", va="center", fontsize=13, weight="bold")
for i, p in enumerate(rw.index):
    y = top - rh * 1.45 - (i + 1) * rh
    hl = p == VAR_POS
    if hl:
        ax.add_patch(plt.Rectangle((x0, y), cw * len(heads) - 0.004, rh, color="#fff1d6", lw=0))
    vals = [str(p), rw.loc[p, "aa"]] + [("–" if pd.isna(rw.loc[p, c]) else
                                       (("yes" if rw.loc[p, c] else "no") if fmt == "{}" else fmt.format(rw.loc[p, c])))
                                      for c, _, _, fmt in cols6]
    for j, val in enumerate(vals):
        ax.text(x0 + j * cw + cw / 2, y + rh / 2, val, ha="center", va="center", fontsize=13.5,
                weight="bold" if hl else "normal")
src_y = top - rh * 1.45 - 8.4 * rh
for k, (colr, labl) in enumerate([("grey", "sequence"), (BLUE, "human variants"), (VIOLET, "ESM-1v"), (AQUA, "structure"),
                                  (YELLOW, "annotations, disorder, conservation")]):
    ax.add_patch(plt.Rectangle((k * 0.175, src_y), 0.016, 0.034, color=TINT.get(colr, TINT["grey"]), ec="#9b9a96", lw=0.6))
    ax.text(k * 0.175 + 0.022, src_y + 0.017, labl, va="center", fontsize=13, color=INK2)
vrow = lab[(lab.gene == GENE) & (lab.protein_consequence == f"p.A{VAR_POS}P")].iloc[0]
ax.add_patch(FancyBboxPatch((0, 0.01), 0.975, 0.2, boxstyle="round,pad=0.006,rounding_size=0.012", fc=TINT["grey"], ec="#9b9a96", lw=0.9))
ax.text(0.012, 0.19, f"One row per variant: p.A{VAR_POS}P (ClinVar: {vrow.clinsig_group})", fontsize=15, weight="bold", va="top")
ax.text(0.012, 0.135, f"=  the residue-{VAR_POS} row above, looked up by (gene, position)\n"
        f"+  the swap itself: ESM LLR {vrow.esm_llr:.1f}, BLOSUM62 / Grantham severity, hydrophobicity change\n"
        f"+  human-variant features with residue {VAR_POS}'s own variants hidden: tolerance, nearby gnomAD / ClinVar counts",
        fontsize=13.2, va="top", color=INK2, linespacing=1.45)
ax.annotate("", (0.5, 0.215), (0.5, src_y - 0.01), arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=2.2))
save(fig, "M4_residue_to_variant_row")

# ================================================================== R1: model table
rows_t = [  # model, inputs, CV AUROC, test AUROC, test AUROC where ESM unsure
    ("ESM-1v LLR", "no training: evolution only", "0.906", "0.620", "0.545"),
    ("Tolerance score", "no training: human variants only", "0.837", "0.828", "0.839"),
    ("Logistic regression", "ESM LLR + tolerance", "0.933", "0.826", "0.833"),
    ("Logistic regression", "all 42 features", "0.961", "0.890", "0.879"),
    ("Gradient-boosted trees", "all 42 features", "0.958", "0.860", "0.840"),
    ("MLP neural network", "42 features + 64 ESM components", "0.959", "not tested", "not tested"),
    ("Window CNN", "±50 residues: sequence, raw variant counts", "0.943", "not tested", "not tested"),
    ("Window BiLSTM", "±50 residues: sequence, raw variant counts", "0.950", "not tested", "not tested"),
]
fig, ax = plt.subplots(figsize=(14, 6.8)); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
colx = [0.0, 0.215, 0.535, 0.69, 0.845]; widths = [0.215, 0.32, 0.155, 0.155, 0.155]
head1 = ["Model", "Inputs", "Cross-\nvalidation\n(15,339)", "Untouched\ngenes\n(1,357)", "Untouched genes,\nESM unsure\n(668)"]
rh = 0.09; ytop = 0.97
for j, h in enumerate(head1):
    ax.add_patch(plt.Rectangle((colx[j], ytop - rh * 1.9), widths[j] - 0.004, rh * 1.9, color=TINT["grey"]))
    ax.text(colx[j] + (0.01 if j < 2 else widths[j] / 2), ytop - rh * 0.95, h, ha="left" if j < 2 else "center",
            va="center", fontsize=13, weight="bold", linespacing=1.15)
for i, r in enumerate(rows_t):
    y = ytop - rh * 1.9 - (i + 1) * rh
    best = r[1] == "all 42 features" and r[0].startswith("Logistic")
    if best:
        ax.add_patch(plt.Rectangle((0, y), 0.997, rh, color=TINT[BLUE], lw=0))
    for j, cell in enumerate(r):
        ax.text(colx[j] + (0.01 if j < 2 else widths[j] / 2), y + rh / 2, cell, ha="left" if j < 2 else "center",
                va="center", fontsize=13.5, color=MUTED if cell == "not tested" else INK, weight="bold" if best else "normal")
    ax.plot([0, 0.997], [y, y], color=GRID, lw=0.8)
ax.text(0, ytop - rh * 1.9 - 8.6 * rh, "AUROC: 0.5 = chance, 1.0 = perfect separation of benign and pathogenic ClinVar variants. "
        "Highlighted: the model tested on untouched genes.", fontsize=12.5, color=INK2, va="top")
save(fig, "R1_model_auroc_table")

# ================================================================== R2: logistic regression on its S-curve (test genes)
p_t = pred[PRIMARY].clip(1e-6, 1 - 1e-6).to_numpy(); z = np.log(p_t / (1 - p_t)); y_t = pred.benign.to_numpy()
fig, axes = plt.subplots(1, 2, figsize=(15, 5.6), gridspec_kw={"width_ratios": [1.15, 1]})
ax = axes[0]; ax.set_xlim(-0.01, 1.01); ax.set_ylim(-0.01, 1.01); ax.axis("off")


def pbox(ax, x, y, w, h, title, body, colr, ts=14, bs=12.5):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.015", fc=TINT.get(colr, colr), ec="#9b9a96", lw=0.9))
    ax.text(x + w / 2, y + h * 0.66, title, ha="center", va="center", fontsize=ts, weight="bold")
    if body:
        ax.text(x + w / 2, y + h * 0.3, body, ha="center", va="center", fontsize=bs, color=INK2)


feats = ["ESM LLR", "gnomAD nearby", "ClinVar nearby", "structure", "… 38 more"]
for k, f in enumerate(feats):
    yk = 0.84 - k * 0.165
    pbox(ax, 0.0, yk - 0.055, 0.25, 0.11, f, "", "grey", ts=13)
    ax.annotate("", (0.36, 0.5), (0.25, yk), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.2))
pbox(ax, 0.36, 0.37, 0.2, 0.26, "weighted\nsum  z", "", GREEN)
ax.annotate("", (0.62, 0.5), (0.56, 0.5), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.4))
pbox(ax, 0.62, 0.37, 0.17, 0.26, "S-curve", "1 / (1 + e^(−z))", GREEN)
ax.annotate("", (0.84, 0.5), (0.79, 0.5), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.4))
pbox(ax, 0.84, 0.37, 0.16, 0.26, "P(benign)", "per variant", GREEN)
ax.set_title("Logistic regression: 42 features → one probability")
ax = axes[1]
zz = np.linspace(-8, 8, 300)
ax.plot(zz, 1 / (1 + np.exp(-zz)), color=INK2, lw=2)
rng_j = np.random.default_rng(0)
for is_b, c, labl in [(False, RED, "pathogenic"), (True, BLUE, "benign")]:
    s_ = y_t == is_b; zc = np.clip(z[s_], -8, 8)
    ax.scatter(zc, 1 / (1 + np.exp(-zc)) + rng_j.normal(0, 0.014, s_.sum()), s=12, alpha=0.4, color=c, lw=0,
               label=f"test {labl} ({s_.sum()})")
ax.axhline(0.5, color=MUTED, ls="--", lw=1)
ax.set_xlabel("z (weighted sum of the features)"); ax.set_ylabel("P(benign)")
ax.set_title(f"{len(pred):,} variants in 85 untouched genes")
ax.legend(frameon=False, markerscale=2.5, loc="upper left")
fig.tight_layout(); save(fig, "R2_logistic_regression_curve")

# ================================================================== R2b: ROC curves on the test genes
UNSURE = lab.loc[lab.benign, "esm_llr"].median()                              # from training variants, as in §16
unsure_t = (pred["ESM-1v LLR"] > UNSURE).to_numpy()
fig, axes = plt.subplots(1, 2, figsize=(14, 6.2), sharey=True)
for ax, subset, title in [(axes[0], np.ones(len(pred), bool), f"All test variants ({len(pred):,})"),
                          (axes[1], unsure_t, f"Where ESM-1v is unsure ({unsure_t.sum()})")]:
    for name, c, lw, labl in [("ESM-1v LLR", MUTED, 2.2, "ESM-1v alone"), ("tolerance (own residue hidden)", ORANGE, 2.2, "tolerance alone"),
                              (PRIMARY, BLUE, 3.0, "logistic regression, 42 features")]:
        s_ = pred[name].to_numpy()[subset]
        fpr, tpr, _ = roc_curve(y_t[subset], s_)
        ax.plot(fpr, tpr, color=c, lw=lw, label=f"{labl}: {roc_auc_score(y_t[subset], s_):.2f}")
    ax.plot([0, 1], [0, 1], color=GRID, ls="--", lw=1.2)
    ax.set_xlabel("pathogenic variants called benign (share)")
    ax.set_title(title); ax.legend(frameon=False, loc="lower right", title="AUROC", title_fontsize=13.5)
axes[0].set_ylabel("benign variants called benign (share)")
fig.tight_layout(); save(fig, "R2b_roc_untouched_genes")

# ================================================================== R3: untouched genes, model vs ESM-1v
sc = result["scores"]; qs = result["questions"]
rowsR = [("All variants", "AUROC", "Q1: primary beats ESM-1v on all variants"),
         ("Within genes\n(8 genes)", "within-gene AUROC", "Q2: primary beats ESM-1v within genes"),
         ("Where ESM is unsure\n(668 variants)", "AUROC where ESM is unsure", "Q3: primary beats ESM-1v where ESM is unsure")]
fig, ax = plt.subplots(figsize=(12, 5.2))
for k, (labl, key, q) in enumerate(rowsR):
    yk = len(rowsR) - 1 - k
    a, b = sc["ESM-1v LLR"][key], sc[PRIMARY][key]
    ax.plot([a, b], [yk, yk], color=GRID, lw=5, zorder=1, solid_capstyle="round")
    ax.scatter(a, yk, s=170, color=MUTED, zorder=3); ax.scatter(b, yk, s=210, color=BLUE, zorder=3)
    ax.text(a - 0.012, yk + 0.2, f"{a:.2f}", ha="right", fontsize=13.5, color=INK2)
    ax.text(b + 0.012, yk + 0.2, f"{b:.2f}", ha="left", fontsize=13.5, color=BLUE, weight="bold")
    qq = qs[q]
    verdict = "significant" if qq["passes"] else "not significant"
    ax.text(1.005, yk - 0.05, f"Δ {qq['difference']:+.3f}\n95% CI {qq['lo']:+.3f} to {qq['hi']:+.3f}\n{verdict}",
            va="center", fontsize=12.5, color=INK if qq["passes"] else MUTED, linespacing=1.3)
ax.set_yticks(range(len(rowsR))[::-1], [r[0] for r in rowsR])
ax.set_xlim(0.5, 1.0); ax.set_ylim(-0.6, len(rowsR) - 0.3)
ax.set_xlabel("AUROC on 85 untouched genes (1,357 ClinVar variants)")
ax.scatter([], [], s=110, color=MUTED, label="ESM-1v alone"); ax.scatter([], [], s=130, color=BLUE, label="logistic regression, 42 features")
ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0, -0.42), ncol=2)
ax.set_title("Tested once on genes never used for training or model choice")
fig.tight_layout(); save(fig, "R3_untouched_genes_vs_esm")

# ================================================================== R4: feature groups (Figure 23, poster version)
group_pred = pd.read_parquet(os.path.join(OUT, "baselines", "ablation_variant_features.parquet"))
y_lab_ = lab.benign.to_numpy(); genes_lab = lab.gene.to_numpy()
esm_only = group_pred["only ESM scores"].to_numpy()


def within_auc(s):
    rows_of = pd.Series(np.arange(len(s))).groupby(genes_lab).indices
    vals = [roc_auc_score(y_lab_[i], s[i]) for i in rows_of.values() if 3 <= y_lab_[i].sum() <= len(i) - 3]
    return float(np.mean(vals))


NICE = {"ESM scores": "ESM scores", "human variants nearby": "human variants nearby\n(tolerance, gnomAD / ClinVar counts)",
        "gene-wide variant rates": "gene-wide variant rates", "substitution type": "substitution chemistry",
        "structure (§7)": "structure (RSA, pLDDT,\ninterfaces)", "annotations (§8)": "annotations (domains,\nmotifs, PTMs, …)",
        "disorder + conservation": "disorder + conservation"}
rows4 = []
for g in VARIANT_GROUPS:
    alone = group_pred[f"only {g}"].to_numpy()
    r = {"group": g, "alone": roc_auc_score(y_lab_, alone), "alone_within": within_auc(alone)}
    if g != "ESM scores":
        plus = group_pred[f"ESM + {g}"].to_numpy()
        r["gain"] = roc_auc_score(y_lab_, plus) - roc_auc_score(y_lab_, esm_only)
        stat = lambda rows, a=plus: roc_auc_score(y_lab_[rows], a[rows]) - roc_auc_score(y_lab_[rows], esm_only[rows])
        r["gain_lo"], r["gain_hi"] = gene_bootstrap(stat, genes_lab, n=500)
        r["gain_within"] = within_auc(plus) - within_auc(esm_only)
    rows4.append(r)
abl = pd.DataFrame(rows4).set_index("group")
order = abl.sort_values("alone").index
fig, axes = plt.subplots(1, 2, figsize=(16, 6.4), gridspec_kw={"width_ratios": [1, 1]})
yy = np.arange(len(order))
ax = axes[0]
ax.barh(yy + 0.2, abl.loc[order, "alone"], height=0.38, color=BLUE, label="all variants")
ax.barh(yy - 0.2, abl.loc[order, "alone_within"], height=0.38, color="#86b6ef", label="within genes")
for k, g in enumerate(order):
    ax.text(abl.loc[g, "alone"] + 0.005, k + 0.2, f"{abl.loc[g, 'alone']:.2f}", va="center", fontsize=12)
ax.set_xlim(0.45, 1.0); ax.axvline(0.5, color=GRID, lw=1)
ax.set_yticks(yy, [NICE[g] for g in order], fontsize=13)
ax.set_xlabel("AUROC with this group alone"); ax.legend(frameon=False, loc="lower right")
ax.set_title("A. Each feature group on its own")
ax = axes[1]
order2 = [g for g in order if g != "ESM scores"]
yy2 = np.arange(len(order2))
gain = abl.loc[order2, "gain"].to_numpy(float)
ax.barh(yy2 + 0.2, gain, height=0.38, color=BLUE, label="all variants (95% CI)",
        xerr=[gain - abl.loc[order2, "gain_lo"].to_numpy(float), abl.loc[order2, "gain_hi"].to_numpy(float) - gain],
        error_kw=dict(ecolor=INK2, lw=1.3, capsize=3))
ax.barh(yy2 - 0.2, abl.loc[order2, "gain_within"].to_numpy(float), height=0.38, color="#86b6ef", label="within genes")
ax.axvline(0, color=INK2, lw=1)
ax.set_yticks(yy2, [NICE[g] for g in order2], fontsize=13)
ax.set_xlabel("AUROC gained when added to the ESM scores"); ax.legend(frameon=False, loc="lower right")
ax.set_title("B. Added to ESM: what each group contributes")
fig.tight_layout(); save(fig, "R4_feature_groups")

# ================================================================== R5: the final model's weights (Figure 26, poster version)
final = make_pipeline(SimpleImputer(strategy="median", add_indicator=True), StandardScaler(), LogisticRegression(max_iter=5000))
final.fit(lab[ALL_VAR_FEATURES].to_numpy(), lab.benign.to_numpy())
names = final[:-1].get_feature_names_out()
coef = pd.Series(final[-1].coef_[0], index=names)
weights = pd.Series({ALL_VAR_FEATURES[int(n[1:])]: w for n, w in coef.items() if n.startswith("x")})
PLAIN = {"esm_llr": "ESM LLR of the substitution", "gnomad_within_15": "gnomAD variants within ±15 residues",
         "clinvar_within_15": "ClinVar variants within ±15 residues", "msa_conservation": "conservation (MobiDB)",
         "clinvar_per_residue": "ClinVar variants per residue (whole gene)", "iface_any": "binding interface",
         "rsa": "solvent accessibility (RSA)", "blosum62": "BLOSUM62 dissimilarity", "disorder": "disorder (AIUPred)",
         "clinvar_within_5": "ClinVar variants within ±5 residues", "tolerance": "tolerance score",
         "esm_entropy": "ESM entropy", "gnomad_per_residue": "gnomAD variants per residue (whole gene)",
         "binding": "binding region (AIUPred)", "gnomad_within_5": "gnomAD variants within ±5 residues",
         "esm_mean_llr": "ESM mean LLR", "esm_wt_lp": "ESM wild-type log-prob", "grantham": "Grantham distance",
         "hydrophobicity_change": "hydrophobicity change", "plddt": "AlphaFold pLDDT", "gnomad_density": "gnomAD density",
         "clinvar_density": "ClinVar density"}
top = weights.reindex(weights.abs().sort_values(ascending=False).index[:15])
fig, ax = plt.subplots(figsize=(13, 7.4))
yy = np.arange(len(top))[::-1]
ax.barh(yy, top, color=[BLUE if w > 0 else RED for w in top], height=0.7)
for y_i, (f, w) in zip(yy, top.items()):
    ax.text(w + (0.03 if w > 0 else -0.03), y_i, f"{w:+.2f}", va="center", ha="left" if w > 0 else "right", fontsize=12, color=INK2)
ax.set_yticks(yy, [PLAIN.get(f, f) for f in top.index], fontsize=13)
ax.axvline(0, color=INK2, lw=1)
ax.set_xlim(top.min() - 0.35, top.max() + 0.35)
ax.set_xlabel("weight per standard deviation of the feature")
ax.text(top.max() + 0.3, -1.25, "→ pushes towards benign", ha="right", color=BLUE, fontsize=13, weight="bold")
ax.text(top.min() - 0.3, -1.25, "pushes towards pathogenic ←", ha="left", color=RED, fontsize=13, weight="bold")
ax.set_ylim(-1.8, len(top) - 0.4)
ax.set_title("What the final model relies on: its 15 largest weights")
fig.tight_layout(); save(fig, "R5_model_weights")

json.dump({"group_stats": group_stats, "feature_groups": abl.round(4).to_dict(orient="index"),
           "weights_top15": top.round(3).to_dict(), "unsure_threshold": float(UNSURE)},
          open(os.path.join(POSTER, "poster_numbers.json"), "w"), indent=1, default=float)
print(json.dumps({k: {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in d.items() if not isinstance(vv, dict)}
                  for k, d in group_stats.items()}, indent=1))
print({k: round(v, 3) for k, v in group_stats["tolerance_v4_train"]["means"].items()}, {k: round(v, 3) for k, v in group_stats["esm_llr"]["means"].items()})
print(abl.round(3).to_string()); print(top.round(2).to_string())
print(f"done in {time.time() - t0:.0f} s")
