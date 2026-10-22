"""New explanatory figures for the training guide (report.html -> PDF). Run with the LehrPDI environment:

    python make_figs.py
Reads the saved results in results/tolerance_model/ (no notebook kernel needed) and writes figs/*.png.
"""
import os, sys, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from scipy.stats import binom, spearmanr
from sklearn.metrics import roc_auc_score, roc_curve

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)                                   # results/tolerance_model
FIG = os.path.join(HERE, "figs")
MODEL_DIR = os.path.abspath(os.path.join(OUT, "..", "..", "..", "Model Training"))
plt.rcParams.update({"font.family": "Arial", "font.size": 9.5, "axes.titlesize": 10.5, "axes.labelsize": 9.5,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#8a8984",
                     "axes.labelcolor": "#2b2b2a", "xtick.color": "#52514e", "ytick.color": "#52514e",
                     "savefig.dpi": 200, "savefig.bbox": "tight", "figure.facecolor": "white"})
# reference palette (dataviz skill), by role
BLUE, RED, VIOLET, AQUA, YELLOW, ORANGE, GREEN = "#2a78d6", "#e34948", "#6250d6", "#1baf7a", "#eda100", "#eb6834", "#008300"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e5e4e0"
TINT = {BLUE: "#e3eefb", RED: "#fbe6e6", VIOLET: "#ebe8fa", AQUA: "#dff4ec", YELLOW: "#fcf1d6", ORANGE: "#fde9e0",
        GREEN: "#e2f1e2", "grey": "#f0efec"}
GENE, VAR_POS = "GLRA1", 300


def save(fig, name):
    fig.savefig(os.path.join(FIG, name)); plt.close(fig); print("wrote", name)


def box(ax, x, y, w, h, title, body="", colour="grey", title_size=10, body_size=8.6, title_colour=INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.012",
                                fc=TINT.get(colour, colour), ec="#9b9a96", lw=0.8))
    ax.text(x + 0.012, y + h - 0.02, title, fontsize=title_size, weight="bold", va="top", color=title_colour)
    if body:
        ax.text(x + 0.012, y + h - 0.02 - 0.045 * title_size / 10 * (title.count("\n") + 1), body, fontsize=body_size,
                va="top", color=INK2, linespacing=1.35)


def arrow(ax, x0, y0, x1, y1, colour=MUTED, lw=1.2, style="-|>", text=None, text_offset=(0, 0.015), fs=8.5):
    ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle=style, color=colour, lw=lw, shrinkA=1, shrinkB=1))
    if text:
        ax.text((x0 + x1) / 2 + text_offset[0], (y0 + y1) / 2 + text_offset[1], text, fontsize=fs, color=colour,
                ha="center", va="bottom")


def canvas(w, h):
    fig, ax = plt.subplots(figsize=(w, h)); ax.set_xlim(-0.012, 1.012); ax.set_ylim(-0.012, 1.012); ax.axis("off")
    return fig, ax


# ------------------------------------------------------------------------------------------- data
tol = pd.read_parquet(os.path.join(OUT, "tolerance_v4.parquet"))
row = tol.set_index("gene").loc[GENE]
var = pd.read_parquet(os.path.join(OUT, "variant_scores.parquet"),
                      columns=["gene", "id", "position", "reference", "alternate", "protein_consequence", "database",
                               "total_count", "clinsig_group", "stars", "esm_llr"])
res = pd.read_parquet(os.path.join(OUT, "residue_scores.parquet"))
pred = pd.read_parquet(os.path.join(OUT, "final_test", "predictions.parquet"))
results = json.load(open(os.path.join(OUT, "final_test", "results.json")))
split_of = json.load(open(os.path.join(OUT, "split_manifest.json")))["split"]

# =========================================================================== F1: the two designs
fig, ax = canvas(11.5, 7.2)


def design_row(y0, letter, title, inp, model, output, target, target_colour, status):
    ax.text(0.0, y0 + 0.455, f"{letter}.  {title}", fontsize=12, weight="bold", va="top", color=INK)
    box(ax, 0.00, y0 + 0.16, 0.17, 0.24, "Input", inp, "grey", body_size=8.6)
    box(ax, 0.21, y0 + 0.16, 0.15, 0.24, "Model", model, VIOLET, body_size=8.6)
    box(ax, 0.40, y0 + 0.16, 0.17, 0.24, "Output", output, GREEN, body_size=8.6)
    arrow(ax, 0.17, y0 + 0.28, 0.21, y0 + 0.28); arrow(ax, 0.36, y0 + 0.28, 0.40, y0 + 0.28)
    box(ax, 0.21, y0 + 0.0, 0.36, 0.11, "The answer it learns from (the target)", target, target_colour, body_size=8.6)
    arrow(ax, 0.485, y0 + 0.11, 0.485, y0 + 0.16, target_colour, lw=1.6)
    box(ax, 0.61, y0 + 0.0, 0.39, 0.40, "Where it stands", status, "grey", body_size=8.6)


design_row(0.53, "A", "Your original design: the model predicts TOLERANCE",
           "the protein sequence\nonly, read in context\n(no variant data,\nno structure)",
           "ESM language\nmodel + a small\n'head' network",
           "a predicted tolerance\nfor EVERY residue,\neven where no\nvariant was seen",
           "the observed tolerance track: gnomAD density minus ClinVar\ndensity (§3), centred on each gene's median", BLUE,
           "Built in §11–§13 (data, model, training loop).\nTested on this laptop in 14c and 14f with ESM kept frozen.\n\n"
           "Best: within-gene Spearman 0.125 on unseen genes,\nagainst a noise ceiling of 0.39 (the most any model\ncould reach "
           "with this noisy target).\n\nWith the 'pathogenic-only' target recipe (14g) the ceiling\nrises to 0.61 and the "
           "simplest model to 0.32.\n\nThe fine-tuning run on Biowulf (§15) has not been done.")
design_row(0.02, "B", "Tested in 14h and §16: the model predicts the CLINVAR LABEL",
           "42 numbers per variant:\nESM scores, OBSERVED\ntolerance nearby (own\nresidue hidden), structure,\n"
           "annotations, disorder,\nchemistry of the swap",
           "logistic\nregression\n(a weighted\nsum)",
           "P(benign) for\nEACH variant",
           "the ClinVar classification: pathogenic / likely pathogenic\nvs benign / likely benign (at least 1 review star)", RED,
           "Cross-validated in 14h (15,339 variants, 716 genes).\nTested once on 85 untouched genes in §16 (1,357 variants).\n\n"
           "Test AUROC 0.89, against 0.62 for ESM-1v alone;\n0.88 vs 0.55 where ESM is unsure.\n\n"
           "Here tolerance is one of the INPUTS: evidence the model\nweighs, not the answer it produces. Its output can still be\n"
           "turned into a per-residue tolerance score (section 9).")
ax.plot([0, 1], [0.505, 0.505], color=GRID, lw=1.2)
save(fig, "f01_two_designs.png")

# =========================================================================== F2: data-use map
rows = [  # data item, [builds tolerance, A: networks, A: trees 14g, B: variant model, labels/evaluation]
    ("Protein sequence (MANE transcript)", ["track length", "input (tokens)", "amino acid, position", "chemistry of the swap", ""]),
    ("gnomAD variants: where they are", ["gnomAD density", "", "", "counts nearby, tolerance", ""]),
    ("gnomAD allele count", ["filter: seen ≥ 2×", "", "", "", ""]),
    ("gnomAD allele frequency", ["", "", "", "", ""]),
    ("gnomAD zygosity (X / Y genes)", ["filter: X-linked rule", "", "", "", ""]),
    ("ClinVar variants: where they are", ["ClinVar density", "", "", "counts nearby, tolerance", ""]),
    ("ClinVar classification (P / B / VUS)", ["tested only (14g)", "", "", "", "THE LABEL (P vs B)"]),
    ("ClinVar review stars", ["", "", "", "", "filter: ≥ 1 star"]),
    ("AIUPred disorder and binding", ["chooses the smoothing", "", "3 features", "3 features", ""]),
    ("ESM-1v scores (log-probabilities)", ["", "", "3 scores", "4 scores incl. LLR", "baseline to beat"]),
    ("ESM-1v representations (1,280 / residue)", ["", "input (frozen)", "64 summaries", "", ""]),
    ("ESM-2 weights, fine-tuned with LoRA", ["", "planned (Biowulf)", "", "", ""]),
    ("Structure: RSA, pLDDT, interfaces (§7)", ["", "", "7 features", "6 features", ""]),
    ("Annotations + conservation (§8)", ["", "", "17 features", "17 features", ""]),
]
cols = ["Builds the\ntolerance\nscore (§3)", "Design A:\nnetworks predict\ntolerance", "Design A:\ntrees predict\ntolerance (14g)",
        "Design B:\nvariant model\n(14h, §16)", "Labels and\nevaluation"]
CW, LEFT = 1.3, 3.1
fig, ax = plt.subplots(figsize=(12, 7.4)); ax.axis("off")
ax.set_xlim(-LEFT, len(cols) * CW); ax.set_ylim(len(rows) + 2.2, -1.6)
for j, c in enumerate(cols):
    ax.text(j * CW + CW / 2, -0.75, c, ha="center", va="center", fontsize=8.8, weight="bold", color=INK)
for i, (name, cells) in enumerate(rows):
    if i % 2 == 0:
        ax.add_patch(plt.Rectangle((-LEFT, i), LEFT + len(cols) * CW, 1, color="#f7f6f3", lw=0))
    ax.text(-LEFT + 0.08, i + 0.5, name, va="center", fontsize=9, color=INK)
    for j, text in enumerate(cells):
        if not text:
            ax.text(j * CW + CW / 2, i + 0.5, "–", ha="center", va="center", color="#c9c8c3", fontsize=10)
            continue
        kind = ("filter" if text.startswith("filter") else "planned" if text.startswith(("planned", "tested only"))
                else "label" if text.startswith(("THE LABEL", "baseline")) else "used")
        c = {"used": BLUE, "filter": ORANGE, "planned": MUTED, "label": RED}[kind]
        ax.scatter(j * CW + 0.1, i + 0.5, s=50, color=c if kind != "planned" else "white", edgecolor=c, lw=1.5, zorder=3)
        ax.text(j * CW + 0.2, i + 0.5, text, va="center", fontsize=8.1, color=INK2)
legend = [(BLUE, "used as information (an input, or to build the score)"), (ORANGE, "used only to decide which variants count"),
          (RED, "the answer, or what accuracy is measured against"), (MUTED, "planned, or tried only as a test")]
for k, (c, lab) in enumerate(legend):
    lx, ly = -LEFT + 0.1 + (k % 2) * 4.3, len(rows) + 0.55 + (k // 2) * 0.6
    ax.scatter(lx, ly, s=45, color=c if c != MUTED else "white", edgecolor=c, lw=1.5, clip_on=False)
    ax.text(lx + 0.12, ly, lab, va="center", fontsize=8.4, color=INK2)
ax.text(-LEFT + 0.1, len(rows) + 1.85, "'–' = not used there. Note: gnomAD allele FREQUENCY was added to the variant table "
        "but is not used by any model.", fontsize=8.4, color=INK2, va="center", style="italic")
save(fig, "f02_data_use_map.png")

# =========================================================================== F3: building the tolerance score
from matplotlib.lines import Line2D  # noqa: E402
L = len(row.sequence); x = np.arange(1, L + 1)
gv = var[(var.gene == GENE)]
gn = gv[(gv.database == "GnomAD") & (gv.total_count >= 2)]
cv = gv[gv.database == "ClinVar"]
gn_cnt = np.bincount(gn.position - 1, minlength=L); cv_cnt = np.bincount(cv.position - 1, minlength=L)
norm = lambda a: (np.asarray(a) - np.min(a)) / (np.ptp(a) or 1)
dis = np.asarray(row.disorder_scores)
fig, axes = plt.subplots(3, 1, figsize=(11, 7.6), sharex=True, gridspec_kw={"height_ratios": [1.1, 1, 1.15]})
ax = axes[0]
ax.vlines(x[gn_cnt > 0], 0, gn_cnt[gn_cnt > 0], color=BLUE, lw=1.1)
ax.vlines(x[cv_cnt > 0], 0, -cv_cnt[cv_cnt > 0], color=RED, lw=1.1)
ax.axhline(0, color=MUTED, lw=0.8); ax.set_ylim(-cv_cnt.max() - 1.5, gn_cnt.max() + 3.5)
ax.set_ylabel("variants\nper residue")
ax.legend(handles=[Line2D([], [], color=BLUE, lw=2, label=f"gnomAD (up): {len(gn)} variants seen in ≥ 2 people of the healthy population"),
                   Line2D([], [], color=RED, lw=2, label=f"ClinVar (down): {len(cv)} variants reported to the clinical database, any classification")],
          frameon=True, framealpha=0.95, edgecolor="white", fontsize=8.4, loc="upper left", ncol=1)
ax.set_title(f"Step 1. Where the missense variants are in {GENE} ({L} residues)", loc="left")
ax = axes[1]
ax.plot(x, norm(row.gnomad_y_v4), color=BLUE, lw=1.6)
ax.plot(x, norm(row.clinvar_y_v4), color=RED, lw=1.6)
ax.set_ylabel("smoothed density\n(scaled 0 to 1)")
ax.set_title("Step 2. Each variant is weighted by how drastic its swap is, then spread over its neighbours "
             "(blue: gnomAD, red: ClinVar)", loc="left")
ax = axes[2]
t = np.asarray(row.tolerance_v4)
ax.fill_between(x, 0, 1, where=dis >= 0.5, transform=ax.get_xaxis_transform(), color=TINT[YELLOW], lw=0, step="mid")
ax.plot(x, t, color=INK, lw=1.4)
ax.fill_between(x, t, 0, where=t > 0, color=BLUE, alpha=0.25, interpolate=True)
ax.fill_between(x, t, 0, where=t < 0, color=RED, alpha=0.25, interpolate=True)
ax.axhline(0, color=MUTED, lw=0.8); ax.axhline(np.median(t), color=INK2, lw=1, ls="--")
ax.text(52, 0.46, f"dashed line: this gene's median ({np.median(t):+.2f}). Subtracting it 'centres' the target,\n"
        "so every gene's track is measured from its own typical level.", fontsize=8.3, color=INK2, va="top")
ax.annotate(f"residue {VAR_POS}: p.A{VAR_POS}P\nis a pathogenic variant", (VAR_POS, t[VAR_POS - 1]), (VAR_POS + 22, -0.5),
            fontsize=8.3, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8), color=INK2)
ax.set_ylim(-0.72, 0.5)
ax.set_ylabel("tolerance\n= blue − red"); ax.set_xlabel(f"{GENE} residue")
ax.set_title("Step 3. Tolerance = gnomAD density − ClinVar density. Above 0: tolerated; below 0: intolerant. "
             "Yellow bands: disordered", loc="left")
fig.tight_layout(); save(fig, "f03_tolerance_build.png")

# kernel explainer
sys.path.insert(0, os.path.abspath(os.path.join(OUT, "..", "..")))
from convolve import substitution_weight  # noqa: E402
fig, axes = plt.subplots(1, 2, figsize=(11, 3.0), gridspec_kw={"width_ratios": [1.2, 1]})
ax = axes[0]
for w, c, lab in [(31, BLUE, "ordered residue: w = 31\n(spread over 63 residues)"),
                  (5, ORANGE, "disordered residue: w = 5\n(spread over 11 residues)")]:
    k = binom.pmf(np.arange(2 * w + 1), 2 * w, 0.5)
    ax.plot(np.arange(-w, w + 1), k / k.max(), color=c, lw=2, label=lab)
ax.axvline(0, color=MUTED, lw=0.8, ls=":")
ax.set_xlabel("distance from the variant (residues)"); ax.set_ylabel("relative weight")
ax.legend(frameon=False, fontsize=8, loc="upper left")
ax.set_title("How one variant's weight is spread to its neighbours", loc="left")
ax = axes[1]
subs = [("A", "S"), ("A", "V"), ("A", "G"), ("A", "E"), ("A", "P"), ("G", "W"), ("C", "W")]
wb = [substitution_weight(a, b, "blosum62") for a, b in subs]
wg = [substitution_weight(a, b, "grantham") for a, b in subs]
yy = np.arange(len(subs))[::-1]
ax.barh(yy + 0.18, wb, height=0.34, color=BLUE, label="BLOSUM62 (used in ordered regions)")
ax.barh(yy - 0.18, wg, height=0.34, color=ORANGE, label="Grantham (used in disordered regions)")
ax.set_yticks(yy, [f"{a}→{b}" for a, b in subs]); ax.set_xlim(0, 1.05)
ax.set_xlabel("weight of the variant (0 = harmless swap, 1 = most drastic)")
ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.45, -0.3), ncol=2)
ax.set_title("How much one variant counts: the severity of the swap", loc="left")
fig.tight_layout(); save(fig, "f03b_kernel.png")

# =========================================================================== F4: ESM scores
import esm  # noqa: E402  (alphabet only, no weights)
alphabet = esm.data.Alphabet.from_architecture("ESM-1b")
AA = list("AVLIMFWCGPSTYNQHDEKR")
lp = np.load(os.path.join(MODEL_DIR, "psd_esm_scores.npz"))[GENE][:, [alphabet.get_idx(a) for a in AA]]
lp = lp - np.logaddexp.reduce(lp, axis=1, keepdims=True)                     # over the 20 amino acids
lo_, hi_ = VAR_POS - 20, VAR_POS + 19
win = np.arange(lo_, hi_ + 1)
seq = row.sequence
fig = plt.figure(figsize=(11, 7.0))
ax = fig.add_axes([0.07, 0.40, 0.83, 0.52])
im = ax.imshow(np.clip(lp[win - 1].T, -14, 0), aspect="auto", cmap="Blues", vmin=-14, vmax=0,
               extent=[lo_ - 0.5, hi_ + 0.5, len(AA) - 0.5, -0.5])
for p in win:
    ax.add_patch(plt.Rectangle((p - 0.5, AA.index(seq[p - 1]) - 0.5), 1, 1, fill=False, ec=RED, lw=1.3))
ax.add_patch(plt.Rectangle((VAR_POS - 0.5, AA.index("P") - 0.5), 1, 1, fill=False, ec=ORANGE, lw=2.6))
ax.set_yticks(range(len(AA)), AA, fontsize=8)
ax.set_xticks(win[::5], [f"{p}\n{seq[p - 1]}" for p in win[::5]], fontsize=8)
ax.set_ylabel("each amino acid ESM-1v could place here")
ax.set_title(f"ESM-1v's log-probability for all 20 amino acids at {GENE} residues {lo_}–{hi_}. Red outline: the amino acid "
             f"that is there; orange: proline at {VAR_POS}", loc="left", fontsize=10)
cax = fig.add_axes([0.915, 0.40, 0.015, 0.52])
cb = fig.colorbar(im, cax=cax); cb.set_label("log-probability (0 = certain, −14 = essentially never)", fontsize=8.3)
lp_wt, lp_p = lp[VAR_POS - 1, AA.index("A")], lp[VAR_POS - 1, AA.index("P")]
wt_lp = lp[win - 1, [AA.index(seq[p - 1]) for p in win]]
ent = -(np.exp(lp[win - 1]) * lp[win - 1]).sum(1)
mean_llr = (lp[win - 1].sum(1) - wt_lp) / 19 - wt_lp
tax = fig.add_axes([0.0, 0.03, 0.27, 0.27]); tax.axis("off")
tax.text(0.0, 1.0, f"The variant score (LLR) for p.A{VAR_POS}P\n\n"
         f"LLR = log P(proline) − log P(alanine)\n      = {lp_p:.1f} − ({lp_wt:.1f}) = {lp_p - lp_wt:.1f}\n\n"
         "Very negative: ESM finds this swap\nhighly unexpected. 0 would mean\n'as expected as what is there'.",
         fontsize=8.8, va="top", color=INK, bbox=dict(fc=TINT[ORANGE], ec="#c9c8c3", lw=0.6, boxstyle="round,pad=0.5"))
for k, (vals, lab, c) in enumerate([(wt_lp, "esm_wt_lp: log-probability of\nthe amino acid that is there", VIOLET),
                                     (ent, "esm_entropy: how spread out\nESM's guesses are", AQUA),
                                     (mean_llr, "esm_mean_llr: average LLR\nof all 19 possible swaps", RED)]):
    a2 = fig.add_axes([0.33 + k * 0.23, 0.05, 0.19, 0.22])
    a2.plot(win, vals, color=c, lw=1.5, marker="o", ms=2.3)
    a2.axvline(VAR_POS, color=ORANGE, lw=1, ls=":")
    a2.set_title(lab, loc="left", fontsize=8.6); a2.tick_params(labelsize=7.5)
save(fig, "f04_esm_scores.png")

# =========================================================================== F5: ESM, three ways
fig, ax = canvas(11.5, 4.8)
box(ax, 0.0, 0.33, 0.16, 0.38, "ESM protein\nlanguage model",
    "learned from tens of\nmillions of natural\nprotein sequences by\nfilling in hidden\namino acids", VIOLET, body_size=8.4)
ax.text(0.0, 0.97, "ESM-1v (650 million weights) was used for 1 and 2; ESM-2 is the model planned for 3.\n"
        "ESM's own weights were never changed in any model trained so far.", fontsize=8.6, color=INK2, va="top")
ways = [
    (0.62, "1. Scores (log-probabilities)",
     "§4: one pass per protein. Per residue: esm_wt_lp,\nesm_entropy, esm_mean_llr. Per variant: esm_llr\n(the LLR of that exact swap).",
     "4 features of the final variant model (§16);\n3 features in 14g; and the 'ESM-1v alone'\nbaseline every result is compared with.", GREEN),
    (0.33, "2. Representations (1,280 numbers per residue)",
     "14c: ESM's internal description of each residue in\nits context, read out with ESM frozen (weights\nunchanged); cached in probe_reps/.",
     "the input of the linear probe and the network\nheads that predict tolerance (14c, 14f);\ncompressed to 64 numbers for 14g.", GREEN),
    (0.04, "3. Weights, fine-tuned (LoRA)",
     "§12: ESM-2 650M with small trainable add-on\nmatrices in its attention layers, so the whole\nmodel can adapt to the tolerance target.",
     "planned: the Biowulf run (§15), Design A.\nNot run yet; the §13 smoke test used a tiny\nESM-2 (8 million weights) on this laptop.", "grey"),
]
for y, t_, how, use, colr in ways:
    box(ax, 0.22, y, 0.37, 0.25, t_, how, "grey", body_size=8.4)
    box(ax, 0.63, y, 0.37, 0.25, "Where it went", use, colr, body_size=8.4)
    arrow(ax, 0.16, 0.52, 0.22, y + 0.125, VIOLET); arrow(ax, 0.59, y + 0.125, 0.63, y + 0.125)
save(fig, "f05_esm_three_ways.png")

# =========================================================================== F6: residue table -> variant row
rg = res[(res.gene == GENE) & res.pos.between(VAR_POS - 3, VAR_POS + 3)].set_index("pos")
cols6 = [("tolerance_v4_centered", "tolerance\n(centred)", BLUE, "{:+.2f}"), ("esm_mean_llr", "ESM mean\nLLR", VIOLET, "{:.1f}"),
         ("esm_entropy", "ESM\nentropy", VIOLET, "{:.2f}"), ("rsa", "RSA", AQUA, "{:.2f}"), ("plddt", "pLDDT", AQUA, "{:.0f}"),
         ("iface_any", "interface", AQUA, "{}"), ("ann_transmembrane", "trans-\nmembrane", YELLOW, "{}"),
         ("disorder", "disorder", YELLOW, "{:.2f}"), ("msa_conservation", "conser-\nvation", YELLOW, "{:.1f}")]
fig, ax = canvas(11.5, 5.2)
ax.text(0, 0.985, f"One row per residue (residue_df): every source lined up on the same numbering. {GENE} residues "
        f"{VAR_POS - 3}–{VAR_POS + 3}", fontsize=10.5, weight="bold", va="top")
x0, cw, rh, top = 0.0, 0.0885, 0.062, 0.89
heads = [("pos", "residue", "grey"), ("aa", "amino\nacid", "grey")] + [(c, h, colr) for c, h, colr, _ in cols6]
for j, (_, h, colr) in enumerate(heads):
    ax.add_patch(plt.Rectangle((x0 + j * cw, top - rh * 1.3), cw - 0.004, rh * 1.3, color=TINT.get(colr, TINT["grey"])))
    ax.text(x0 + j * cw + cw / 2, top - rh * 0.65, h, ha="center", va="center", fontsize=8.3, weight="bold")
for i, p in enumerate(rg.index):
    y = top - rh * 1.3 - (i + 1) * rh
    hl = p == VAR_POS
    if hl:
        ax.add_patch(plt.Rectangle((x0, y), cw * len(heads) - 0.004, rh, color="#fff1d6", lw=0))
    vals = [str(p), rg.loc[p, "aa"]] + [("–" if pd.isna(rg.loc[p, c]) else
                                       (("yes" if rg.loc[p, c] else "no") if fmt == "{}" else fmt.format(rg.loc[p, c])))
                                      for c, _, _, fmt in cols6]
    for j, v in enumerate(vals):
        ax.text(x0 + j * cw + cw / 2, y + rh / 2, v, ha="center", va="center", fontsize=8.8,
                weight="bold" if hl else "normal", color=INK)
src_y = top - rh * 1.3 - 8.3 * rh
leg = [("grey", "sequence"), (BLUE, "human variants (§3)"), (VIOLET, "ESM-1v (§4)"), (AQUA, "structure (§7)"),
       (YELLOW, "annotations / disorder / conservation (§8)")]
for k, (colr, lab) in enumerate(leg):
    ax.add_patch(plt.Rectangle((k * 0.19, src_y), 0.018, 0.035, color=TINT.get(colr, TINT["grey"]), ec="#9b9a96", lw=0.5))
    ax.text(k * 0.19 + 0.025, src_y + 0.017, lab, va="center", fontsize=8.3, color=INK2)
vrow = var[(var.gene == GENE) & (var.protein_consequence == f"p.A{VAR_POS}P") & (var.database == "ClinVar")].iloc[0]
box(ax, 0.0, 0.02, 0.975, 0.2, f"One row per variant: p.A{VAR_POS}P (ClinVar: {vrow.clinsig_group}, {int(vrow.stars)} review star)",
    f"=  the whole residue-{VAR_POS} row above, looked up by (gene, position)\n"
    f"+  what depends on the exact swap: ESM LLR {vrow.esm_llr:.1f}, BLOSUM62 and Grantham severity, change in hydrophobicity\n"
    f"+  the human-variant features recomputed with residue {VAR_POS}'s own variants hidden: tolerance, counts within ±5 / ±15, gene-wide rates",
    "grey", body_size=8.6)
arrow(ax, 0.5, src_y - 0.015, 0.5, 0.22, ORANGE, lw=1.8)
save(fig, "f06_residue_to_variant.png")

# =========================================================================== F7: training ideas
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), gridspec_kw={"width_ratios": [1, 1.15]})
ax = axes[0]; ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
steps = [(0.5, 0.86, "1. Predict", "with the current weights"), (0.82, 0.5, "2. Compare", "with the answer:\nthe 'loss'"),
         (0.5, 0.14, "3. Adjust", "nudge every weight a\nlittle to shrink the loss"), (0.18, 0.5, "4. Repeat", "one pass over all the\ndata = one epoch")]
for x_, y_, t_, b_ in steps:
    ax.add_patch(FancyBboxPatch((x_ - 0.15, y_ - 0.11), 0.3, 0.22, boxstyle="round,pad=0.01,rounding_size=0.03",
                                fc=TINT[GREEN], ec="#9b9a96", lw=0.8))
    ax.text(x_, y_ + 0.05, t_, ha="center", fontsize=10, weight="bold")
    ax.text(x_, y_ - 0.02, b_, ha="center", fontsize=8.3, color=INK2, va="top")
for (a, b), (sx, sy, ex, ey) in zip([(0, 1), (1, 2), (2, 3), (3, 0)],
                                    [(0.65, 0.86, 0.82, 0.62), (0.82, 0.38, 0.65, 0.14), (0.35, 0.14, 0.18, 0.38), (0.18, 0.62, 0.35, 0.86)]):
    ax.annotate("", (ex, ey), (sx, sy), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.3, connectionstyle="arc3,rad=-0.3"))
ax.set_title("What 'training' means (for every model in this notebook)", loc="left")
ax = axes[1]
rng = np.random.default_rng(1)
n_cl = 24; fold = np.arange(n_cl) % 5; rng.shuffle(fold)
sizes = rng.integers(1, 5, n_cl)
for k in range(5):
    left = 0
    for c in range(n_cl):
        held = fold[c] == k
        ax.barh(4 - k, sizes[c], left=left, color=RED if held else BLUE, alpha=0.85 if held else 0.35, edgecolor="white", lw=1)
        left += sizes[c]
ax.set_yticks(range(5)[::-1], [f"round {k + 1}" for k in range(5)]); ax.set_xticks([])
ax.set_xlabel("gene clusters (bar width = number of genes in the cluster)")
ax.set_title("Cross-validation by cluster: fit on the blue clusters, score the red, rotate.\n"
             "Every cluster is scored once, by a model that never saw it or its relatives.", loc="left", fontsize=9.5)
ax.spines[["left", "bottom"]].set_visible(False)
fig.tight_layout(); save(fig, "f07_training_ideas.png")

# =========================================================================== F8: Spearman toy
rng = np.random.default_rng(4)
xs = np.arange(1, 41); target = np.sin(xs / 5) + rng.normal(0, 0.35, 40)
good = np.sin(xs / 5) * 0.4 + rng.normal(0, 0.12, 40); bad = rng.normal(0, 0.3, 40)
fig, axes = plt.subplots(1, 3, figsize=(11.5, 2.9), gridspec_kw={"width_ratios": [1.6, 1, 1]})
ax = axes[0]
ax.plot(xs, target, color=INK, lw=1.6, label="target (observed tolerance)")
ax.plot(xs, good, color=BLUE, lw=1.4, label=f"prediction A (Spearman {spearmanr(good, target).statistic:.2f})")
ax.plot(xs, bad, color=ORANGE, lw=1.2, label=f"prediction B (Spearman {spearmanr(bad, target).statistic:.2f})")
ax.set_xlabel("residue"); ax.legend(frameon=False, fontsize=7.8, loc="lower left", ncol=1)
ax.set_title("Two predictions for one gene", loc="left")
for ax, p_, c, nm in [(axes[1], good, BLUE, "A"), (axes[2], bad, ORANGE, "B")]:
    ax.scatter(pd.Series(target).rank(), pd.Series(p_).rank(), s=14, color=c)
    ax.set_xlabel("rank of the target"); ax.set_ylabel(f"rank of prediction {nm}")
    ax.set_title(f"Ranks agree? Spearman = {spearmanr(p_, target).statistic:.2f}", loc="left")
fig.tight_layout(); save(fig, "f08_spearman.png")

# =========================================================================== F9: AUROC toy
scores = np.array([-14, -12.5, -11, -9.5, -8.6, -7.2, -6.0, -5.1, -4.4, -3.3, -2.6, -1.8])
benign = np.array([0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1], bool)
fig, axes = plt.subplots(1, 2, figsize=(11.5, 3.2), gridspec_kw={"width_ratios": [1.5, 1]})
ax = axes[0]
for s, b in zip(scores, benign):
    ax.scatter(s, 0, s=120, color=BLUE if b else RED, zorder=3, edgecolor="white", lw=1.5)
ax.axvline(-6.6, color=INK2, ls="--", lw=1)
ax.text(-6.5, 0.32, "a threshold: call everything to the right 'benign'", fontsize=8.3, color=INK2)
ax.set_yticks([]); ax.set_ylim(-0.6, 0.6); ax.set_xlabel("score (higher = model thinks more benign)")
ax.spines["left"].set_visible(False)
ax.scatter([], [], color=BLUE, s=60, label="truly benign"); ax.scatter([], [], color=RED, s=60, label="truly pathogenic")
ax.legend(frameon=False, fontsize=8.5, loc="lower left", ncol=2)
pairs = [(b, p) for b in scores[benign] for p in scores[~benign]]
right = sum(b > p for b, p in pairs)
ax.set_title(f"12 variants sorted by score: {right} of {len(pairs)} benign–pathogenic pairs are in the right order", loc="left")
ax = axes[1]
fpr, tpr, _ = roc_curve(benign, scores)
ax.plot(fpr, tpr, color=BLUE, lw=2, drawstyle="steps-post"); ax.fill_between(fpr, tpr, step="post", color=BLUE, alpha=0.12)
ax.plot([0, 1], [0, 1], color=MUTED, ls="--", lw=1)
ax.text(0.45, 0.25, f"AUROC = {roc_auc_score(benign, scores):.2f}\n(= {right}/{len(pairs)})", fontsize=10, color=INK)
ax.set_xlabel("share of pathogenic called benign"); ax.set_ylabel("share of benign called benign")
ax.set_title("Slide the threshold: the ROC curve", loc="left")
fig.tight_layout(); save(fig, "f09_auroc.png")

# =========================================================================== F10: pooled vs within genes
g1 = (np.array([-9, -8.5, -8, -7.2, -6.5, -6]), np.array([-5.5, -5, -4.6, -4.2, -3.8, -3.5]))   # gene 1: (pathogenic, benign)
g2 = (np.array([-4.8, -4.2, -3.9, -3.4, -3.1, -2.8]), np.array([-2.2, -1.9, -1.5, -1.2, -0.9, -0.6]))
fig, ax = plt.subplots(figsize=(11.5, 2.8))
for yrow, (nm, (p_, b_)) in enumerate([("gene 1", g1), ("gene 2", g2)]):
    ax.scatter(p_, np.full(6, 2 - yrow), s=90, color=RED, edgecolor="white", lw=1.2, zorder=3)
    ax.scatter(b_, np.full(6, 2 - yrow), s=90, color=BLUE, edgecolor="white", lw=1.2, zorder=3)
    yy_ = np.r_[np.zeros(6), np.ones(6)]
    ax.text(0.2, 2 - yrow, f"{nm}: AUROC {roc_auc_score(yy_, np.r_[p_, b_]):.2f}", va="center", fontsize=9)
allp, allb = np.r_[g1[0], g2[0]], np.r_[g1[1], g2[1]]
ax.scatter(allp, np.zeros(12), s=90, color=RED, edgecolor="white", lw=1.2, zorder=3)
ax.scatter(allb, np.zeros(12), s=90, color=BLUE, edgecolor="white", lw=1.2, zorder=3)
ax.text(0.2, 0, f"pooled: AUROC {roc_auc_score(np.r_[np.zeros(12), np.ones(12)], np.r_[allp, allb]):.2f}", va="center", fontsize=9, weight="bold")
ax.set_yticks([2, 1, 0], ["gene 1", "gene 2", "both genes\ntogether"]); ax.set_xlim(-9.6, 3.6)
ax.set_xlabel("score (higher = more benign)"); ax.spines["left"].set_visible(False)
ax.set_title("Perfect within each gene, worse when pooled: gene 2's whole scale is shifted (like RYR1 vs ABL1 for ESM-1v)", loc="left")
fig.tight_layout(); save(fig, "f10_pooled_vs_within.png")

# =========================================================================== F11: the bootstrap (real: Q3)
train_genes = {g for g, s in split_of.items() if s in ("train", "val")}
tr = var[(var.database == "ClinVar") & var.clinsig_group.isin(["pathogenic", "benign"]) & (var.stars >= 1)
         & var.gene.isin(train_genes)].drop_duplicates("id")
unsure_llr = tr.loc[tr.clinsig_group == "benign", "esm_llr"].median()
yv = pred.clinsig_group.eq("benign").to_numpy(); keep = (pred["ESM-1v LLR"] > unsure_llr).to_numpy()
y_, g_ = yv[keep], pred.gene.to_numpy()[keep]
sa, sb = pred["logistic: all features"].to_numpy()[keep], pred["ESM-1v LLR"].to_numpy()[keep]
rows_of = pd.Series(np.arange(len(g_))).groupby(g_).indices; keys = np.array(list(rows_of)); rng = np.random.default_rng(0)
diffs = []
for _ in range(1000):
    r_ = np.concatenate([rows_of[k] for k in rng.choice(keys, len(keys))])
    if len(set(y_[r_])) == 2:
        diffs.append(roc_auc_score(y_[r_], sa[r_]) - roc_auc_score(y_[r_], sb[r_]))
lo95, hi95 = np.percentile(diffs, [2.5, 97.5]); obs = roc_auc_score(y_, sa) - roc_auc_score(y_, sb)
print(f"Q3 reproduced: {obs:.3f} ({lo95:.3f} to {hi95:.3f}); results.json: "
      f"{results['questions']['Q3: primary beats ESM-1v where ESM is unsure']}")
fig, ax = plt.subplots(figsize=(11.5, 2.9))
ax.hist(diffs, bins=40, color=BLUE, alpha=0.55)
ax.axvline(0, color=INK, lw=1.2); ax.axvline(obs, color=BLUE, lw=2)
ax.axvspan(lo95, hi95, color=GREEN, alpha=0.08)
for v_, t_ in [(lo95, f"2.5%: {lo95:+.2f}"), (hi95, f"97.5%: {hi95:+.2f}")]:
    ax.axvline(v_, color=GREEN, ls="--", lw=1.2); ax.text(v_, ax.get_ylim()[1] * 0.92, " " + t_, fontsize=8.5, color=GREEN)
ax.text(obs, ax.get_ylim()[1] * 0.6, f" observed {obs:+.3f}", fontsize=9, color=BLUE)
ax.text(0.006, ax.get_ylim()[1] * 0.8, "0 = no difference", fontsize=8.5, color=INK, ha="left")
ax.set_xlabel("AUROC(model) − AUROC(ESM-1v) where ESM is unsure, in one resampled test set")
ax.set_ylabel("resamples")
ax.set_title("Q3, resampled 1,000 times by whole genes: the whole 95% range (green) is above 0, so the question passes", loc="left")
fig.tight_layout(); save(fig, "f11_bootstrap.png")

# =========================================================================== F12: what improved
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), gridspec_kw={"width_ratios": [1, 1.1]})
ax = axes[0]
items = [("ESM-1v zero-shot (no training)", 0.020), ("disorder alone (no training)", 0.052), ("best from-scratch network (BiLSTM)", 0.054),
         ("ESM-1v + linear probe (14c)", 0.105), ("ESM-1v + BiLSTM head (14f, best)", 0.125), ("trees + structure/annotations (14g)", 0.126)]
yy = np.arange(len(items))[::-1]
ax.hlines(yy, 0, [v for _, v in items], color=GRID, lw=2)
ax.scatter([v for _, v in items], yy, s=60, color=[MUTED, MUTED, ORANGE, BLUE, BLUE, AQUA], zorder=3)
for y_i, (_, v) in zip(yy, items):
    ax.text(v + 0.008, y_i, f"{v:.3f}", va="center", fontsize=8.5)
ax.axvline(0.39, color=RED, ls="--", lw=1.2); ax.text(0.397, 4.6, "noise ceiling 0.39\n(the most any model\ncould reach)", color=RED, fontsize=8.3, va="top")
ax.scatter([0.320], [-1.3], s=60, color=VIOLET, zorder=3); ax.text(0.33, -1.3, "0.320, with its own ceiling 0.61", va="center", fontsize=8.5)
ax.set_yticks(list(yy) + [-1.3], [n for n, _ in items] + ["linear probe, pathogenic-only target (14g)"], fontsize=8.5)
ax.set_xlim(0, 0.66); ax.set_xlabel("within-gene Spearman with the tolerance target (validation genes)")
ax.set_title("Design A: predicting the tolerance score", loc="left")
ax = axes[1]
rowsB = [("all variants", "14h cross-validation", 0.906, 0.961), ("within genes", "14h cross-validation", 0.885, 0.920),
         ("where ESM is unsure", "14h cross-validation", 0.532, 0.867), ("all variants", "§16 test genes", 0.620, 0.890),
         ("within genes", "§16 test genes", 0.749, 0.904), ("where ESM is unsure", "§16 test genes", 0.545, 0.879)]
yy = np.array([6.5, 5.5, 4.5, 2.5, 1.5, 0.5])
for y_i, (m, src, a, b) in zip(yy, rowsB):
    ax.plot([a, b], [y_i, y_i], color=GRID, lw=3, zorder=1)
    ax.scatter(a, y_i, s=55, color=MUTED, zorder=3); ax.scatter(b, y_i, s=70, color=BLUE, zorder=3)
    ax.text(b + 0.012, y_i, f"{a:.3f} → {b:.3f}", va="center", fontsize=8.3)
ax.set_yticks(yy, [m for m, *_ in rowsB], fontsize=8.5)
ax.text(0.505, 7.25, "14h: cross-validation, training + validation genes", fontsize=8.8, weight="bold")
ax.text(0.505, 3.25, "§16: the untouched test genes", fontsize=8.8, weight="bold")
ax.scatter([], [], color=MUTED, s=45, label="ESM-1v alone"); ax.scatter([], [], color=BLUE, s=60, label="final model (all 42 features)")
ax.legend(frameon=False, fontsize=8.3, loc="lower left", bbox_to_anchor=(0, -0.28), ncol=2)
ax.set_xlim(0.5, 1.07); ax.set_ylim(0, 7.6); ax.set_xlabel("AUROC, benign vs pathogenic ClinVar variants")
ax.set_title("Design B: predicting the ClinVar label", loc="left")
fig.tight_layout(); save(fig, "f12_improvements.png")
print("done")

# =========================================================================== F13: the test genes, a closer look (2 x 2)
from sklearn.metrics import confusion_matrix, brier_score_loss  # noqa: E402
P = "logistic: all features"
tp_ = pred.assign(benign=pred.clinsig_group.eq("benign"))
y_t, p_t = tp_.benign.to_numpy(), tp_[P].to_numpy()
tr_lab = tr.assign(benign=tr.clinsig_group.eq("benign"))
big3 = ["RYR1", "ABL1", "VHL"]
fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2))
ax = axes[0, 0]
STAR = ["#b7d3f6", "#5598e7", "#1c5cab"]
for k, (name, t_) in enumerate([("training", tr_lab), ("test", tp_)]):
    for j, b in enumerate([False, True]):
        share = t_[t_.benign == b].stars.value_counts(normalize=True).reindex([1, 2, 3], fill_value=0)
        left = 0
        for s_, c in zip([1, 2, 3], STAR):
            ax.barh(k * 2.4 + j, share[s_], left=left, color=c, edgecolor="white", lw=1.5,
                    label=f"{s_} review star{'s' if s_ > 1 else ''}" if (k, j) == (0, 0) else None)
            if share[s_] > 0.06:
                ax.text(left + share[s_] / 2, k * 2.4 + j, f"{share[s_]:.0%}", ha="center", va="center", fontsize=8,
                        color="white" if s_ > 1 else INK)
            left += share[s_]
ax.set_yticks([0, 1, 2.4, 3.4], ["training: pathogenic", "training: benign", "test: pathogenic", "test: benign"])
ax.set_xlabel("share of variants"); ax.legend(frameon=False, fontsize=8.3, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.15))
ax.set_title("How confident the labels are (ClinVar review stars)", loc="left")
ax = axes[0, 1]
per = tp_.groupby("gene").benign.agg(["size", "sum"]).sort_values("size", ascending=False).head(10)[::-1]
ax.barh(range(len(per)), per["size"] - per["sum"], color=RED, label="pathogenic")
ax.barh(range(len(per)), per["sum"], left=per["size"] - per["sum"], color=BLUE, label="benign")
ax.set_yticks(range(len(per)), per.index); ax.set_xlabel("test variants")
ax.legend(frameon=False, fontsize=8.5, loc="lower right")
ax.set_title(f"The 10 largest test genes: RYR1, ABL1 and VHL hold {tp_.gene.isin(big3).mean():.0%} of all test variants",
             loc="left", fontsize=10)
ax = axes[1, 0]
cm = confusion_matrix(y_t, p_t >= 0.5)
ax.imshow(cm / cm.sum(axis=1, keepdims=True), cmap="Blues", vmin=0, vmax=1)
for i in range(2):
    for j in range(2):
        frac = cm[i, j] / cm[i].sum()
        ax.text(j, i, f"{cm[i, j]} variants\n({frac:.0%} of the row)", ha="center", va="center", fontsize=9.5,
                color="white" if frac > 0.6 else INK)
ax.set_xticks([0, 1], ["called pathogenic\n(P(benign) < 0.5)", "called benign\n(P(benign) ≥ 0.5)"])
ax.set_yticks([0, 1], ["truly\npathogenic", "truly\nbenign"])
ax.spines[["left", "bottom"]].set_visible(False)
ax.set_title(f"Decisions at a cut-off of 0.5: accuracy {(cm[0, 0] + cm[1, 1]) / cm.sum():.0%}", loc="left")
ax = axes[1, 1]
ax.plot([0, 1], [0, 1], color=MUTED, ls="--", lw=1)
for name, m, c in [("RYR1", tp_.gene == "RYR1", RED), ("ABL1", tp_.gene == "ABL1", BLUE), ("VHL", tp_.gene == "VHL", ORANGE),
                   ("other 82 genes", ~tp_.gene.isin(big3), INK2)]:
    m = m.to_numpy()
    ax.scatter(p_t[m].mean(), y_t[m].mean(), s=m.sum() / 2 + 30, color=c, alpha=0.85, zorder=3, edgecolor="white", lw=1.5)
    right_side = name == "VHL"
    ax.annotate(f"{name} ({m.sum()} variants)", (p_t[m].mean(), y_t[m].mean()), fontsize=8.5,
                xytext=(12, 6) if right_side else (-10, 10), textcoords="offset points", ha="left" if right_side else "right")
ax.set_xlim(0, 1.05); ax.set_ylim(0, 1.05)
ax.set_xlabel("average predicted P(benign)"); ax.set_ylabel("share that really is benign")
ax.set_title(f"Calibration by gene group (on the dashed line = honest). Brier score {brier_score_loss(y_t, p_t):.3f}",
             loc="left", fontsize=10)
fig.tight_layout(); save(fig, "f13_test_closer_look.png")
