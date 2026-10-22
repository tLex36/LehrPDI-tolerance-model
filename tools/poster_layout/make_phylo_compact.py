"""Poster-sized version of I0 (conservation across the species tree), drawn at its exact slot size so the
text stays readable: same data and content as notebook 06 §3, compact layout.

    python make_phylo_compact.py [width_in height_in]     # default 7.5 x 5.85
Writes poster/svg/I0_conservation_phylogeny_poster.svg and a 300-dpi PNG next to this script (assets/).
"""
import os, re, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/Users/awlehr/LehrPDI/figures")
import evolution_conservation as ec
from species_tree import build_tree, draw_species_tree, label_nodes

W, H = (float(sys.argv[1]), float(sys.argv[2])) if len(sys.argv) > 2 else (7.5, 5.85)
HERE = os.path.dirname(os.path.abspath(__file__))
POSTER = "/Users/awlehr/LehrPDI/figures/results/tolerance_model/poster"
CACHE = os.path.join(POSTER, "conservation_cache")
os.makedirs(os.path.join(HERE, "assets"), exist_ok=True)
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
plt.rcParams.update({"font.family": "Arial", "font.size": 10, "svg.fonttype": "none", "figure.facecolor": "white",
                     "axes.edgecolor": "#8a8984", "xtick.color": "#52514e", "xtick.labelsize": 8.5})

mane = pd.read_csv("/Users/awlehr/LehrPDI/mane_transcripts.csv")
mane = mane[mane.MANE_status == "MANE Select"]
ensg = dict(zip(mane.symbol, mane.Ensembl_Gene))
psd = pd.read_parquet("/Users/awlehr/LehrPDI/figures/results/tolerance_model/tolerance_v4.parquet", columns=["gene"]).gene.unique()
GENE_SETS = {"Postsynaptic\nproteins": [ensg[g] for g in psd if g in ensg],
             "Olfactory\nreceptors": [ensg[s] for s in mane.symbol if re.fullmatch(r"OR\d+[A-Z]+\d+[A-Z]?", str(s))],
             "Ribosomal\nproteins": [ensg[s] for s in mane.symbol if re.fullmatch(r"RP[LS]\d+[AXYL]?\d*|RPLP[0-2]|RPSA", str(s))]}
REGION_SETS = {"Regulatory\nelements": lambda: ec.sample_regions("regulatory", 80, os.path.join(CACHE, "regulatory_regions.json")),
               "Coding\nregions": lambda: ec.sample_cds(mane.Ensembl_nuc.tolist(), 80, os.path.join(CACHE, "cds_regions.json"))}

value, coverage = {}, {}
for col, genes in GENE_SETS.items():
    per_gene = ec.protein_similarity(genes, os.path.join(CACHE, "similarity"))
    value[col] = ec.summarise(per_gene.where(per_gene > 0).median().to_frame().T)
    coverage[col] = ec.summarise(((per_gene > 0).mean() * 100).to_frame().T)
aligned = ec.species_with_alignment()
for col, sample in REGION_SETS.items():
    value[col] = ec.summarise(ec.dna_identity(sample(), os.path.join(CACHE, "dna"), aligned))

tree = build_tree(keep={"Human", *ec.SPECIES}); tree.ladderize(reverse=True)
tips = [t.name for t in tree.get_terminals()]
table = pd.DataFrame(value).reindex(tips); table.loc["Human"] = 100.0
cov = pd.DataFrame(coverage).reindex(tips); cov.loc["Human"] = 100.0
n, m = table.shape
n_protein = len(GENE_SETS)

fig = plt.figure(figsize=(W, H))
gs = fig.add_gridspec(2, 2, width_ratios=[1.16, 1.84], height_ratios=[n, 0.55], wspace=0.02, hspace=0.06,
                      left=0.01, right=0.975, top=0.86, bottom=0.07)
ax_tree, ax_hm, ax_cb = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])
with plt.rc_context({"font.size": 10}):
    draw_species_tree(tree, ax=ax_tree, title="")
for t_ in ax_tree.texts:                                                # species names
    t_.set_fontsize(10)
label_nodes(tree, ax_tree, names=["Mammalia", "Amniota", "Vertebrata", "Ecdysozoa"], fontsize=8.5)
x0_, x1_ = ax_tree.get_xlim(); ax_tree.set_xlim(x0_ - 0.3 * (x1_ - x0_), x1_ * 1.8)   # room for 'Ecdysozoa'
ax_tree.xaxis.label.set_fontsize(9); ax_tree.tick_params(axis="x", labelsize=8.5, pad=1)
ax_tree.set_xlabel("Million years ago", fontsize=9, labelpad=1)
cmap = plt.get_cmap("Blues").copy(); cmap.set_bad("white")
shown = table.copy()
for col in GENE_SETS:
    shown.loc[cov[col] < 10, col] = np.nan
im = ax_hm.imshow(shown.to_numpy(float), cmap=cmap, vmin=0, vmax=100, aspect="auto", extent=(-0.5, m - 0.5, n + 0.5, 0.5))
for i, sp in enumerate(tips, start=1):
    for j, col in enumerate(table.columns):
        val = table.loc[sp, col]
        cv = cov.loc[sp, col] if j < n_protein else np.nan
        hatched = (j < n_protein and cv < 10) or np.isnan(val)
        if hatched:
            ax_hm.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, hatch="////", ec="#c9c8c3", lw=0))
            txt = ("none" if cv == 0 else f"{cv:.0f}% have one") if j < n_protein and not np.isnan(cv) else "no alignment"
            ax_hm.text(j, i, txt, ha="center", va="center", fontsize=6.5, color=MUTED)
        else:
            ink = "white" if val > 60 else INK
            if j < n_protein and sp != "Human":
                ax_hm.text(j - 0.08, i, f"{val:.0f}", ha="right", va="center", fontsize=10, color=ink, weight="bold")
                ax_hm.text(j - 0.02, i + 0.02, f"{cv:.0f}%", ha="left", va="center", fontsize=6.5, color=ink)
            else:
                ax_hm.text(j, i, f"{val:.0f}", ha="center", va="center", fontsize=10, color=ink, weight="bold")
ax_hm.set_xticks(np.arange(-0.5, m, 1), minor=True); ax_hm.set_yticks(np.arange(0.5, n + 1, 1), minor=True)
ax_hm.grid(which="minor", color="white", lw=1.5); ax_hm.tick_params(which="both", length=0)
ax_hm.axvline(n_protein - 0.5, color="white", lw=5)
ax_hm.set_xticks(range(m), table.columns, fontsize=8.5, color=INK2, linespacing=1.05); ax_hm.xaxis.tick_top()
ax_hm.tick_params(axis="x", pad=1.5)
ax_hm.set_yticks([]); ax_hm.set_ylim(n + 0.5, 0.5); ax_tree.set_ylim(n + 0.5, 0.5)
for s_ in ax_hm.spines.values():
    s_.set_visible(False)
ax_hm.text((n_protein - 1) / 2, -1.55, "Protein: % similarity", ha="center", fontsize=9.5, weight="bold", color=INK)
ax_hm.text(n_protein + (m - n_protein - 1) / 2, -1.55, "DNA: % identical", ha="center", fontsize=9.5, weight="bold", color=INK)
cb = fig.colorbar(im, cax=ax_cb, orientation="horizontal")
cb.set_label("% conserved vs human  (small %: share of genes with an orthologue)", fontsize=8.5, labelpad=1)
cb.ax.tick_params(labelsize=8, pad=1, length=2); cb.outline.set_visible(False)
svg = os.path.join(POSTER, "svg", "I0_conservation_phylogeny_poster.svg")
fig.savefig(svg)                                                        # exact size: no bbox cropping
fig.savefig(os.path.join(HERE, "assets", "I0_conservation_phylogeny_poster.png"), dpi=300)
print("wrote", svg, f"({W} x {H} in)")
