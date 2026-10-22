"""Build the data files of the interactive explorer (docs/) from the saved results of notebook 05.

    python tools/explorer/build_data.py --results /path/to/LehrPDI/figures/results/tolerance_model

Writes, under docs/data/:
    genes/<GENE>.js   one file per gene: sequence, per-residue features, annotation ranges and the variants
    genes.js          the gene index (length, split, variant counts)
    test.js           the official predictions on the untouched test genes (final_test/predictions.parquet)
    summary.js        result tables (cross-validation, test, feature groups, structure classes, conservation)

Every file is plain JSON wrapped in a function call (TOL.gene({...}), TOL.index([...]), ...), so the page also works
when index.html is opened straight from disk (file://), where browsers block fetch().

Check: the explorer recomputes the tolerance score in the browser from the variants written here. This script first
rebuilds tolerance_v4 in Python from the same variants (tolerance.build_tolerance_columns) and stops unless it matches
the saved scores exactly, so the files hold exactly the variants that fed the official score.
"""
import argparse, json, os, re, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
ap = argparse.ArgumentParser()
ap.add_argument("--results", default=os.path.join(REPO, "figures", "results", "tolerance_model"),
                help="results/tolerance_model folder with the parquet files of notebook 05")
ap.add_argument("--out", default=os.path.join(REPO, "docs", "data"))
ap.add_argument("--skip-conservation", action="store_true", help="leave the conservation table out")
ap.add_argument("--mane", default=None, help="mane_transcripts.csv (default: three folders above --results)")
args = ap.parse_args()
RES, OUT = args.results, args.out
args.mane = args.mane or os.path.join(RES, "..", "..", "..", "mane_transcripts.csv")
sys.path.insert(0, os.path.join(REPO, "figures"))
from tolerance import build_tolerance_columns                      # noqa: E402  (same recipe as notebook 05 §3)
sys.path.insert(0, HERE)
from jsdata import write_js, enc                                    # noqa: E402
os.makedirs(os.path.join(OUT, "genes"), exist_ok=True)


def bare(pc):
    """'p.E154D' -> 'E154D'"""
    return pc[2:] if pc.startswith("p.") else pc


def runs(mask):
    """Boolean per-residue mask -> [[start, end], ...], 1-based inclusive."""
    m = np.asarray(mask, bool)
    if not m.any():
        return []
    d = np.diff(np.r_[0, m.astype(int), 0])
    starts, ends = np.where(d == 1)[0] + 1, np.where(d == -1)[0]
    return [[int(s), int(e)] for s, e in zip(starts, ends)]


# ------------------------------------------------------------------ inputs
seq_df = pd.read_parquet(os.path.join(RES, "tolerance_v4.parquet"))
res_df = pd.read_parquet(os.path.join(RES, "residue_scores.parquet"))
var = pd.read_parquet(os.path.join(RES, "variant_scores.parquet"),
                      columns=["gene", "id", "position", "protein_consequence", "database", "clinsig_group", "stars",
                               "total_count", "sex_chrom", "n_hemi", "n_hom"])
manifest = json.load(open(os.path.join(RES, "split_manifest.json")))
pred = pd.read_parquet(os.path.join(RES, "final_test", "predictions.parquet"))
results = json.load(open(os.path.join(RES, "final_test", "results.json")))
poster_numbers = json.load(open(os.path.join(RES, "poster", "poster_numbers.json")))

# ------------------------------------------------------------------ which variants fed the score (notebook 05 §3)
# gnomAD: allele count >= 2; on X / Y outside the pseudoautosomal regions, only if a man (hemizygous) or a homozygous
# woman carries it (X_MODE = "hemizygous", X_MIN_MALE_CARRIERS = 1). ClinVar: every missense entry, any classification.
is_gnomad, is_clinvar = var.database.eq("GnomAD"), var.database.eq("ClinVar")
sex_linked = var.sex_chrom.fillna("").ne("")
male_seen = (var.n_hemi.fillna(0) + var.n_hom.fillna(0)) >= 1
keep_min = lambda n: is_gnomad & (var.total_count.fillna(0) >= n) & (~sex_linked | male_seen)
used2, used1 = keep_min(2), keep_min(1)                           # official filter; the 14g "+ singletons" variant

check = build_tolerance_columns(seq_df[["gene", "sequence", "disorder_scores"]], var[used2 | is_clinvar],
                                name="chk", gnomad_min_count=0, verbose=False)
worst = max(float(np.max(np.abs(np.asarray(a) - np.asarray(b)))) if len(a) else 0.0
            for a, b in zip(check.tolerance_chk, seq_df.tolerance_v4))
print(f"tolerance rebuilt from the exported variants: max |difference| from tolerance_v4 = {worst:.2e}")
if worst > 1e-9:
    sys.exit("the exported variants do not reproduce tolerance_v4: stopping")

# ------------------------------------------------------------------ per-gene files
CLS = {"pathogenic": "p", "benign": "b", "vus": "v"}               # anything else (conflicting, risk factor, ...) -> o
ANN = {"domain": "ann_domain", "repeat": "ann_repeat", "tm": "ann_transmembrane", "signal": "ann_signal",
       "ligand": "ann_ligand_site", "active": "ann_active_site", "disulfide": "ann_disulfide", "ptm": "ann_ptm",
       "motif": "ann_motif", "motifPred": "ann_motif_predicted", "lip": "ann_lip", "lipPred": "ann_lip_predicted",
       "bindingRegion": "ann_binding_region", "coiledCoil": "ann_coiled_coil", "conserved": "ann_conserved_site",
       "disorderCurated": "ann_disorder_curated", "iface": "iface_any", "ifaceHomo": "iface_homo",
       "ifaceHetero": "iface_hetero"}
res_by_gene = dict(tuple(res_df.sort_values(["gene", "pos"]).groupby("gene")))
var["used"] = np.select([used2, used1], ["u", "s"], "x")
vg = dict(tuple(var.groupby("gene")))
lab_all = var[is_clinvar & var.clinsig_group.isin(["pathogenic", "benign"]) & (var.stars >= 1)].drop_duplicates("id")
split = manifest["split"]
pred_by_gene = dict(tuple(pred.groupby("gene")))

index = []
for row in seq_df.itertuples():
    g, seq = row.gene, row.sequence
    L = len(seq)
    dis = np.asarray(row.disorder_scores, float)
    dis2 = np.round(dis, 2)                                         # 2 decimals, but never across the 0.5 threshold
    flip = (dis >= 0.5) != (dis2 >= 0.5)
    dis2[flip] = np.where(dis[flip] >= 0.5, 0.5, 0.49)
    r = res_by_gene.get(g)
    if r is not None and (len(r) != L or not np.array_equal(r.pos.to_numpy(), np.arange(1, L + 1))):
        sys.exit(f"{g}: residue table does not match the sequence")
    col = (lambda c, scale: enc(r[c], scale)) if r is not None else (lambda c, scale: "," * (L - 1))
    d = vg.get(g, var.iloc[:0]).copy()
    gn = d[d.database.eq("GnomAD")].drop_duplicates("id").sort_values("position")
    cv = d[d.database.eq("ClinVar")].drop_duplicates(["id", "clinsig_group"]).sort_values("position")
    gene = {
        "gene": g, "split": split.get(g, "excluded"),
        "uniprot": None if r is None else r.uniprot_id.iloc[0],
        "seq": seq,
        # per-residue values as integer strings: disorder, binding, rsa x 100; plddt x 1; esm, cons x 10
        "dis": enc(dis2, 100),
        "bind": col("binding", 100), "rsa": col("rsa", 100), "plddt": col("plddt", 1),
        "esm": col("esm_mean_llr", 10), "cons": col("msa_conservation", 10),
        "ann": {} if r is None else {k: runs(r[c].fillna(False).to_numpy(bool)) for k, c in ANN.items()},
        # tokens: E154Du = ref, position, alt, use flag (u: used, s: singleton, x: excluded by the X rule / count 0)
        "gn": " ".join(f"{bare(pc)}{u}" for pc, u in zip(gn.protein_consequence, gn.used)),
        # tokens: E154Dv1 = ref, position, alt, class (p / b / v / o), review stars
        "cv": " ".join(f"{bare(pc)}{CLS.get(c, 'o')}{int(s) if pd.notna(s) else 0}"
                       for pc, c, s in zip(cv.protein_consequence, cv.clinsig_group, cv.stars)),
        "check": [round(float(np.sum(row.tolerance_v4)), 6), round(float(np.sum(np.abs(row.tolerance_v4))), 6)],
    }
    bad = [t for t in (gene["gn"] + " " + gene["cv"]).split() if not re.match(r"^[A-Z]\d+[A-Z*][usxpbvo]\d?$", t)]
    if bad:
        sys.exit(f"{g}: unexpected variant tokens {bad[:3]}")
    if g in pred_by_gene:                                            # official test predictions (§16)
        p = pred_by_gene[g].sort_values("position")
        gene["test"] = [[bare(pc), int(cl == "benign"), round(float(a), 4), round(float(b), 3), round(float(c), 4)]
                        for pc, cl, a, b, c in zip(p.protein_consequence, p.clinsig_group, p["logistic: all features"],
                                                   p["ESM-1v LLR"], p["tolerance (own residue hidden)"])]
    write_js(os.path.join(OUT, "genes", f"{g}.js"), "TOL.gene", gene)
    lab = lab_all[lab_all.gene.eq(g)]
    index.append([g, L, gene["split"], int((gn.used == "u").sum()), int(len(cv)),
                  int(lab.clinsig_group.eq("pathogenic").sum()), int(lab.clinsig_group.eq("benign").sum()),
                  int(cv.clinsig_group.eq("vus").sum()), bool(row.tolerance_ok), gene["uniprot"]])
write_js(os.path.join(OUT, "genes.js"), "TOL.index",
         {"columns": ["gene", "length", "split", "gnomad_used", "clinvar", "labelled_pathogenic", "labelled_benign",
                      "vus", "tolerance_ok", "uniprot"], "rows": index})
print(f"wrote {len(index)} gene files")

# ------------------------------------------------------------------ test predictions (§16)
t = pred.sort_values(["gene", "position"])
write_js(os.path.join(OUT, "test.js"), "TOL.test", {
    "columns": ["gene", "variant", "benign", "stars", "model", "esm", "tolerance", "esm_tolerance", "trees"],
    "rows": [[g, bare(pc), int(c == "benign"), int(s), round(float(a), 5), round(float(b), 4), round(float(tl), 5),
              round(float(et), 5), round(float(tr), 5)]
             for g, pc, c, s, a, b, tl, et, tr in zip(t.gene, t.protein_consequence, t.clinsig_group, t.stars,
                                                     t["logistic: all features"], t["ESM-1v LLR"],
                                                     t["tolerance (own residue hidden)"],
                                                     t["logistic: ESM LLR + tolerance"],
                                                     t["gradient-boosted trees: all features"])],
    "esm_unsure_threshold": poster_numbers["esm_unsure_threshold"],
    "results": results})

# ------------------------------------------------------------------ result tables shown on the page
splits = res_df.drop_duplicates("gene").split.value_counts().to_dict()
residues = res_df.split.value_counts().to_dict()
summary = {
    "split": {s: {"genes": int(splits.get(s, 0)), "residues": int(residues.get(s, 0))} for s in ["train", "val", "test"]},
    "variants": {"gnomad_total": int(is_gnomad.sum()), "gnomad_used": int(used2.sum()),
                 "clinvar_total": int(is_clinvar.sum()),
                 "labelled_trainval": 15339, "labelled_trainval_pathogenic": 5449, "labelled_trainval_benign": 9890,
                 "labelled_test": 1357, "labelled_test_pathogenic": 602, "labelled_test_benign": 755,
                 "test_genes_with_labels": results["variants"]["genes"]},
    # notebook 05, 14h: 5-fold cross-validation by gene cluster on the 15,339 training + validation variants
    "cv": [["ESM-1v LLR alone", "no training: evolution only", 0.906, 0.885, 0.532, 0.620],
           ["Tolerance score alone", "no training: human variants only", 0.837, 0.728, 0.747, 0.829],
           ["Logistic regression", "ESM LLR + tolerance", 0.933, 0.889, 0.717, 0.826],
           ["Logistic regression", "all 42 features", 0.961, 0.920, 0.867, 0.890],
           ["Gradient-boosted trees", "all 42 features", 0.958, 0.920, 0.854, 0.860],
           ["MLP neural network", "42 features + 64 ESM components", 0.959, 0.920, 0.862, None],
           ["Window CNN", "±50 residues: sequence, variant counts", 0.943, 0.905, 0.798, None],
           ["Window BiLSTM", "±50 residues: sequence, variant counts", 0.950, 0.913, 0.818, None]],
    "feature_groups": poster_numbers["feature_groups"],
    # held-out ClinVar pathogenic / benign variants (>= 1 star) by structural class, all genes (notebook 06 §8)
    "struct_classes": [["core", "core (RSA < 0.2)", 220, 327], ["boundary", "boundary (0.2 to 0.4)", 172, 63],
                       ["surface", "surface, not interface (> 0.4)", 292, 52],
                       ["interface", "interface (exposed, touches another chain)", 56, 118],
                       ["ligand", "ligand-binding residue", 14, 85]],
    # notebook 05, 14f / 14g: Design A (the model predicts the tolerance score), within-gene Spearman, validation genes
    "design_a": {"esm_zero_shot": 0.020, "from_scratch": [0.03, 0.05], "esm_heads": [0.09, 0.125], "best": 0.125,
                 "trees_all_inputs": 0.126, "ceiling": 0.39,
                 "pathogenic_only": {"ceiling_before": 0.37, "ceiling_after": 0.61, "probe_before": 0.10,
                                     "probe_after": 0.32, "auroc_before": 0.79, "auroc_after": 0.84}},
}

# ------------------------------------------------------------------ conservation across species (notebook 06 §3)
if not args.skip_conservation:
    import evolution_conservation as ec
    from Bio import Phylo
    from io import StringIO
    from species_tree import build_tree
    cache = os.path.join(RES, "poster", "conservation_cache")
    mane = pd.read_csv(args.mane)
    mane = mane[mane.MANE_status == "MANE Select"]
    ensg = dict(zip(mane.symbol, mane.Ensembl_Gene))
    sets = {"Postsynaptic proteins": [ensg[g] for g in seq_df.gene if g in ensg],
            "Olfactory receptors": [ensg[s] for s in mane.symbol if re.fullmatch(r"OR\d+[A-Z]+\d+[A-Z]?", str(s))],
            "Ribosomal proteins": [ensg[s] for s in mane.symbol
                                   if re.fullmatch(r"RP[LS]\d+[AXYL]?\d*|RPLP[0-2]|RPSA", str(s))]}
    value, coverage = {}, {}
    for name, genes in sets.items():
        per_gene = ec.protein_similarity(genes, os.path.join(cache, "similarity"))
        value[name] = ec.summarise(per_gene.where(per_gene > 0).median().to_frame().T)
        coverage[name] = ec.summarise(((per_gene > 0).mean() * 100).to_frame().T)
    aligned = sorted({sp for f in os.listdir(os.path.join(cache, "dna")) if f.endswith(".json")
                      for sp in json.load(open(os.path.join(cache, "dna", f)))})
    for name, kind in [("Regulatory elements", "regulatory"), ("Coding regions", "cds")]:
        regions = [tuple(x) for x in json.load(open(os.path.join(cache, f"{kind}_regions.json")))]
        value[name] = ec.summarise(ec.dna_identity(regions, os.path.join(cache, "dna"), aligned))
    tree = build_tree(keep={"Human", *ec.SPECIES}); tree.ladderize(reverse=True)
    tips = [c.name for c in tree.get_terminals()]
    tab = pd.DataFrame(value).reindex(tips); tab.loc["Human"] = 100.0
    cov = pd.DataFrame(coverage).reindex(tips); cov.loc["Human"] = 100.0
    nwk = StringIO(); Phylo.write(tree, nwk, "newick")
    clade_names = {}
    for c in tree.find_clades():
        if getattr(c, "clade_name", None):
            clade_names[c.clade_name] = [t.name for t in c.get_terminals()]
    summary["conservation"] = {
        "species": tips, "columns": list(tab.columns),
        "values": [[None if pd.isna(x) else round(float(x), 2) for x in tab.loc[s]] for s in tips],
        "coverage": [[None if pd.isna(cov.loc[s].get(c, np.nan)) else round(float(cov.loc[s][c]), 0)
                      for c in tab.columns] for s in tips],
        "newick": nwk.getvalue().strip(), "clades": clade_names,
        "n_genes": {k: len(v) for k, v in sets.items()}, "n_regions": 80, "ensembl_release": 116}
    print("conservation table:\n", tab.round(0).to_string())

write_js(os.path.join(OUT, "summary.js"), "TOL.summary", summary)
print("splits:", summary["split"])
print("variants:", summary["variants"])
