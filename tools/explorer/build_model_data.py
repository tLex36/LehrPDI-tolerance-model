"""Explorer data that need notebook 05's working state: the held-out experiment (§6g), the final model's weights and
how each feature moves each labelled variant's prediction. Rebuilds that state from the notebook's own cells, as
poster_figs.py does (about 3 minutes and 4-5 GB of memory; no kernel needed).

    python tools/explorer/build_model_data.py --notebooks /path/to/LehrPDI/figures/notebooks

Writes docs/data/model.js (held-out histograms, weights) and docs/data/model/<GENE>.js (one file per training or
validation gene with labelled variants).
"""
import argparse, json, os, sys, time
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
ap = argparse.ArgumentParser()
ap.add_argument("--notebooks", default=os.path.join(REPO, "figures", "notebooks"),
                help="folder holding 05_tolerance_score_modelling.ipynb, with its data in place")
ap.add_argument("--out", default=os.path.join(REPO, "docs", "data"))
args = ap.parse_args()
OUT_DIR = os.path.abspath(args.out)
os.makedirs(os.path.join(OUT_DIR, "model"), exist_ok=True)
sys.path.insert(0, HERE)
from jsdata import write_js                                        # noqa: E402  (same JSON-in-a-call format)

# ------------------------------------------------------------------ notebook state (as in poster_figs.py)
os.chdir(args.notebooks)
import matplotlib
matplotlib.use("Agg")
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
AA20 = "ACDEFGHIKLMNPQRSTVWY"
ANN_COLS = [c for c in residue_df.columns if c.startswith("ann_") and c != "ann_disorder_curated"] + ["n_annotations"]
exec(src("# 14h setup:"), globals())                                # lab: 15,339 variants x 42 features
print(f"state rebuilt in {time.time() - t0:.0f} s: {len(lab):,} labelled variants", flush=True)

# ------------------------------------------------------------------ 1. the held-out experiment (§6g)
# Tolerance rebuilt from 90% of the variants, scored at the 10% held out. Groups as in the poster (M2): each variant
# once; a variant in both databases counts as ClinVar.
v = val_scored.sort_values("database").drop_duplicates("id").copy()
v["group"] = np.where(v.database == "GnomAD", "gnomad", v.clinsig_group.astype(str).replace({"vus": "vus"}))
v.loc[(v.database != "GnomAD") & ~v.group.isin(["benign", "vus", "pathogenic"]), "group"] = "other"
# the poster panel's subset: variants in one database only, at residues where another substitution was observed
poster = val_scored[~val_scored.in_both & val_scored.esm_obs_min_llr.notna()].copy()
poster["group"] = np.where(poster.database == "GnomAD", "gnomad", poster.clinsig_group.astype(str))
poster.loc[(poster.database != "GnomAD") & ~poster.group.isin(["benign", "vus", "pathogenic"]), "group"] = "other"
SCORES = {"tolerance": ("tolerance_v4_train", -1.0, 1.0, 0.01), "esm": ("esm_llr", -26.0, 10.0, 0.1)}
GROUPS = ["gnomad", "benign", "vus", "pathogenic", "other"]


def histograms(df):
    out = {"n": {}, "scores": {}}
    for key, (colname, lo, hi, step) in SCORES.items():
        edges = np.round(np.arange(lo, hi + step / 2, step), 4)
        block = {"lo": lo, "hi": hi, "step": step, "groups": {}}
        for g in GROUPS:
            x = df.loc[df.group == g, colname].dropna().to_numpy(float)
            if not len(x):
                continue
            counts, _ = np.histogram(np.clip(x, lo, hi - 1e-9), bins=edges)
            block["groups"][g] = {"n": int(len(x)), "mean": round(float(x.mean()), 4),
                                  "median": round(float(np.median(x)), 4), "counts": counts.astype(int).tolist()}
        # exact AUROCs (higher score = more benign side), for checking the page's binned estimate
        exact = {}
        for i, a in enumerate(GROUPS):
            for b in GROUPS[i + 1:]:
                xa = df.loc[df.group == a, colname].dropna().to_numpy(float)
                xb = df.loc[df.group == b, colname].dropna().to_numpy(float)
                if len(xa) and len(xb):
                    exact[f"{a}|{b}"] = round(float(roc_auc_score(np.r_[np.ones(len(xa)), np.zeros(len(xb))],
                                                                  np.r_[xa, xb])), 4)
        block["exact_auroc"] = exact
        out["scores"][key] = block
    out["n"] = {g: int((df.group == g).sum()) for g in GROUPS}
    return out


heldout = {"all": histograms(v), "poster": histograms(poster),
           "n_heldout_variants": int(v.id.nunique()), "n_genes": int(v.gene.nunique())}
for k, blk in heldout["all"]["scores"].items():
    print(k, {g: (d["n"], d["mean"]) for g, d in blk["groups"].items()})

# ------------------------------------------------------------------ 2. the final model and its weights (§16c)
final = make_pipeline(SimpleImputer(strategy="median", add_indicator=True), StandardScaler(),
                      LogisticRegression(max_iter=5000))
X = lab[ALL_VAR_FEATURES].to_numpy(float)
final.fit(X, lab.benign.to_numpy())
imputer, scaler, lr = final[0], final[1], final[2]
names = list(final[:-1].get_feature_names_out())                    # x0..x41, then missingindicator_x..
coef = lr.coef_[0]
GROUP_OF = {f: g for g, fs in VARIANT_GROUPS.items() for f in fs}
PLAIN = {
    "esm_llr": "ESM-1v LLR of this substitution", "esm_wt_lp": "ESM log-probability of the original amino acid",
    "esm_entropy": "ESM entropy at the residue (spread of its guesses)", "esm_mean_llr": "ESM mean LLR of all 19 swaps",
    "tolerance": "tolerance score (own residue hidden)", "gnomad_density": "smoothed gnomAD density",
    "clinvar_density": "smoothed ClinVar density", "gnomad_within_5": "gnomAD variants within ±5 residues",
    "gnomad_within_15": "gnomAD variants within ±15 residues", "clinvar_within_5": "ClinVar variants within ±5 residues",
    "clinvar_within_15": "ClinVar variants within ±15 residues",
    "gnomad_per_residue": "gnomAD variants per residue (whole gene)",
    "clinvar_per_residue": "ClinVar variants per residue (whole gene)",
    "blosum62": "BLOSUM62 dissimilarity of the swap", "grantham": "Grantham distance of the swap",
    "hydrophobicity_change": "change in hydrophobicity (Kyte-Doolittle)", "rsa": "solvent accessibility (RSA)",
    "plddt": "AlphaFold confidence (pLDDT)", "buried_sasa": "surface buried by partners (few genes)",
    "iface_any": "touches another chain (interface)", "iface_homo": "touches a copy of itself",
    "iface_hetero": "touches a different protein", "disorder": "disorder (AIUPred)",
    "binding": "disordered binding region (AIUPred)", "ann_disorder_curated": "curated disorder (MobiDB)",
    "msa_conservation": "conservation across related sequences", "n_annotations": "number of annotations at the residue",
}
features = []
for i, n in enumerate(names):
    if n.startswith("x"):
        f = ALL_VAR_FEATURES[int(n[1:])]
        label = PLAIN.get(f, f.replace("ann_", "annotation: ").replace("_", " "))
        features.append({"key": f, "label": label, "group": GROUP_OF[f], "weight": round(float(coef[i]), 4),
                         "missing_flag": False})
    else:                                                           # missingindicator_x7 -> "rsa is missing"
        f = ALL_VAR_FEATURES[int(n.split("_x")[-1])]
        features.append({"key": f"{f}_missing", "label": f"{PLAIN.get(f, f)} missing",
                         "group": "missing-value flags", "weight": round(float(coef[i]), 4), "missing_flag": True})
w = pd.Series({f["key"]: f["weight"] for f in features if not f["missing_flag"]})
top = w.reindex(w.abs().sort_values(ascending=False).index[:15]).round(3)
ref = json.load(open(os.path.join(OUT, "poster", "poster_numbers.json")))["weights_top15"]
if any(not abs(top.get(k, np.nan) - val) <= 2e-3 for k, val in ref.items()):   # poster rounded to 3 dp
    sys.exit(f"weights differ from the poster figure:\n{top}\n{ref}")
print("top-15 weights match the poster (Figure 26)")

# ------------------------------------------------------------------ 3. each labelled variant, feature by feature
Xs = scaler.transform(imputer.transform(X))                          # what the model actually multiplies
contrib = Xs * coef                                                 # one column per model input
z = contrib.sum(1) + lr.intercept_[0]
p = 1 / (1 + np.exp(-z))
assert np.allclose(p, final.predict_proba(X)[:, 1])
print(f"in-sample AUROC {roc_auc_score(lab.benign, p):.3f} (cross-validated: 0.961)")
TOP_K = 12
split_of = json.load(open(os.path.join(OUT, "split_manifest.json")))["split"]
written = 0
for gene, idx in lab.groupby("gene").indices.items():
    rows = []
    for i in idx[np.argsort(lab.position.to_numpy()[idx])]:
        order = np.argsort(-np.abs(contrib[i]))[:TOP_K]
        raw = [(int(not np.isfinite(X[i, ALL_VAR_FEATURES.index(features[j]["key"][:-len("_missing")])]))
                if features[j]["missing_flag"] else
                None if not np.isfinite(X[i, j]) else float(np.format_float_positional(
                    X[i, j], precision=3, unique=False, fractional=False, trim="-"))) for j in order]
        pc = lab.protein_consequence.iloc[i]
        rows.append({"v": pc[2:] if pc.startswith("p.") else pc, "benign": int(lab.benign.iloc[i]),
                     "stars": int(lab.stars.iloc[i]), "p": round(float(p[i]), 4), "z": round(float(z[i]), 3),
                     "top": [[int(j), round(float(contrib[i, j]), 3), r] for j, r in zip(order, raw)],
                     "rest": round(float(contrib[i].sum() - contrib[i, order].sum()), 3)})
    write_js(os.path.join(OUT_DIR, "model", f"{gene}.js"), "TOL.model",
             {"gene": gene, "split": split_of.get(gene), "variants": rows})
    written += 1

write_js(os.path.join(OUT_DIR, "model.js"), "TOL.modelSummary", {
    "heldout": heldout, "features": features, "intercept": round(float(lr.intercept_[0]), 4),
    "n_training_variants": int(len(lab)), "in_sample_auroc": round(float(roc_auc_score(lab.benign, p)), 4),
    "top_k": TOP_K, "groups": list(VARIANT_GROUPS) + ["missing-value flags"]})
print(f"wrote model.js and {written} per-gene model files in {time.time() - t0:.0f} s")
