"""One-time fetch of the evolution-conservation data (cached; safe to rerun - only missing items are fetched)."""
import os, re, sys, time
import pandas as pd
sys.path.insert(0, "/Users/awlehr/LehrPDI/figures")
import evolution_conservation as ec

CACHE = "/Users/awlehr/LehrPDI/figures/results/tolerance_model/poster/conservation_cache"
mane = pd.read_csv("/Users/awlehr/LehrPDI/mane_transcripts.csv")
ensg = dict(zip(mane.symbol, mane.Ensembl_Gene))
psd = pd.read_parquet("/Users/awlehr/LehrPDI/figures/results/tolerance_model/tolerance_v4.parquet", columns=["gene"]).gene.unique()
sets = {"postsynaptic": [ensg[g] for g in psd if g in ensg],
        "olfactory": [ensg[s] for s in mane.symbol if re.fullmatch(r"OR\d+[A-Z]+\d+[A-Z]?", str(s))],
        "ribosomal": [ensg[s] for s in mane.symbol if re.fullmatch(r"RP[LS]\d+[AXYL]?\d*|RPLP[0-2]|RPSA", str(s))]}
t0 = time.time()
for name, genes in sets.items():
    df = ec.protein_identity(genes, os.path.join(CACHE, "orthologues"))
    print(f"{name}: {len(df)} of {len(genes)} genes ({time.time() - t0:.0f} s)", flush=True)
species = ec.species_with_alignment()
print("species with a LASTZ_NET alignment to human:", species, flush=True)
for kind in ["regulatory", "genome"]:
    regions = ec.sample_regions(kind, 80, os.path.join(CACHE, f"{kind}_regions.json"), seed=0)
    print(f"{kind}: {len(regions)} regions sampled ({time.time() - t0:.0f} s)", flush=True)
    df = ec.dna_identity(regions, os.path.join(CACHE, "dna"), species, threads=16)
    print(f"{kind}: done, {len(df)} regions ({time.time() - t0:.0f} s)", flush=True)
print("ALL DONE", flush=True)
