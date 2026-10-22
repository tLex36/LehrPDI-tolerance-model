"""Fetch for the revised conservation figure: protein SIMILARITY (REST, perc_pos) for the three gene sets,
and DNA identity of 80 random coding exons (CDS). Cached; safe to rerun - only missing items are fetched."""
import os, re, sys, time
import pandas as pd
sys.path.insert(0, "/Users/awlehr/LehrPDI/figures")
import evolution_conservation as ec

CACHE = "/Users/awlehr/LehrPDI/figures/results/tolerance_model/poster/conservation_cache"
mane = pd.read_csv("/Users/awlehr/LehrPDI/mane_transcripts.csv")
mane = mane[mane.MANE_status == "MANE Select"]
ensg = dict(zip(mane.symbol, mane.Ensembl_Gene))
psd = pd.read_parquet("/Users/awlehr/LehrPDI/figures/results/tolerance_model/tolerance_v4.parquet", columns=["gene"]).gene.unique()
sets = {"postsynaptic": [ensg[g] for g in psd if g in ensg],
        "olfactory": [ensg[s] for s in mane.symbol if re.fullmatch(r"OR\d+[A-Z]+\d+[A-Z]?", str(s))],
        "ribosomal": [ensg[s] for s in mane.symbol if re.fullmatch(r"RP[LS]\d+[AXYL]?\d*|RPLP[0-2]|RPSA", str(s))]}
t0 = time.time()
species = ec.species_with_alignment()
regions = ec.sample_cds(mane.Ensembl_nuc.tolist(), 80, os.path.join(CACHE, "cds_regions.json"), seed=0, threads=6)
print(f"cds: {len(regions)} regions sampled ({time.time() - t0:.0f} s)", flush=True)
df = ec.dna_identity(regions, os.path.join(CACHE, "dna"), species, threads=6)
print(f"cds: done, {len(df)} regions ({time.time() - t0:.0f} s)", flush=True)
for name, genes in sets.items():
    df = ec.protein_similarity(genes, os.path.join(CACHE, "similarity"), threads=6)
    print(f"{name}: {len(df)} of {len(genes)} genes ({time.time() - t0:.0f} s)", flush=True)
print("ALL DONE", flush=True)
