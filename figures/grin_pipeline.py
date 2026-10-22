"""NMDA receptor (GRIN) structure + variant pipeline.

    from grin_pipeline import run_grin_structures, grin_region_density
    gres = run_grin_structures()                     # core/surface/subunit interface
    tab  = grin_region_density(gres, mis)            # densities incl. TM and CTD

Structures are full-length human cryo-EM receptors, so they resolve ATD+LBD+TM but not the
disordered CTD. RSA is only meaningful outside the membrane, so core/surface are computed for
extracellular residues only; TM and CTD come from UniProt topology.
"""
from netutil import pkg_path
import json
import numpy as np, pandas as pd
from structure import annotate_pdb
from structure_pipeline import consolidate
from conservation import fetch_uniprot, uniprot_topology
from variants import residue_table

# gene -> [(pdb, chain, {label: partner chains})]; chains verified with structure.find_contacts
GRIN_STRUCTURES = {
    "GRIN1":  [("7EU7", "A", {"subunit": ["B", "D"]})],   # human GluN1/GluN2A
    "GRIN2A": [("7EU7", "B", {"subunit": ["A", "C"]})],
    "GRIN2B": [("9D37", "B", {"subunit": ["A", "C"]})],     # human GluN1/GluN2B/GluN2D
    "GRIN2D": [("8E96", "B", {"subunit": ["A", "C"]})],     # human GluN1/GluN2D
}
# The parquet's GRIN positions match UniProt canonical numbering exactly (verified: 100% of
# reference residues agree for all seven genes), so no Ensembl/MANE round trip is needed here.

def run_grin_structures(structures=GRIN_STRUCTURES, accessions=None, verbose=True):
    accessions = accessions or pkg_path("grin_accessions.json")
    acc = json.load(open(accessions)) if isinstance(accessions, str) else accessions
    results = {}
    for gene, entries in structures.items():
        human_seq = fetch_uniprot({"x": acc["_human"][gene]})["x"]
        annots = []
        for pdb, chain, partners in entries:
            if verbose: print(f"{gene}: {pdb} chain {chain} ...", end=" ", flush=True)
            a = annotate_pdb(pdb, chain, partners=partners)
            if verbose: print(f'{a["entry"]}, {len(a["residues"])} residues')
            annots.append(a)
        r = consolidate(annots, human_seq)
        r["structures"] = [f'{a["pdb"]}:{a["chain"]} ({a["entry"]})' for a in annots]
        r["human_sequence"] = human_seq
        results[gene] = r
    return results

def grin_regions(results, gene, mane, topo, min_rsa_core=0.20, max_rsa_surface=0.40):
    """{region: set of MANE residue numbers}: core/surface (extracellular only, from structure),
    subunit interfaces, TM and CTD (from UniProt topology transferred to MANE numbering)."""
    seq = mane[gene]
    unp = results[gene]["human_sequence"] if gene in results else None
    # topology is in UniProt numbering -> map to MANE
    from structure_pipeline import map_all_positions
    if unp is None or unp == seq:
        to_mane = {i: i for i in range(1, len(seq) + 1)}
    else:
        to_mane = map_all_positions(unp, seq)
    def rng(pairs):
        return {to_mane[i] for s, e in pairs for i in range(s, e + 1) if i in to_mane}
    regions = {"TM": rng(topo["transmembrane"]),
               "CTD": rng([topo["ctd"]]) if topo["ctd"] else set()}
    extracellular = rng(topo["extracellular"])
    if gene in results:
        t = residue_table(results, gene, seq)              # index = MANE resnum
        ecd = t.index.isin(extracellular)
        regions["core"] = set(t.index[ecd & (t["class"] == "core")])
        regions["surface"] = set(t.index[ecd & (t["class"] == "surface")])
        for c in [c for c in t.columns if c.startswith("iface_")]:
            if t[c].any():
                regions[c.replace("iface_", "")] = set(t.index[t[c]])
    return regions

GRIN_REGION_ORDER = ["core", "surface", "subunit", "TM", "CTD"]

def grin_region_density(results, missense, genes=None, accessions="grin_accessions.json",
                        source="gnomAD", long=False):
    """Missense density per region per GRIN gene (variants per residue)."""
    acc = json.load(open(accessions)) if isinstance(accessions, str) else accessions
    genes = genes or list(results)
    mane = fetch_uniprot({g: acc["_human"][g] for g in genes})
    rows = []
    for g in genes:
        topo = uniprot_topology(acc["_human"][g])
        regions = grin_regions(results, g, mane, topo)
        v = missense[missense.gene == g]
        if source: v = v[v.source == source]
        for name, resi in regions.items():
            if not resi: continue
            sub = v[v.position.isin(resi)]
            rows.append({"gene": g, "region": name, "n_residues": len(resi), "n_variants": len(sub),
                         "density": len(sub) / len(resi),
                         "mean_blosum": sub.blosum62.mean() if len(sub) else np.nan,
                         "severity_density": sub.severity.clip(lower=0).sum() / len(resi)})
    df = pd.DataFrame(rows)
    if long: return df
    cols = sorted(df.region.unique(),
                  key=lambda r: GRIN_REGION_ORDER.index(r) if r in GRIN_REGION_ORDER else 99)
    return df.pivot(index="gene", columns="region", values="density").reindex(
        index=list(genes), columns=cols).round(3)
