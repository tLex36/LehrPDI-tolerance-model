"""Loop annotate_pdb over NLGN crystal structures and consolidate onto human numbering.

    results = run_structure_pipeline()             # {gene: {...}}
    results["NLGN1"]["ranges"]["core"]             # [(s, e), ...] human UniProt numbering
    results["NLGN1"]["ranges"]["surface"]
    results["NLGN1"]["ranges"]["interface"]        # surface residues in any interface
    results["NLGN1"]["ranges"]["iface_NRXN1B"]     # per partner
    results["NLGN1"]["table"]                      # per-residue DataFrame
"""
from netutil import pkg_path
import json
import numpy as np, pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
from structure import annotate_pdb, fetch_uniprot_seq

# gene -> [(pdb, chain, partners)]; partner chains verified with structure.find_contacts
STRUCTURES = {
    "NLGN1":  [("3B3Q", "A", {"dimer": ["B"], "NRXN1B": ["E"]}),      # mouse NLGN1 + human NRXN1β
               ("3BIW", "A", {"dimer": ["D"], "NRXN1B": ["E"]}),      # rat NLGN1 + rat NRXN1β
               ("5OJ6", "A", {"MDGA1": ["B"]})],                       # human NLGN1 + chicken MDGA1
    "NLGN2":  [("5V5V", "A", {"dimer": ["B"], "MDGA1": ["H"]})],      # rat NLGN2 + human MDGA1
    "NLGN3":  [("8GS3", "A", {"dimer": ["B"]}),                       # human NLGN3 cryo-EM (ECD resolved: 41-633)
               ("7CEE", "A", {"dimer": ["B"]})],                       # mouse NLGN3 crystal
    "NLGN4X": [("2WQZ", "A", {"dimer": ["B"], "NRXN1B": ["D"]}),      # human NLGN4X + rat NRXN1β
               ("3BE8", "A", {"dimer": ["B"]})],                       # human NLGN4X apo
}

def map_all_positions(seq_from, seq_to) -> dict:
    """{1-based resnum in seq_from: 1-based resnum in seq_to} via global alignment (aligned residues only)."""
    a = PairwiseAligner(); a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score, a.extend_gap_score = -10, -0.5; a.mode = "global"
    aln = a.align(seq_from, seq_to)[0]
    m = {}
    for (f0, f1), (t0, t1) in zip(*aln.aligned):
        for k in range(f1 - f0):
            m[f0 + k + 1] = t0 + k + 1
    return m

def to_ranges(nums):
    nums = sorted(set(int(n) for n in nums)); out = []
    for n in nums:
        if out and n == out[-1][1] + 1: out[-1] = (out[-1][0], n)
        else: out.append((n, n))
    return out

def consolidate(annotations, human_seq, core=0.20, surface=0.40):
    """Merge per-structure annotations onto human numbering.
    RSA: median across structures -> class.  Interfaces: union across structures, per label."""
    per_struct = []
    for a in annotations:
        df = a["residues"].copy()
        if a["species"] != "HUMAN":
            m = map_all_positions(a["sequence"], human_seq)
            df["human_resnum"] = df.unp_resnum.map(m)
            df = df.dropna(subset=["human_resnum"])
        else:
            df["human_resnum"] = df.unp_resnum
        df["human_resnum"] = df.human_resnum.astype(int)
        df["structure"] = f'{a["pdb"]}:{a["chain"]}'
        per_struct.append(df)
    long = pd.concat(per_struct, ignore_index=True)
    iface_cols = sorted({c for c in long.columns if c.startswith("iface_")})

    t = long.groupby("human_resnum").agg(rsa=("rsa", "median"), n_structures=("structure", "nunique"))
    t["aa"] = [human_seq[i - 1] for i in t.index]
    t["class"] = np.where(t.rsa < core, "core", np.where(t.rsa > surface, "surface", "boundary"))
    for c in iface_cols:
        t[c] = long.groupby("human_resnum")[c].any().reindex(t.index).fillna(False).astype(bool)
    t["interface"] = t[iface_cols].any(axis=1) if iface_cols else False
    t = t[["aa", "rsa", "class", "n_structures", *iface_cols, "interface"]]

    ranges = {c: to_ranges(t.index[t["class"] == c]) for c in ("core", "boundary", "surface")}
    ranges["interface"] = to_ranges(t.index[t.interface & (t["class"] != "core")])   # surface/boundary residues that bind
    for c in iface_cols:
        ranges[c] = to_ranges(t.index[t[c]])
    return {"table": t, "ranges": ranges, "per_structure": long}

def run_structure_pipeline(structures=STRUCTURES, accessions=None, verbose=True):
    accessions = accessions or pkg_path("nlgn_accessions.json")
    acc = json.load(open(accessions)) if isinstance(accessions, str) else accessions
    results = {}
    for gene, entries in structures.items():
        human_seq = fetch_uniprot_seq(acc[gene]["Human"])
        annots = []
        for pdb, chain, partners in entries:
            if verbose: print(f"{gene}: {pdb} chain {chain} ...", end=" ", flush=True)
            a = annotate_pdb(pdb, chain, partners=partners)
            if verbose: print(f'{a["entry"]}, {len(a["residues"])} residues')
            annots.append(a)
        if not annots:
            if verbose: print(f"{gene}: no structures configured - skipped")
            continue
        r = consolidate(annots, human_seq)
        r["structures"] = [f'{a["pdb"]}:{a["chain"]} ({a["entry"]})' for a in annots]
        r["human_sequence"] = human_seq
        results[gene] = r
    return results

if __name__ == "__main__":
    res = run_structure_pipeline()
    for g, r in res.items():
        t = r["table"]
        print(f"\n== {g}  structures: {r['structures']}")
        print("   residues covered:", t.index.min(), "-", t.index.max(), f"({len(t)})")
        print("   class:", t["class"].value_counts().to_dict())
        for k, v in r["ranges"].items():
            print(f"   {k:14} {sum(e - s + 1 for s, e in v):3} residues in {len(v):3} ranges")

# ---------------------------------------------------------------- RSA lookups
def _expand(spec):
    """471 | (500, 505) | [471, (500, 505), ...] -> flat list of residue numbers.
    A bare 2-tuple/2-list of ints is a range, not two separate residues."""
    if isinstance(spec, int):
        spec = [spec]
    elif isinstance(spec, (tuple, list)) and len(spec) == 2 and all(isinstance(x, int) for x in spec):
        spec = [tuple(spec)]                            # (start, end) -> one range item
    out = []
    for item in spec:
        if isinstance(item, (tuple, list)) and len(item) == 2:
            out.extend(range(int(item[0]), int(item[1]) + 1))
        else:
            out.append(int(item))
    return out

def rsa_table(results, selection, columns=None, missing="warn"):
    """RSA and structural class for chosen residues.

    selection accepts:
      {"NLGN3": [471, 604]}                       -> those residues of that gene
      {"NLGN3": {"variants": [471], "loop": (770, 790)}}  -> named sets, kept as a `set` column
      {"NLGN1": [(450, 470)], "NLGN3": [471]}     -> several genes at once
    Ranges are (start, end) inclusive, in human UniProt numbering.

    Returns a tidy DataFrame: gene, set, resnum, aa, rsa, class, n_structures, interface flags.
    Residues with no structural coverage are reported with NaN (missing='warn' prints them).
    """
    rows, absent = [], []
    for gene, spec in selection.items():
        if gene not in results:
            absent.append((gene, "no structure annotation")); continue
        t = results[gene]["table"]
        named = spec if isinstance(spec, dict) else {"": spec}
        for label, sub in named.items():
            for r in _expand(sub):
                if r in t.index:
                    row = t.loc[r].to_dict()
                    row.update({"gene": gene, "set": label, "resnum": r})
                else:
                    row = {"gene": gene, "set": label, "resnum": r, "aa": None,
                           "rsa": float("nan"), "class": None, "n_structures": 0}
                    absent.append((gene, r))
                rows.append(row)
    df = pd.DataFrame(rows)
    front = [c for c in ("gene", "set", "resnum", "aa", "rsa", "class", "n_structures") if c in df]
    rest = [c for c in df.columns if c not in front]
    df = df[front + (columns or rest)]
    if missing == "warn" and absent:
        print(f"{len(absent)} residue(s) outside the resolved structure: {absent[:10]}"
              + (" ..." if len(absent) > 10 else ""))
    return df.reset_index(drop=True)

def rsa_summary(results, selection, by=("gene", "set")):
    """Mean/median RSA and class composition per named set."""
    t = rsa_table(results, selection, missing="ignore")
    g = t.groupby(list(by), dropna=False)
    out = g.agg(n=("resnum", "size"), mean_rsa=("rsa", "mean"), median_rsa=("rsa", "median"),
                n_core=("class", lambda s: (s == "core").sum()),
                n_boundary=("class", lambda s: (s == "boundary").sum()),
                n_surface=("class", lambda s: (s == "surface").sum()))
    return out.round(3).reset_index()
