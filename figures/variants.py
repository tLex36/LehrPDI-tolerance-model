"""Missense-variant density by structural residue class, BLOSUM62-weighted.

    from variants import load_missense, variant_class_table
    mis = load_missense("../Model Training/missense_df.parquet")
    tab = variant_class_table(results, mis)        # results = structure_pipeline output

Positions in the parquet follow MANE Select protein numbering; the structure annotation is in
UniProt numbering. mane_sequences() fetches the MANE proteins and the tables are joined by
aligning the two sequences, so no numbering assumptions are made.
"""
from netutil import pkg_path
import json, os, time, urllib.request
from netutil import urlopen_retry
import numpy as np, pandas as pd
from Bio.Align import substitution_matrices
from structure_pipeline import map_all_positions, to_ranges
from structure import fetch_uniprot_seq

BLOSUM62 = substitution_matrices.load("BLOSUM62")

# MANE Select Ensembl protein ids (from mane_transcripts.csv)
MANE_PROT = {"NLGN1": "ENSP00000511841", "NLGN2": "ENSP00000305288", "NLGN3": "ENSP00000351591",
             "NLGN4X": "ENSP00000370485", "NLGN4Y": "ENSP00000510011"}

def mane_sequences(genes=None, ids=MANE_PROT, cache=None):
    """{gene: MANE Select protein sequence}, cached locally."""
    cache = cache or pkg_path("mane_proteins.json")
    try:    seqs = json.load(open(cache))
    except Exception: seqs = {}
    for g in (genes or ids):
        if g in seqs: continue
        url = f"https://rest.ensembl.org/sequence/id/{ids[g]}?type=protein;content-type=text/x-fasta"
        fa = urlopen_retry(url).decode()
        seqs[g] = "".join(fa.split("\n")[1:]).strip().rstrip("*")
    json.dump(seqs, open(cache, "w"))
    return {g: seqs[g] for g in (genes or ids) if g in seqs}

def load_missense(path, genes=("NLGN1", "NLGN2", "NLGN3", "NLGN4X", "NLGN4Y")):
    """Read the parquet, keep NLGN missense rows, add BLOSUM62 score and severity weight."""
    df = pd.read_parquet(path)
    df = df[df.gene.isin(genes)].copy()
    aas = set(BLOSUM62.alphabet)
    df = df[df.reference.isin(aas) & df.alternate.isin(aas)]
    df["blosum62"] = [BLOSUM62[r][a] for r, a in zip(df.reference, df.alternate)]
    df["severity"] = -df.blosum62                      # radical substitution -> high weight
    df["source"] = np.where(df.database.str.lower().str.startswith("clin"), "ClinVar", "gnomAD")
    return df

def residue_table(results, gene, mane_seq, uniprot_acc=None):
    """Per-residue structural annotation in MANE numbering for one gene."""
    t = results[gene]["table"].copy()                  # index = human UniProt resnum
    unp = results[gene].get("human_sequence") or fetch_uniprot_seq(uniprot_acc)
    if unp == mane_seq:
        t["mane_resnum"] = t.index
    else:
        m = map_all_positions(unp, mane_seq)
        t["mane_resnum"] = [m.get(i) for i in t.index]
        t = t.dropna(subset=["mane_resnum"])
        t["mane_resnum"] = t.mane_resnum.astype(int)
    bad = [i for i, r in zip(t.mane_resnum, t.aa) if mane_seq[i - 1] != r]
    if bad:
        print(f"  {gene}: {len(bad)} residues disagree with the MANE sequence after mapping")
    return t.set_index("mane_resnum")

def variant_class_table(results, missense, genes=None, mane=None, classes=("core", "boundary", "surface"),
                        split_interface=True):
    """Counts and BLOSUM62-weighted densities per (gene, source, residue class).

    Columns: n_residues, n_variants, density (variants/residue), mean_blosum (mean BLOSUM62 of
    observed substitutions; lower = more radical changes seen), severity_density (sum of
    severity per residue), and for gnomAD alleles / allele_density (total_count-weighted).
    """
    genes = genes or [g for g in results if g in (mane or MANE_PROT)]
    mane = mane or mane_sequences(genes)
    rows = []
    for g in genes:
        t = residue_table(results, g, mane[g])
        if split_interface:                            # non-core binding residues as their own class
            t["group"] = np.where(t.interface & (t["class"] != "core"), "interface", t["class"])
        else:
            t["group"] = t["class"]
        v = missense[missense.gene == g].merge(t[["group"]], left_on="position", right_index=True)
        for src in ("gnomAD", "ClinVar"):
            vs = v[v.source == src]
            for grp, n_res in t.group.value_counts().items():
                sub = vs[vs.group == grp]
                row = {"gene": g, "source": src, "group": grp, "n_residues": int(n_res),
                       "n_variants": len(sub),
                       "density": len(sub) / n_res,
                       "mean_blosum": sub.blosum62.mean() if len(sub) else np.nan,
                       "severity_density": sub.severity.clip(lower=0).sum() / n_res}
                if src == "gnomAD":
                    row["alleles"] = sub.total_count.fillna(0).sum()
                    row["allele_density"] = row["alleles"] / n_res
                rows.append(row)
    df = pd.DataFrame(rows)
    order = ["core", "boundary", "surface", "interface"]
    df["group"] = pd.Categorical(df.group, [c for c in order if c in set(df.group)], ordered=True)
    return df.sort_values(["gene", "source", "group"]).reset_index(drop=True)

def pooled(tab):
    """Sum across genes -> one row per (source, group), plus the ClinVar:gnomAD density ratio."""
    g = (tab.groupby(["source", "group"], observed=True)
            .agg(n_residues=("n_residues", "sum"), n_variants=("n_variants", "sum"),
                 severity=("severity_density", lambda s: np.nan))
            .reset_index())
    # recompute weighted quantities from the raw sums
    raw = tab.set_index(["source", "group"])
    g["density"] = g.n_variants / g.n_residues
    piv = g.pivot(index="group", columns="source", values="density")
    piv["clinvar_per_gnomad"] = piv["ClinVar"] / piv["gnomAD"]
    return g, piv.round(3)

def plot_variant_classes(tab, metric="density", ax=None, cmap=("#4C78A8", "#E45756"), title=None):
    """Grouped bars: one group per residue class, one bar per source, faceted by gene."""
    import matplotlib.pyplot as plt
    genes = list(dict.fromkeys(tab.gene))
    groups = [g for g in ("core", "boundary", "surface", "interface") if g in set(tab.group)]
    if ax is None:
        fig, ax = plt.subplots(figsize=(1.6 * len(genes) * len(groups) / 2 + 2, 3.6))
    w, x0 = 0.38, 0
    ticks, labels = [], []
    for gene in genes:
        for grp in groups:
            for k, src in enumerate(("gnomAD", "ClinVar")):
                r = tab[(tab.gene == gene) & (tab.group == grp) & (tab.source == src)]
                if len(r):
                    ax.bar(x0 + k * w, float(r[metric].iloc[0]), w, color=cmap[k],
                           label=src if (gene == genes[0] and grp == groups[0]) else None)
            ticks.append(x0 + w / 2); labels.append(grp)
            x0 += 1.1
        x0 += 0.5
    ax.set_xticks(ticks); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel({"density": "variants per residue", "severity_density": "BLOSUM62 severity per residue",
                   "mean_blosum": "mean BLOSUM62"}.get(metric, metric))
    n = len(groups)
    for i, gene in enumerate(genes):                   # gene label under each block
        ax.text(ticks[i * n] + (ticks[i * n + n - 1] - ticks[i * n]) / 2, -0.22, gene,
                transform=ax.get_xaxis_transform(), ha="center", va="top", fontweight="bold", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    if title: ax.set_title(title)
    return ax

# ---------------------------------------------------------------- ClinVar classifications
CLINSIG_GROUPS = {                                      # collapse ClinVar labels to 4 buckets
    "pathogenic": ["pathogenic", "likely pathogenic", "pathogenic/likely pathogenic",
                   "pathogenic; risk factor", "likely pathogenic; risk factor"],
    "benign": ["benign", "likely benign", "benign/likely benign"],
    "vus": ["uncertain significance", "uncertain significance; risk factor",
            "conflicting classifications of pathogenicity", "conflicting interpretations of pathogenicity"],
}
def clinsig_group(s):
    """ClinVar label -> pathogenic / benign / vus / other. Compound labels fall back to their
    first term, so 'Pathogenic; drug response' -> pathogenic."""
    s = str(s).strip().lower()
    for key in (s, s.split(";")[0].strip()):
        for grp, labels in CLINSIG_GROUPS.items():
            if key in labels: return grp
    return "other"

# ClinVar review status -> gold stars (0-4), for filtering to well-supported classifications
REVIEW_STARS = {"practice guideline": 4, "reviewed by expert panel": 3,
                "criteria provided, multiple submitters, no conflicts": 2,
                "criteria provided, single submitter": 1,
                "criteria provided, conflicting classifications": 1,
                "criteria provided, conflicting interpretations": 1}
def review_stars(s):
    return REVIEW_STARS.get(str(s).strip().lower(), 0)

def parse_pergene_clinvar(paths, mane, gene_from_name=lambda p: p.split("/")[-1].split("_")[0]):
    """Parse ClinVar per-gene web exports (tab-separated, 'Protein change' + 'Germline
    classification'). Rows whose reference residue disagrees with the MANE sequence are dropped
    -- those are multi-gene CNV entries that carry other genes' protein changes."""
    import glob, re
    rows = []
    for path in sum([glob.glob(p) for p in ([paths] if isinstance(paths, str) else paths)], []):
        gene = gene_from_name(path)
        if gene not in mane: continue
        d = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        d = d[d["Protein change"].notna()]
        for _, r in d.iterrows():
            for pc in str(r["Protein change"]).split(", "):
                m = re.fullmatch(r"([A-Z])(\d+)([A-Z])", pc.strip())
                if not m: continue
                pos, ref, alt = int(m.group(2)), m.group(1), m.group(3)
                seq = mane[gene]
                if pos > len(seq) or seq[pos - 1] != ref:     # not this gene's variant
                    continue
                rows.append({"gene": gene, "position": pos, "reference": ref, "alternate": alt,
                             "clinsig": r["Germline classification"],
                             "review": r.get("Germline review status"),
                             "conditions": r.get("Condition(s)")})
    df = pd.DataFrame(rows)
    return df.drop_duplicates(["gene", "position", "reference", "alternate"])

def parse_variant_summary(path, genes=("NLGN1","NLGN2","NLGN3","NLGN4X","NLGN4Y"), assembly="GRCh38",
                          ids=None):
    """Parse ClinVar's variant_summary.txt(.gz) (full release). Returns rows keyed by the same
    'chrom-pos-ref-alt' id the parquet uses. Pass ids (e.g. set(missense.id)) to keep exactly
    those variants instead of filtering by gene symbol."""
    cols = ["GeneSymbol", "ClinicalSignificance", "ReviewStatus", "PhenotypeList", "Assembly",
            "Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF"]
    ids = set(ids) if ids is not None else None
    out = []
    for chunk in pd.read_csv(path, sep="\t", dtype=str, chunksize=500_000, usecols=cols, low_memory=False):
        chunk = chunk[chunk.Assembly == assembly]
        chunk = chunk.assign(id=chunk.Chromosome + "-" + chunk.PositionVCF + "-"
                                + chunk.ReferenceAlleleVCF + "-" + chunk.AlternateAlleleVCF)
        c = chunk[chunk.id.isin(ids)] if ids is not None else chunk[chunk.GeneSymbol.isin(genes)]
        if len(c): out.append(c)
    d = pd.concat(out, ignore_index=True)
    d["stars"] = d.ReviewStatus.map(review_stars)
    return d.rename(columns={"GeneSymbol": "gene", "ClinicalSignificance": "clinsig",
                             "ReviewStatus": "review", "PhenotypeList": "conditions"})[
        ["gene", "id", "clinsig", "review", "stars", "conditions"]].drop_duplicates("id")

def add_clinsig(missense, pergene=None, summary=None):
    """Attach ClinVar classifications to the missense table. `pergene` joins on
    (gene, position, ref, alt); `summary` joins on the genomic id. Adds clinsig,
    clinsig_group, clinsig_source (plus review, stars from `summary`); unmatched rows get
    NaN / 'unclassified'. Labels attach by variant, so gnomAD rows that are also in ClinVar
    get them too -- the `database` column still says where each row came from."""
    df = missense.copy()
    df["clinsig"] = pd.Series(pd.NA, index=df.index, dtype="object")
    df["clinsig_source"] = pd.Series(pd.NA, index=df.index, dtype="object")
    if summary is not None:
        s = summary.set_index("id")
        m = df.id.map(s.clinsig)
        df.loc[m.notna(), "clinsig"] = m[m.notna()]
        df.loc[m.notna(), "clinsig_source"] = "variant_summary"
        df["review"] = df.id.map(s.review)
        df["stars"] = df.id.map(s.stars).astype("Int64")
    if pergene is not None:
        key = ["gene", "position", "reference", "alternate"]
        m = df.merge(pergene[key + ["clinsig"]], on=key, how="left", suffixes=("", "_pg"))
        fill = df.clinsig.isna() & m.clinsig_pg.notna().values
        df.loc[fill, "clinsig"] = m.clinsig_pg.values[fill]
        df.loc[fill, "clinsig_source"] = "pergene_export"
    df["clinsig_group"] = np.where(df.clinsig.isna(), "unclassified",
                                   [clinsig_group(s) for s in df.clinsig])
    return df

# ---------------------------------------------------------------- gnomAD allele frequencies
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
# the API caps each request at cost 25 and one gene's variant list costs 11 -> 2 genes/request
_GENE_QUERY = ('g{i}: gene(gene_symbol: "{gene}", reference_genome: GRCh38) '
               '{{ variants(dataset: {dataset}) {{ variant_id joint {{ ac an }} }} }}')

def _gnomad_post(genes, dataset, tries=6, backoff=5):
    """One request for up to 2 genes -> {gene: [[variant_id, ac, an], ...] or None if not found}.
    Retries on rate limiting (429), server errors and timeouts."""
    import requests
    query = "{" + " ".join(_GENE_QUERY.format(i=i, gene=g, dataset=dataset)
                           for i, g in enumerate(genes)) + "}"
    for attempt in range(tries):
        try:
            r = requests.post(GNOMAD_API, json={"query": query}, timeout=300)
        except requests.RequestException:
            time.sleep(backoff * 2 ** attempt); continue
        if r.status_code == 429 or r.status_code >= 500:
            time.sleep(int(r.headers.get("Retry-After", backoff * 2 ** attempt))); continue
        r.raise_for_status()
        data = r.json().get("data") or {}
        out = {}
        for i, g in enumerate(genes):
            hit = data.get(f"g{i}")
            out[g] = None if hit is None else [[v["variant_id"], v["joint"]["ac"], v["joint"]["an"]]
                                               for v in hit["variants"] if v["joint"]]
        return out
    raise RuntimeError(f"gnomAD API failed {tries} times for {genes}")

def fetch_gnomad_counts(genes, dataset="gnomad_r4", cache_dir="gnomad_cache", pause=0.5,
                        verbose=True):
    """Joint (exome + genome) allele count and allele number for every gnomAD variant in each
    gene, cached as one JSON per gene so an interrupted run picks up where it stopped.
    Returns (DataFrame[id, gene, gnomad_ac, gnomad_an], [genes the API could not find])."""
    cache = pkg_path(cache_dir)
    os.makedirs(cache, exist_ok=True)
    path = lambda g: os.path.join(cache, f"{g}_{dataset}.json")
    todo = [g for g in dict.fromkeys(genes) if not os.path.exists(path(g))]
    if verbose: print(f"{len(genes) - len(todo)} genes cached, {len(todo)} to fetch")
    missing = []
    for n in range(0, len(todo), 2):
        for g, rows in _gnomad_post(todo[n:n + 2], dataset).items():
            if rows is None: missing.append(g); continue      # not cached, so retried next run
            json.dump(rows, open(path(g), "w"))
        if verbose and (n // 2) % 25 == 0:
            print(f"{n + 2}/{len(todo)} genes fetched")
        time.sleep(pause)
    frames = []
    for g in dict.fromkeys(genes):
        if os.path.exists(path(g)):
            f = pd.DataFrame(json.load(open(path(g))), columns=["id", "gnomad_ac", "gnomad_an"])
            frames.append(f.assign(gene=g))
    counts = pd.concat(frames, ignore_index=True).drop_duplicates("id")
    if verbose and missing: print(f"not found in gnomAD: {missing}")
    return counts[["id", "gene", "gnomad_ac", "gnomad_an"]], missing

def add_allele_freq(missense, counts, common=0.01):
    """Attach gnomAD AC/AN and allele frequency of the ALTERNATE allele by variant id.
    Adds gnomad_ac, gnomad_an, allele_freq, af_percent, common (AF > `common`) and
    ref_is_minor (AF > 1 - `common`: the reference residue in the sequence is the rare one).
    ClinVar rows that are also in gnomAD get frequencies too; ids absent from gnomAD get NaN."""
    df = missense.copy()
    c = counts.set_index("id")
    df["gnomad_ac"] = df.id.map(c.gnomad_ac).astype("Int64")
    df["gnomad_an"] = df.id.map(c.gnomad_an).astype("Int64")
    df["allele_freq"] = (df.gnomad_ac / df.gnomad_an.where(df.gnomad_an > 0)).astype(float)
    df["af_percent"] = 100 * df.allele_freq
    df["common"] = df.allele_freq > common
    df["ref_is_minor"] = df.allele_freq > 1 - common
    return df

# ---------------------------------------------------------------- region definitions
import re as _re
TM_MOTIF = _re.compile(r"DYSTE(?=LSVT)")               # ECD ends here; TM helix follows
TM_LENGTH = 22

def sequence_regions(seq, tm_motif=TM_MOTIF, tm_length=TM_LENGTH):
    """{'TM': set(resnums), 'ICD': set(resnums)} (1-based) from the TM motif. Empty if not found."""
    m = tm_motif.search(seq)
    if not m: return {"TM": set(), "ICD": set()}
    tm0 = m.end()                                       # 0-based index of the first TM residue
    return {"TM": set(range(tm0 + 1, tm0 + tm_length + 1)),
            "ICD": set(range(tm0 + tm_length + 1, len(seq) + 1))}

def region_membership(results, gene, mane, transfer_from=None, classes=("core", "surface"),
                      ifaces=("iface_NRXN1B", "iface_MDGA1", "iface_dimer")):
    """{region: set of MANE residue numbers} for one gene.

    Structural classes and interfaces come from `results` (crystal structures, UniProt numbering,
    mapped to MANE). Regions a gene has no structure for can be filled from a donor paralog via
    `transfer_from=(donor_gene, donor_results)` -- residues are transferred by aligning the two
    MANE sequences. Transferred regions are suffixed '*' to mark them as inferred.
    TM and ICD always come from the sequence motif.
    """
    seq = mane[gene]
    regions = sequence_regions(seq)
    if gene in results:
        t = residue_table(results, gene, seq)           # index = MANE resnum
        for c in classes:
            regions[c] = set(t.index[t["class"] == c])
        for f in ifaces:
            if f in t.columns and t[f].any():
                regions[f.replace("iface_", "")] = set(t.index[t[f]])
    if transfer_from:
        donor, dres = transfer_from
        dseq = mane[donor]
        m = map_all_positions(dseq, seq)                # donor resnum -> this gene's resnum
        dreg = region_membership(dres, donor, mane)
        for k, v in dreg.items():
            if k in ("TM", "ICD") or k in regions: continue
            mapped = {m[i] for i in v if i in m}
            if mapped: regions[k + "*"] = mapped
    return regions

REGION_ORDER = ["core", "surface", "NRXN1B", "MDGA1", "dimer", "TM", "ICD"]

def region_density_table(results, missense, mane=None, genes=("NLGN1","NLGN2","NLGN3","NLGN4X","NLGN4Y"),
                         transfers=None, source="gnomAD", order=REGION_ORDER, long=False,
                         merge_transferred=True):
    """Missense density per region per gene.

    transfers: {gene: donor_gene} -- fill missing structural regions from a paralog by homology
      (marked with '*' in the output).
    source: 'gnomAD', 'ClinVar', or None for both sources combined.
    Returns a wide table (genes x regions) of variants per residue; long=True gives the full
    table with counts, residue counts, and BLOSUM62-weighted severity.
    """
    mane = mane or mane_sequences(list(genes))
    transfers = transfers or {}
    rows = []
    for g in genes:
        if g not in mane: continue
        tf = (transfers[g], results) if g in transfers and transfers[g] in results else None
        regions = region_membership(results, g, mane, transfer_from=tf)
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
    if merge_transferred:                               # one column per region; '*' tracked separately
        df["inferred"] = df.region.str.endswith("*")
        df["region"] = df.region.str.rstrip("*")
    def rank(r):                                        # keep the requested order
        base = r.rstrip("*")
        return (order.index(base) if base in order else len(order), r.endswith("*"))
    cols = sorted(df.region.unique(), key=rank)
    wide = df.pivot(index="gene", columns="region", values="density").reindex(index=list(genes), columns=cols)
    if merge_transferred:
        flags = df.pivot(index="gene", columns="region", values="inferred").reindex(index=list(genes), columns=cols)
        return wide.round(3), flags.notna() & flags.eq(True)
    return wide.round(3)
