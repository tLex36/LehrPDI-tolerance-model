"""How much of the human sequence is conserved in each species, for five kinds of sequence.

Protein columns (Ensembl Compara orthologues):
    protein_similarity (REST API): % of the human protein's residues identical or similar (positive BLOSUM62
    score) in its orthologue (Compara's perc_pos; the best orthologue when there are several).
    protein_identity (BioMart, one table per species): % identical (Compara's perc_id).
    A human gene with no orthologue counts as 0.
DNA columns (Ensembl whole-genome alignments, REST API):
    share of human bases that are identical in the species' aligned sequence, from Ensembl's pairwise
    LASTZ_NET alignment of that species with human. Unaligned bases count as 0. Species for which Ensembl
    has no pairwise alignment with human are left blank (NaN).

Every answer is cached under cache_dir, so only the first run is slow (~15 min).

    from evolution_conservation import SPECIES, protein_identity, sample_regions, dna_identity, summarise
"""
import os, json, re, time, random, threading
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
import requests

SERVER = "https://rest.ensembl.org"
BIOMART = "https://useast.ensembl.org/biomart/martservice"     # same release (116) as the REST server
# tree tip name (species_tree.TREE) -> Ensembl name
SPECIES = {
    "Chimpanzee": "pan_troglodytes", "Rhesus macaque": "macaca_mulatta", "Marmoset": "callithrix_jacchus",
    "Mouse": "mus_musculus", "Rat": "rattus_norvegicus", "Guinea pig": "cavia_porcellus",
    "Dog": "canis_lupus_familiaris", "Horse": "equus_caballus", "Pig": "sus_scrofa", "Cow": "bos_taurus",
    "Dolphin": "tursiops_truncatus", "Elephant": "loxodonta_africana", "Chicken": "gallus_gallus",
    "Xenopus tropicalis": "xenopus_tropicalis", "Zebrafish": "danio_rerio",
    "Drosophila": "drosophila_melanogaster", "C. elegans": "caenorhabditis_elegans",
}
CHROM_LENGTHS = {  # GRCh38
    "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555, "5": 181538259, "6": 170805979,
    "7": 159345973, "8": 145138636, "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
    "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345, "17": 83257441, "18": 80373285,
    "19": 58617616, "20": 64444167, "21": 46709983, "22": 50818468, "X": 156040895}
NONCODING_TYPES = {"enhancer", "CTCF_binding_site", "TF_binding_site", "open_chromatin_region"}

_lock = threading.Lock()
_last = [0.0]


def get(path, params=None, tries=6):
    """GET from the Ensembl REST API as JSON, at most ~12 requests/s, retrying politely."""
    for k in range(tries):
        with _lock:                                          # rate limit shared by all threads
            wait = _last[0] + 1 / 12 - time.time()
            if wait > 0:
                time.sleep(wait)
            _last[0] = time.time()
        try:
            r = requests.get(SERVER + path, params=params, headers={"Content-Type": "application/json"}, timeout=60)
        except requests.RequestException:
            time.sleep(2 ** k); continue
        if r.status_code == 429:                             # too many requests: wait as told
            time.sleep(float(r.headers.get("Retry-After", 2)) + 0.5); continue
        if r.status_code in (400, 404):                      # e.g. "no alignment available"
            try:
                return {"error": r.json().get("error", r.text)}
            except ValueError:
                return {"error": r.text}
        if r.ok:
            return r.json()
        time.sleep(2 ** k)
    raise RuntimeError(f"Ensembl request failed: {path}")


def _cached(path, fetch):
    if os.path.exists(path):
        return json.load(open(path))
    value = fetch()
    json.dump(value, open(path, "w"))
    return value


# ------------------------------------------------------------------------------------------------ proteins
def _biomart_prefix(species):
    """BioMart's name for a species' homolog attributes: canis_lupus_familiaris -> clfamiliaris."""
    *first, last = species.split("_")
    return "".join(w[0] for w in first) + last


def _orthologue_table(species, cache_dir, tries=5):
    """One row per (human gene, orthologue in `species`): gene, % identity, homology type (cached TSV).
    Human genes without an orthologue appear once with empty fields."""
    path = os.path.join(cache_dir, f"{species}.tsv")
    if not os.path.exists(path):
        p = _biomart_prefix(species)
        query = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" '
                 'formatter="TSV" header="0" uniqueRows="0" completionStamp="1" datasetConfigVersion="0.6">'
                 '<Dataset name="hsapiens_gene_ensembl" interface="default"><Attribute name="ensembl_gene_id"/>'
                 f'<Attribute name="{p}_homolog_perc_id"/><Attribute name="{p}_homolog_orthology_type"/>'
                 '</Dataset></Query>')
        for k in range(tries):
            try:
                r = requests.post(BIOMART, data={"query": query}, timeout=600)
                text = r.text.rstrip()
                if r.ok and text.endswith("[success]"):    # BioMart marks a complete download
                    break
            except requests.RequestException:
                pass
            time.sleep(10 * (k + 1))
        else:
            raise RuntimeError(f"BioMart download failed: {species}")
        open(path, "w").write(text[:-len("[success]")])
    return pd.read_csv(path, sep="\t", header=None, names=["gene", "perc_id", "type"])


def protein_identity(ensembl_genes, cache_dir):
    """DataFrame: one row per human gene, one column per species, % identity (0 = no orthologue).
    Genes Ensembl doesn't recognise are dropped."""
    os.makedirs(cache_dir, exist_ok=True)
    ids = list(dict.fromkeys(g.split(".")[0] for g in ensembl_genes))     # unique, order kept
    cols, known = {}, None
    for sp in SPECIES.values():
        t = _orthologue_table(sp, cache_dir)
        known = set(t.gene) if known is None else known & set(t.gene)
        cols[sp] = t[t.type.fillna("").str.startswith("ortholog")].groupby("gene").perc_id.max()
    ids = [g for g in ids if g in known]
    return pd.DataFrame({sp: c.reindex(ids).fillna(0.0) for sp, c in cols.items()}, index=ids)


def _orthologue_similarity(ensg, cache_dir):
    """{ensembl species: % of the human protein's residues identical OR similar (positive BLOSUM62) in its
    best orthologue} for one human gene, from the REST API (Compara's perc_pos; BioMart doesn't provide it)."""
    path = os.path.join(cache_dir, f"{ensg}.json")
    if os.path.exists(path):
        return json.load(open(path))
    d = get(f"/homology/id/human/{ensg}", {"type": "orthologues", "sequence": "none", "aligned": "0",
                                            "target_species": list(SPECIES.values())})
    if "error" in d or not d.get("data"):                       # Ensembl doesn't know this gene
        best = {"_missing": True}
    else:
        best = {}
        for h in d["data"][0]["homologies"]:
            sp = h["target"]["species"]
            best[sp] = max(best.get(sp, 0.0), float(h["source"]["perc_pos"]))
    json.dump(best, open(path, "w"))
    return best


def protein_similarity(ensembl_genes, cache_dir, threads=8):
    """DataFrame: one row per human gene, one column per species, % similarity (0 = no orthologue).
    Genes Ensembl doesn't recognise are dropped; a gene whose requests keep failing is left out (and not
    cached), so running again fills it in."""
    os.makedirs(cache_dir, exist_ok=True)
    ids = list(dict.fromkeys(g.split(".")[0] for g in ensembl_genes))

    def one(g):
        try:
            return _orthologue_similarity(g, cache_dir)
        except RuntimeError as err:
            print("skipped:", err, flush=True)
            return None
    with ThreadPoolExecutor(threads) as pool:
        results = list(pool.map(one, ids))
    failed = sum(r is None for r in results)
    if failed:
        print(f"{failed} of {len(ids)} genes failed; run again to fill them in", flush=True)
    rows = {g: r for g, r in zip(ids, results) if r is not None and not r.get("_missing")}
    return pd.DataFrame([{sp: r.get(sp, 0.0) for sp in SPECIES.values()} for r in rows.values()], index=list(rows))


# ------------------------------------------------------------------------------------------------ DNA
def sample_regions(kind, n, cache_path, seed=0, max_len=1000):
    """n human regions (chrom, start, end), cached. kind = "genome": random 1-kb windows with < 5% N;
    "regulatory": random non-promoter elements of the Ensembl Regulatory Build (enhancers, CTCF / TF
    binding sites, open chromatin), centre-cropped to at most max_len bp."""
    if os.path.exists(cache_path):
        return [tuple(r) for r in json.load(open(cache_path))][:n]
    rng = random.Random(seed)
    chroms, weights = list(CHROM_LENGTHS), list(CHROM_LENGTHS.values())
    regions = []
    while len(regions) < n:
        c = rng.choices(chroms, weights)[0]
        if kind == "genome":
            s = rng.randrange(1, CHROM_LENGTHS[c] - 1000)
            seq = get(f"/sequence/region/human/{c}:{s}..{s + 999}:1").get("seq", "")
            if seq and seq.upper().count("N") < 50:
                regions.append((c, s, s + 999))
        else:
            s = rng.randrange(1, CHROM_LENGTHS[c] - 1_000_000)
            feats = get(f"/overlap/region/human/{c}:{s}-{s + 999_999}", {"feature": "regulatory"})
            feats = [f for f in feats if isinstance(f, dict) and f.get("description") in NONCODING_TYPES] \
                if isinstance(feats, list) else []
            if feats:
                f = rng.choice(feats)
                mid, half = (f["start"] + f["end"]) // 2, min(f["end"] - f["start"] + 1, max_len) // 2
                regions.append((c, mid - half, mid + half - 1))
    json.dump(regions, open(cache_path, "w"))
    return regions


def sample_cds(transcripts, n, cache_path, seed=0, max_len=1000, threads=8):
    """n human coding regions (chrom, start, end), cached: one random coding exon (its CDS part) from each of
    n random transcripts (e.g. MANE Select, one per protein-coding gene), centre-cropped to at most max_len bp.
    Only chr1-22 and X."""
    if os.path.exists(cache_path):
        return [tuple(r) for r in json.load(open(cache_path))][:n]
    rng = random.Random(seed)
    order = rng.sample(list(transcripts), min(len(transcripts), 2 * n))   # spares for genes on other chromosomes

    def one(k_t):
        k, enst = k_t
        enst = enst.split(".")[0]
        feats = get(f"/overlap/id/{enst}", {"feature": "cds"})
        segs = [f for f in feats if isinstance(f, dict) and f.get("Parent") == enst
                and str(f.get("seq_region_name")) in CHROM_LENGTHS] if isinstance(feats, list) else []
        if not segs:
            return None
        f = random.Random(seed * 100003 + k).choice(segs)                 # reproducible per position in `order`
        mid, half = (f["start"] + f["end"]) // 2, min(f["end"] - f["start"] + 1, max_len) // 2
        return (str(f["seq_region_name"]), mid - half, mid + half - 1)
    with ThreadPoolExecutor(threads) as pool:
        found = list(pool.map(one, enumerate(order)))
    regions = [r for r in found if r is not None][:n]
    if len(regions) < n:
        raise RuntimeError(f"only {len(regions)} coding regions found")
    json.dump(regions, open(cache_path, "w"))
    return regions


def _identical_bases(blocks, species):
    """Count human bases identical in `species` over pairwise alignment blocks."""
    same = 0
    for b in blocks:
        al = b.get("alignments", [])
        human = [a["seq"].upper() for a in al if a["species"] == "homo_sapiens"]
        other = [a["seq"].upper() for a in al if a["species"] == species]
        if human and other:
            h = human[0]
            same += max(sum(1 for x, y in zip(h, o) if x != "-" and x == y and x in "ACGT") for o in other)
    return same


def _region_identity(region, species, cache_dir):
    """{ensembl species: share of the region's human bases identical in that species} for one region."""
    c, s, e = region
    length = e - s + 1

    def fetch():
        out = {}
        for sp in species:                                         # pairwise alignment with human
            d = get(f"/alignment/region/homo_sapiens/{c}:{s}-{e}:1",
                    {"method": "LASTZ_NET", "species_set": ["homo_sapiens", sp]})
            if isinstance(d, dict) and "no alignment" not in d.get("error", "").lower():
                raise RuntimeError(f"{c}:{s}-{e} {sp}: {d}")       # a real error: don't cache it as 0
            out[sp] = 0.0 if isinstance(d, dict) else _identical_bases(d, sp) / length
        return out
    return _cached(os.path.join(cache_dir, f"{c}_{s}_{e}.json"), fetch)


def species_with_alignment():
    """The species of SPECIES for which Ensembl has a pairwise (LASTZ_NET) genome alignment with human."""
    lastz = get("/info/compara/species_sets/LASTZ_NET", {"compara": "vertebrates"})
    with_human = {x for s in lastz if "homo_sapiens" in s.get("species_set", []) for x in s["species_set"]}
    return [sp for sp in SPECIES.values() if sp in with_human]


def dna_identity(regions, cache_dir, species, threads=8):
    """DataFrame: one row per region, one column per species, % of human bases identical
    (NaN for species without a pairwise genome alignment to human). A region whose requests keep failing
    is left out (and not cached), so running again fills it in."""
    os.makedirs(cache_dir, exist_ok=True)

    def one(r):
        try:
            return _region_identity(r, species, cache_dir)
        except RuntimeError as err:
            print("skipped:", err, flush=True)
            return None
    with ThreadPoolExecutor(threads) as pool:
        results = list(pool.map(one, regions))
    ok = [(r, x) for r, x in zip(regions, results) if x is not None]
    if len(ok) < len(regions):
        print(f"{len(regions) - len(ok)} of {len(regions)} regions failed; run again to fill them in", flush=True)
    df = pd.DataFrame([x for _, x in ok], index=[f"{c}:{s}-{e}" for (c, s, e), _ in ok])
    return df.reindex(columns=list(SPECIES.values())) * 100           # species with no alignment -> NaN


def summarise(per_item):
    """Mean over genes / regions for each species (columns), as a Series indexed by tree tip names."""
    names = {v: k for k, v in SPECIES.items()}
    return per_item.mean(skipna=False).rename(index=names)
