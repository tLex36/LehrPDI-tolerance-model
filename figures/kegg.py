"""KEGG pathway membership and over-representation testing, cached on disk.

    from kegg import pathway_genes, gene_pathways, enrich
    paths = pathway_genes("mmu")                 # {pathway id: {"name":..., "genes": {...}}}
    hits  = enrich(sig_genes, background=all_genes, organism="mmu")

KEGG's REST API is free for academic use but rate-limited, so every download is cached
under figures/kegg_cache/ and reused.
"""
import json as _json
import os as _os
import re as _re

import pandas as pd

from netutil import pkg_path, urlopen_retry

KEGG_REST = "https://rest.kegg.jp"


def _cached(name, url, cache_dir="kegg_cache"):
    """Fetch a KEGG REST endpoint as text, caching the raw response."""
    d = pkg_path(cache_dir)
    _os.makedirs(d, exist_ok=True)
    path = _os.path.join(d, name)
    if _os.path.exists(path):
        return open(path).read()
    txt = urlopen_retry(url)
    txt = txt.decode() if isinstance(txt, bytes) else txt
    open(path, "w").write(txt)
    return txt


def gene_symbols(organism="mmu") -> dict:
    """{KEGG gene id: primary gene symbol} -- e.g. 'mmu:12064' -> 'Bdnf'."""
    txt = _cached(f"{organism}_genes.tsv", f"{KEGG_REST}/list/{organism}")
    out = {}
    for line in txt.splitlines():
        f = line.split("\t")
        if len(f) < 4:
            continue
        names = f[3].split(";")[0]                 # drop the description
        sym = names.split(",")[0].strip()          # first alias is the primary symbol
        if sym:
            out[f[0]] = sym
    return out


def pathway_names(organism="mmu") -> dict:
    """{pathway id: human-readable name}, organism suffix stripped."""
    txt = _cached(f"{organism}_pathways.tsv", f"{KEGG_REST}/list/pathway/{organism}")
    out = {}
    for line in txt.splitlines():
        f = line.split("\t")
        if len(f) < 2:
            continue
        out[f[0]] = _re.sub(r"\s+-\s+[A-Z][a-z]+ [a-z]+.*$", "", f[1]).strip()
    return out


def pathway_genes(organism="mmu", min_size=5, max_size=500, drop_global=True) -> dict:
    """{pathway id: {"name": str, "genes": set of gene symbols}}.

    min_size/max_size filter out pathways too small to test or so large they are
    uninformative; drop_global removes the KEGG 'global map' entries (01100 Metabolic
    pathways and friends), which contain thousands of genes and dominate any enrichment.
    """
    txt = _cached(f"{organism}_path_links.tsv", f"{KEGG_REST}/link/{organism}/pathway")
    sym = gene_symbols(organism)
    names = pathway_names(organism)
    members = {}
    for line in txt.splitlines():
        f = line.split("\t")
        if len(f) < 2:
            continue
        pid = f[0].replace("path:", "")
        s = sym.get(f[1])
        if s:
            members.setdefault(pid, set()).add(s)
    out = {}
    for pid, genes in members.items():
        if drop_global and pid[3:].startswith("011"):     # 01100/01110/01200... global maps
            continue
        if min_size <= len(genes) <= max_size:
            out[pid] = {"name": names.get(pid, pid), "genes": genes}
    return out


def gene_pathways(organism="mmu", **kw) -> dict:
    """{gene symbol: [pathway ids]} -- the inverse of pathway_genes."""
    out = {}
    for pid, d in pathway_genes(organism, **kw).items():
        for g in d["genes"]:
            out.setdefault(g, []).append(pid)
    return out


def find_pathways(query, organism="mmu", **kw) -> pd.DataFrame:
    """Pathways whose name matches `query` (case-insensitive substring or regex)."""
    rows = [{"pathway": pid, "name": d["name"], "n_genes": len(d["genes"])}
            for pid, d in pathway_genes(organism, **kw).items()
            if _re.search(query, d["name"], _re.I)]
    return pd.DataFrame(rows).sort_values("n_genes", ascending=False).reset_index(drop=True)


def pathway_members(pathway, organism="mmu", **kw) -> set:
    """Gene symbols in one pathway, given its id ('mmu04721') or a name substring."""
    paths = pathway_genes(organism, **kw)
    if pathway in paths:
        return paths[pathway]["genes"]
    hits = [p for p, d in paths.items() if _re.search(pathway, d["name"], _re.I)]
    if not hits:
        raise KeyError(f"no pathway matching {pathway!r}")
    if len(hits) > 1:
        raise KeyError(f"{pathway!r} matches {len(hits)}: "
                       + ", ".join(f"{p} ({paths[p]['name']})" for p in hits[:6]))
    return paths[hits[0]]["genes"]


def enrich(genes, background=None, organism="mmu", min_size=5, max_size=500,
           alpha=0.05, method="fdr_bh") -> pd.DataFrame:
    """Over-representation test of `genes` against every KEGG pathway.

    genes:      your significant gene symbols
    background: the universe to test against -- ALL genes detected in the experiment.
                Defaults to every gene KEGG annotates, which inflates significance; pass
                your own detected-gene list whenever you have one.
    Returns one row per pathway with the hypergeometric p, Benjamini-Hochberg q, and the
    overlapping genes.
    """
    from scipy.stats import hypergeom
    paths = pathway_genes(organism, min_size=min_size, max_size=max_size)
    annotated = set().union(*(d["genes"] for d in paths.values())) if paths else set()

    bg = set(background) & annotated if background is not None else annotated
    hits = set(genes) & bg
    M, N = len(bg), len(hits)                      # universe, drawn
    rows = []
    for pid, d in paths.items():
        inpath = d["genes"] & bg
        n = len(inpath)                            # successes in universe
        k = len(inpath & hits)                     # successes drawn
        if k == 0 or n == 0:
            continue
        p = hypergeom.sf(k - 1, M, n, N)
        rows.append({"pathway": pid, "name": d["name"], "n_pathway": n, "n_overlap": k,
                     "expected": round(n * N / M, 2) if M else float("nan"),
                     "fold_enrichment": round((k / N) / (n / M), 2) if M and N else float("nan"),
                     "p_value": p, "genes": ",".join(sorted(inpath & hits))})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    try:
        from statsmodels.stats.multitest import multipletests
        out["q_value"] = multipletests(out.p_value, method=method)[1]
    except ImportError:                            # Benjamini-Hochberg by hand
        o = out.sort_values("p_value").reset_index(drop=True)
        m = len(o)
        q = (o.p_value * m / (o.index + 1)).cummin()[::-1].cummin()[::-1]
        out = o.assign(q_value=q.clip(upper=1.0))
    out = out.sort_values("p_value").reset_index(drop=True)
    out.attrs["n_query"] = N
    out.attrs["n_background"] = M
    return out
