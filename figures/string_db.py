"""STRING functional enrichment and PPI-network statistics, cached on disk.

    import string_db as sdb
    enr = sdb.enrichment(genes, species="rat")          # all 14 STRING categories
    ppi = sdb.ppi_enrichment(genes, species="rat")      # is the network over-connected?
    net = sdb.network(genes, species="rat")             # scored edges

STRING differs from the KEGG/GO over-representation in `kegg.py` / `go.py` in two ways
that matter for an IP-MS experiment:

  * it tests ~14 annotation sources at once (GO, KEGG, Reactome, Pfam, InterPro,
    WikiPathways, UniProt keywords, PMID co-citation, TISSUES, COMPARTMENTS, and its own
    network-neighbourhood clusters), so it finds structure that a single ontology misses;
  * `ppi_enrichment` asks a question neither can: are these proteins more interconnected
    in prior literature than a random set of the same size? For a pulldown that is the
    single most relevant sanity check -- a real complex is over-connected, a list of
    sticky abundant proteins is not.

Requests are POSTed (a 1300-protein list does not fit in a URL), rate-limited to STRING's
requested 1 req/s, and cached under figures/string_cache/ keyed by payload hash.
"""
import hashlib as _hashlib
import io as _io
import json as _json
import os as _os
import time as _time

import numpy as np
import pandas as pd

from netutil import pkg_path, urlpost_retry

STRING_API = "https://string-db.org/api"

# STRING asks that scripted users identify themselves.
CALLER_IDENTITY = "LehrPDI_figures"

SPECIES = {"rat": 10116, "mouse": 10090, "human": 9606}

_LAST_CALL = [0.0]
_MIN_INTERVAL = 1.0                     # STRING asks for <= 1 request per second


def _species(s):
    return s if isinstance(s, int) else SPECIES[str(s).lower()]


def _post(method, fields, cache_dir="string_cache", refresh=False):
    """POST to a STRING endpoint, caching the raw TSV by payload hash."""
    d = pkg_path(cache_dir)
    _os.makedirs(d, exist_ok=True)
    key = _hashlib.sha1(
        _json.dumps([method, sorted(fields.items())], sort_keys=True).encode()
    ).hexdigest()[:20]
    path = _os.path.join(d, f"{method}_{key}.tsv")
    if _os.path.exists(path) and not refresh:
        return open(path).read()

    wait = _MIN_INTERVAL - (_time.time() - _LAST_CALL[0])
    if wait > 0:
        _time.sleep(wait)
    raw = urlpost_retry(f"{STRING_API}/tsv/{method}", fields)
    _LAST_CALL[0] = _time.time()

    txt = raw.decode() if isinstance(raw, bytes) else raw
    if txt.lstrip().startswith("{") and '"Error"' in txt:
        raise RuntimeError(f"STRING {method}: {txt[:300]}")
    open(path, "w").write(txt)
    return txt


def _frame(txt):
    if not txt.strip():
        return pd.DataFrame()
    return pd.read_csv(_io.StringIO(txt), sep="\t")


def map_ids(genes, species="rat", **kw) -> pd.DataFrame:
    """Resolve gene symbols to STRING protein ids.

    Returns one row per *successfully mapped* query -- unmapped symbols are simply
    absent, so compare len() against the input to see the mapping rate."""
    txt = _post("get_string_ids",
                {"identifiers": "%0d".join(map(str, genes)),
                 "species": _species(species), "limit": 1, "echo_query": 1,
                 "caller_identity": CALLER_IDENTITY}, **kw)
    return _frame(txt)


def _bg_ids(background, species, **kw):
    if background is None:
        return None
    m = map_ids(background, species=species, **kw)
    return "%0d".join(m.stringId.tolist())


def enrichment(genes, species="rat", background=None, categories=None,
               fdr=0.05, **kw) -> pd.DataFrame:
    """Functional enrichment across every STRING category.

    background: gene symbols for a custom universe. Left None, STRING tests against the
            whole genome -- which for a bait-vs-isotype IP is the right question ("what is
            this pulldown?"). Passing the detected proteome instead asks "what is special
            about these relative to everything we saw", which in a clean IP has very
            little contrast left.
    categories: keep only these (e.g. ["KEGG", "Process", "RCTM"]).
    fdr:    drop rows above this FDR; None keeps everything.

    Adds `fold_enrichment` -- STRING reports counts, not the ratio, so it is computed as
    (hits / n_query) / (term_size / background_size).
    """
    fields = {"identifiers": "%0d".join(map(str, genes)),
              "species": _species(species), "caller_identity": CALLER_IDENTITY}
    bg = _bg_ids(background, species, **kw)
    if bg:
        fields["background_string_identifiers"] = bg

    e = _frame(_post("enrichment", fields, **kw))
    if e.empty:
        return e

    n_query = len(set(map(str, genes)))
    bg_total = len(set(map(str, background))) if background is not None else 20000
    e["fold_enrichment"] = ((e.number_of_genes / n_query) /
                            (e.number_of_genes_in_background / bg_total)).round(2)
    e["n_query"] = n_query
    if categories:
        e = e[e.category.isin(categories)]
    if fdr is not None:
        e = e[e.fdr < fdr]
    return e.sort_values("fdr").reset_index(drop=True)


def ppi_enrichment(genes, species="rat", background=None, **kw) -> dict:
    """Is this set more interconnected than a random set of the same size?

    Returns observed vs expected edge counts and STRING's p-value. A significant result
    means the proteins share more prior evidence of interaction than chance -- the
    strongest single argument that a pulldown recovered a complex rather than a pile of
    abundant protein."""
    fields = {"identifiers": "%0d".join(map(str, genes)),
              "species": _species(species), "caller_identity": CALLER_IDENTITY}
    bg = _bg_ids(background, species, **kw)
    if bg:
        fields["background_string_identifiers"] = bg
    f = _frame(_post("ppi_enrichment", fields, **kw))
    if f.empty:
        return {}
    r = f.iloc[0].to_dict()
    exp = r.get("expected_number_of_edges", 0) or np.nan
    r["edge_ratio"] = round(r.get("number_of_edges", 0) / exp, 2) if exp else np.nan
    return r


def network(genes, species="rat", score=400, **kw) -> pd.DataFrame:
    """Scored edges among the query proteins. `score` is STRING's 0-1000 confidence."""
    txt = _post("network",
                {"identifiers": "%0d".join(map(str, genes)),
                 "species": _species(species), "required_score": score,
                 "caller_identity": CALLER_IDENTITY}, **kw)
    return _frame(txt)


# ----------------------------------------------------------------------- charts

# Okabe-Ito: a published colour-vision-deficiency-safe qualitative palette. Assigned to
# categories in fixed order so a category keeps its colour between figures.
CAT_COLORS = {
    "KEGG":              "#0072B2",
    "Process":           "#D55E00",
    "Component":         "#009E73",
    "Function":          "#CC79A7",
    "RCTM":              "#E69F00",
    "WikiPathways":      "#56B4E9",
    "COMPARTMENTS":      "#8c8c8c",
    "TISSUES":           "#7a5195",
    "InterPro":          "#4d4d4d",
    "Pfam":              "#bc5090",
    "SMART":             "#6b8e23",
    "Keyword":           "#a0522d",
    "PMID":              "#999999",
    "NetworkNeighborAL": "#444444",
}

CAT_LABEL = {"RCTM": "Reactome", "Process": "GO biological process",
             "Component": "GO cellular component", "Function": "GO molecular function",
             "NetworkNeighborAL": "STRING cluster"}


def _wrap(s, n=46):
    out, line = [], ""
    for w in str(s).split():
        if len(line) + len(w) + 1 > n:
            out.append(line); line = w
        else:
            line = f"{line} {w}".strip()
    out.append(line)
    return "\n".join(out)


def enrichment_plot(enr, categories=("KEGG", "Process", "Component", "RCTM"),
                    top=8, ax=None, figsize=None, title=None, max_fdr=0.05,
                    max_term_size=600, label_counts=True):
    """Horizontal bars of -log10(FDR), grouped by STRING category.

    One axis only: bar length is significance, the gene count is a direct label rather
    than a second encoding, so nothing is plotted against a hidden second scale.

    max_term_size: drop terms annotated to more genes than this. STRING's roots
            ("Metabolism", "Cellular process", "Immune System") carry thousands of genes
            and reach tiny FDRs in any large list without saying anything; None keeps them.
    """
    import matplotlib.pyplot as plt

    e = enr[enr.fdr < max_fdr] if max_fdr is not None else enr.copy()
    if max_term_size is not None:
        e = e[e.number_of_genes_in_background <= max_term_size]
    keep = [c for c in categories if c in set(e.category)]

    # Build (position, ...) explicitly: category gaps advance y without consuming a tick
    # slot, so labels can never drift out of step with their bars.
    bars, y = [], 0.0
    for c in keep[::-1]:                                    # first category at the top
        sub = e[e.category == c].nsmallest(top, "fdr").iloc[::-1]
        if sub.empty:
            continue
        for _, r in sub.iterrows():
            bars.append(dict(y=y, cat=c, desc=r.description, fdr=r.fdr,
                             n=r.number_of_genes))
            y += 1
        y += 1                                              # gap between categories
    if not bars:
        raise ValueError("nothing passed max_fdr / max_term_size")

    if ax is None:
        h = max(3.0, 0.32 * len(bars) + 1.1 * len(keep))
        _, ax = plt.subplots(figsize=figsize or (9.5, h))

    for b in bars:
        v = -np.log10(max(b["fdr"], 1e-300))
        ax.barh(b["y"], v, height=.68, color=CAT_COLORS.get(b["cat"], "#777"),
                edgecolor="white", linewidth=.8, zorder=3)
        if label_counts:
            ax.text(v + .4, b["y"], f"{int(b['n'])}", va="center", ha="left",
                    fontsize=7.5, color="#444", zorder=4)

    ax.set_yticks([b["y"] for b in bars])
    ax.set_yticklabels([_wrap(b["desc"]) for b in bars], fontsize=8)
    ax.set_ylim(-.8, y - 1)
    ax.set_xlabel("$-\\log_{10}$ FDR")
    ax.axvline(-np.log10(0.05), color="#888", ls=":", lw=1, zorder=2)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#e8e8e8", lw=.8, zorder=0)
    ax.set_axisbelow(True)

    handles = [plt.Rectangle((0, 0), 1, 1, color=CAT_COLORS.get(c, "#777"))
               for c in keep]
    ax.legend(handles, [CAT_LABEL.get(c, c) for c in keep],
              loc="lower right", frameon=False, fontsize=8,
              title="STRING category", title_fontsize=8)
    if label_counts:
        ax.text(1.0, 1.01, "number at bar end = genes in term", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=7.5, color="#666")
    if title:
        ax.set_title(title, loc="left", fontsize=11, weight="bold", pad=14)
    return ax


def compare_plot(a, b, labels=("before", "after"), categories=("KEGG",), top=14,
                 ax=None, figsize=None, title=None, max_fdr=0.05, max_term_size=600):
    """Dumbbell of -log10(FDR) for the same terms under two gene sets.

    Terms are chosen by their best FDR in either run, so a term that only appears in one
    of them is still shown -- which is the point of the comparison. A term present in one
    run only is drawn as a hollow marker on the axis floor, labelled lost/gained."""
    import matplotlib.pyplot as plt

    def idx(e):
        e = e[e.category.isin(categories)]
        if max_term_size is not None:
            e = e[e.number_of_genes_in_background <= max_term_size]
        # key on the term id so identically-named terms from different categories
        # (GO CC and COMPARTMENTS share many names) are never merged
        return {(r.category, r.term): (r.fdr, r.description) for _, r in e.iterrows()}

    A, B = idx(a), idx(b)
    fdr_of = lambda D, k: D[k][0] if k in D else None
    best = lambda k: min(fdr_of(A, k) or 1, fdr_of(B, k) or 1)
    keys = sorted(set(A) | set(B), key=best)
    keys = [k for k in keys if best(k) < max_fdr][:top][::-1]

    if ax is None:
        _, ax = plt.subplots(figsize=figsize or (9, max(3, .36 * len(keys) + 1.4)))

    CA, CB = "#0072B2", "#D55E00"
    nl = lambda v: -np.log10(max(v, 1e-300))
    # a term absent from one run is drawn on the axis floor, which must therefore be
    # visible -- otherwise the hollow "lost"/"gained" marker is clipped off-plot
    floor = 0.0
    # the two runs sit on slightly offset rows: where a term's FDR is unchanged the
    # markers would otherwise coincide exactly and the "before" dot would be hidden
    DY = .16
    for i, k in enumerate(keys):
        va, vb = fdr_of(A, k), fdr_of(B, k)
        xa, xb = (nl(va) if va else floor), (nl(vb) if vb else floor)
        ax.plot([xa, xb], [i + DY, i - DY], color="#c9ccd1", lw=2, zorder=2,
                solid_capstyle="round")
        if va: ax.scatter(xa, i + DY, s=58, color=CA, zorder=3, edgecolor="white", lw=1.2)
        if vb: ax.scatter(xb, i - DY, s=58, color=CB, zorder=3, edgecolor="white", lw=1.2)
        if not vb:
            ax.scatter(floor, i - DY, s=58, facecolor="white", edgecolor=CB, lw=1.4,
                       zorder=3)
            ax.text(floor + .25, i - DY, "lost", fontsize=7, color=CB, va="center",
                    style="italic")
        if not va:
            ax.scatter(floor, i + DY, s=58, facecolor="white", edgecolor=CA, lw=1.4,
                       zorder=3)
            ax.text(floor + .25, i + DY, "gained", fontsize=7, color=CA, va="center",
                    style="italic")

    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels([_wrap((A.get(k) or B.get(k))[1], 52) for k in keys], fontsize=8)
    ax.set_ylim(-.8, len(keys) - .2)
    ax.set_xlim(left=floor - .3)
    ax.set_xlabel("$-\\log_{10}$ FDR")
    ax.axvline(-np.log10(0.05), color="#888", ls=":", lw=1, zorder=1)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#e8e8e8", lw=.8, zorder=0)
    ax.set_axisbelow(True)
    h = [plt.Line2D([], [], marker="o", ls="", color=c, markersize=8,
                    markeredgecolor="white")
         for c in (CA, CB)]
    ax.legend(h, list(labels), loc="lower right", frameon=False, fontsize=8)
    if title:
        ax.set_title(title, loc="left", fontsize=11, weight="bold", pad=14)
    return ax


def network_plot(genes, species="rat", score=700, ax=None, figsize=(9, 9),
                 label=None, highlight=(), title=None, seed=1, **kw):
    """Spring-layout PPI graph of the query proteins.

    Node size is degree *within this subnetwork*, which is a property of the figure's own
    data rather than an external measure -- it reads as "how central here", nothing more.
    """
    import matplotlib.pyplot as plt
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("network_plot needs networkx:  pip install networkx")

    e = network(genes, species=species, score=score, **kw)
    if e.empty:
        raise ValueError(f"no edges at score >= {score}")
    G = nx.Graph()
    G.add_nodes_from(set(e.preferredName_A) | set(e.preferredName_B))
    for _, r in e.iterrows():
        G.add_edge(r.preferredName_A, r.preferredName_B, weight=float(r.score))

    # a highlighted protein with no edge above the cutoff is simply not in the graph
    hi = set(highlight) & set(G)
    dropped_hi = sorted(set(highlight) - set(G))

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    pos = nx.spring_layout(G, seed=seed, k=.55, weight="weight", iterations=180)
    deg = dict(G.degree())

    for u, v, dd in G.edges(data=True):
        ax.plot(*zip(pos[u], pos[v]), color="#d6d9dd", lw=.4 + 2.0 * (dd["weight"] - .4),
                zorder=1, solid_capstyle="round", alpha=.75)
    xs = [pos[n][0] for n in G]; ys = [pos[n][1] for n in G]
    sz = [40 + 26 * deg[n] for n in G]
    co = ["#D55E00" if n in hi else "#0072B2" for n in G]
    ax.scatter(xs, ys, s=sz, c=co, zorder=2, edgecolor="white", linewidth=1.1)

    show = (set(label) if label is not None else {
        n for n in G if deg[n] >= sorted(deg.values())[-min(len(G), 30)]} | hi) & set(G)
    for n in show:
        ax.annotate(n, pos[n], fontsize=7.5, ha="center", va="center",
                    xytext=(0, 9), textcoords="offset points", zorder=3,
                    color="#111" if n in hi else "#333",
                    weight="bold" if n in hi else "normal")
    ax.set_axis_off()
    if title:
        ax.set_title(title, loc="left", fontsize=11, weight="bold")
    note = (f"{G.number_of_nodes()} nodes · {G.number_of_edges()} edges · "
            f"STRING score ≥ {score/1000:.2f}")
    if dropped_hi:
        note += f"\nno edge at this cutoff: {', '.join(dropped_hi)}"
    ax.text(.01, .01, note, transform=ax.transAxes, fontsize=8, color="#666",
            va="bottom")
    return ax, G


# ------------------------------------------------------- annotated module network

# Functional classes, ordered. Modules are coloured by the class their auto-name (or
# override) matches, so a figure's colours mean something rather than just separating
# blobs: saturated hues for the synaptic signal, muted earth tones for the
# co-purifying background that every brain IP carries.
# Checked in order, first match wins, so the specific classes come before the generic
# ones -- "Oxidative phosphorylation" must reach `mitochondrial` before `signalling`
# can claim it on "phosphoryl".
MODULE_CLASSES = [
    ("synaptic",    "#0072B2", ("synap", "nmda", "glutamat", "postsynap", "presynap",
                                "neurotransmitter", "dendrit", "receptor activation",
                                "receptor complex", "active zone", "psd ", "scaffold",
                                "gaba-a", "gabaergic")),
    ("mitochondrial", "#7a6a56", ("mitochond", "oxidative phosphoryl", "respirat",
                                  "beta-oxidation", "electron transport", "tca cycle",
                                  "citrate cycle")),
    ("RNA & translation", "#5B7FA6", ("ribosom", "translation", "splicing", "hnrnp",
                                      "trna", "rna ", "spliceosome", "mrna")),
    ("ER / folding", "#8c6bb1", ("endoplasmic", "chaperone", "folding", "glycosyl",
                                 "lumen", "lysosom", "melanosome", "unfolded")),
    ("trafficking", "#56B4E9", ("vesicle", "endocyt", "exocyt", "exocyst", "transport",
                                "motor", "microtubule", "cytoskelet", "actin", "clathrin",
                                "snare", "rab ", "endosom", "arf ", "recycling", "dynein",
                                "kinesin", "septin")),
    ("metabolic",   "#B07A3C", ("metabol", "glycoly", "carbon", "nucleotide", "amino acid",
                                "catabol", "biosynth", "oxidoreduct")),
    ("adhesion",    "#009E73", ("adhesion", "junction", "axonogenesis", "growth cone",
                                "projection", "myelin", "cell surface", "axon")),
    ("channels",    "#D55E00", ("channel", "potassium", "transporter", "pump", "atpase")),
    ("signalling",  "#CC79A7", ("gtpase", "signal", "kinase", "calcium", "g alpha",
                                "gaba receptor", "camp", "second messenger", "cyclase",
                                "gpcr", "pka")),
]
_FALLBACK_COLOR = "#9aa0a6"


def classify_module(name):
    """Map a module's name onto a functional class + colour. Falls back to grey."""
    low = str(name).lower()
    for cls, colour, keys in MODULE_CLASSES:
        if any(k in low for k in keys):
            return cls, colour
    return "other", _FALLBACK_COLOR


def name_modules(modules, species="rat", categories=("Process", "KEGG", "Component", "RCTM"),
                 max_term_size=400, **kw):
    """Auto-label each module with the term that best covers it.

    Scored by (fraction of the module annotated) x log2(fold enrichment), which prefers a
    term that describes most of the module over a very specific one matching two genes.
    Returns {index: name}; override anything that comes out wrong.
    """
    out = {}
    for i, genes in enumerate(modules):
        genes = sorted(genes)
        if len(genes) < 3:
            out[i] = f"({len(genes)} proteins)"
            continue
        try:
            e = enrichment(genes, species=species, fdr=0.05, **kw)
        except Exception:
            out[i] = f"module {i+1}"; continue
        x = e[e.category.isin(categories) &
              (e.number_of_genes_in_background <= max_term_size)].copy()
        if x.empty:
            x = e[e.category.isin(categories)].copy()
        if x.empty:
            out[i] = f"module {i+1}"; continue
        x["_cover"] = x["number_of_genes"] / len(genes)
        x["_score"] = x["_cover"] * np.log2(x["fold_enrichment"].clip(lower=1.01))
        out[i] = str(x.nlargest(1, "_score").iloc[0]["description"])
    return out


def _blob(pts, pad=.09, n=240):
    """Smooth closed outline enclosing `pts` (Nx2), padded outwards."""
    from scipy.spatial import ConvexHull
    from scipy.interpolate import splprep, splev
    P = np.asarray(pts, float)
    c = P.mean(0)
    if len(P) == 1:
        th = np.linspace(0, 2 * np.pi, n)
        return c + pad * np.c_[np.cos(th), np.sin(th)]
    if len(P) <= 5:
        # a hull over 3-5 points is a hard-edged triangle or quad; a circle enclosing
        # them reads as the same kind of object as the big modules' smooth blobs
        r = np.linalg.norm(P - c, axis=1).max() + pad
        th = np.linspace(0, 2 * np.pi, n)
        return c + r * np.c_[np.cos(th), np.sin(th)]
    if len(P) == 2:                                   # capsule around a pair
        d = P[1] - P[0]; L = np.hypot(*d) or 1.0
        nrm = np.array([-d[1], d[0]]) / L
        P = np.vstack([P[0] + nrm * pad, P[1] + nrm * pad,
                       P[1] - nrm * pad, P[0] - nrm * pad])
    try:
        h = ConvexHull(P); V = P[h.vertices]
    except Exception:                                 # collinear points
        V = P
    r = V - c
    norm = np.linalg.norm(r, axis=1, keepdims=True)
    V = c + r * (1 + pad / np.clip(norm, 1e-6, None))  # push each vertex out by `pad`
    if len(V) < 4:
        return np.vstack([V, V[:1]])
    try:
        V = np.vstack([V, V[:1]])
        tck, _ = splprep([V[:, 0], V[:, 1]], s=0, per=True, k=3)
        u = np.linspace(0, 1, n)
        return np.c_[splev(u, tck)].T
    except Exception:
        return np.vstack([V, V[:1]])


def _module_layout(G, groups, seed=7, margin=.055, iters=900):
    """Two-stage layout: modules placed as disks, members arranged inside their own.

    A single spring layout over the whole graph cannot do this -- disconnected
    components drift to arbitrary corners and dense modules collapse into one blob, so
    the outlines that are meant to annotate the modules end up overlapping instead.
    Here each module gets a disk sized by sqrt(members); the disks are pushed apart
    until they no longer intersect, then each module is laid out inside its own.
    """
    import networkx as nx
    idx = {g: i for i, grp in enumerate(groups) for g in grp}

    M = nx.Graph(); M.add_nodes_from(range(len(groups)))
    for u, v in G.edges():
        a, b = idx.get(u), idx.get(v)
        if a is None or b is None or a == b:
            continue
        if M.has_edge(a, b):
            M[a][b]["weight"] += 1
        else:
            M.add_edge(a, b, weight=1)

    n = len(groups)
    R = np.array([.052 * np.sqrt(len(groups[i])) + .045 for i in range(n)])

    # The module graph is usually disconnected (a module with no inter-module edge is
    # a separate component), and spring_layout throws disconnected components far
    # apart. Weak edges between every pair keep the whole thing in one frame; the real
    # inter-module edge counts still dominate the attraction.
    W = M.copy()
    for i in range(n):
        for j in range(i + 1, n):
            if not W.has_edge(i, j):
                W.add_edge(i, j, weight=.02)
    c = nx.spring_layout(W, seed=seed, weight="weight", iterations=400,
                         k=2.4 / max(1.0, np.sqrt(n)))
    C = np.array([c[i] for i in range(n)], float)
    span = max(np.ptp(C, axis=0).max(), 1e-6)
    C *= 2.6 * R.sum() / span

    def separate(passes=14):
        for _ in range(passes):
            moved = False
            for i in range(n):
                for j in range(i + 1, n):
                    d = C[j] - C[i]
                    dist = float(np.hypot(*d))
                    need = R[i] + R[j] + margin
                    if dist < need:
                        moved = True
                        if dist < 1e-9:
                            d = np.array([np.cos(i * 2.4), np.sin(i * 2.4)]); dist = 1.0
                        push = (need - dist) / 2 * (d / dist)
                        C[i] -= push; C[j] += push
            if not moved:
                return True
        return False

    # circle packing: pull everything toward the centre, then resolve the collisions
    # that creates. Alternating the two settles into a compact, non-overlapping layout
    # no matter how far apart the spring layout started them.
    for _ in range(iters // 6):
        mid = C.mean(0)
        C += (mid - C) * .06
        for i in range(n):                           # keep real neighbours close
            nb = list(M.neighbors(i)) if i in M else []
            if nb:
                tgt = C[nb].mean(0)
                C[i] += (tgt - C[i]) * .03
        separate()
    separate(passes=400)

    pos = {}
    for i, grp in enumerate(groups):
        grp = sorted(grp)
        sub = G.subgraph(grp)
        if len(grp) == 1:
            local = {grp[0]: np.zeros(2)}
        else:
            local = nx.spring_layout(sub, seed=seed, iterations=400,
                                     k=1.1 / np.sqrt(len(grp)))
        P = np.array([local[g] for g in grp], float)
        P -= P.mean(0)
        scale = np.linalg.norm(P, axis=1).max() or 1.0
        P = P / scale * R[i] * .74                   # leave room inside the outline
        for g, p in zip(grp, P):
            pos[g] = C[i] + p
    return pos


def _declutter(xy, texts, ax, fontsize, iters=2500, pad=4.0, anchors=None, pull=.012):
    """Nudge label anchors apart in display space; returns adjusted display coords.

    Where a dense module holds more labels than fit inside it, they have to end up
    outside its outline, so separation runs until the boxes genuinely clear rather than
    for a fixed small number of passes. A weak pull back toward each label's own node
    keeps them from wandering further than the leader line can explain.
    """
    xy = np.asarray(xy, float).copy()
    home = xy.copy() if anchors is None else np.asarray(anchors, float)
    # deliberately over-estimate: these feed a collision test, and a box reserved
    # slightly too small leaves a few-pixel overlap in the rendered figure
    w = np.array([max(len(t), 3) * fontsize * 0.78 + pad * 2 for t in texts])
    h = np.full(len(texts), fontsize * 1.6 + pad * 2)
    n = len(xy)
    for it in range(iters):
        moved = False
        for i in range(n):
            for j in range(i + 1, n):
                dx = xy[j, 0] - xy[i, 0]; dy = xy[j, 1] - xy[i, 1]
                ox = (w[i] + w[j]) / 2 - abs(dx)
                oy = (h[i] + h[j]) / 2 - abs(dy)
                if ox > 0 and oy > 0:                  # boxes overlap
                    moved = True
                    if oy * 2.2 < ox:                  # cheaper to separate vertically
                        s = (oy / 2 + .6) * (1 if dy >= 0 else -1)
                        xy[i, 1] -= s; xy[j, 1] += s
                    else:
                        s = (ox / 2 + .6) * (1 if dx >= 0 else -1)
                        xy[i, 0] -= s; xy[j, 0] += s
        if pull:                                        # drift back toward the node
            xy += (home - xy) * pull
        if not moved and it > 4:
            break
    return xy


def module_network(genes, species="rat", score=700, ax=None, figsize=(13, 11),
                   resolution=1.1, min_module=3, min_component=3, seed=7,
                   names=None, colors=None, highlight=(), label=None,
                   label_top=4, fontsize=7.6, title=None, legend=True,
                   hull_pad=.09, annotate_hulls=True, verbose=True, **kw):
    """PPI network with Louvain modules outlined, named and coloured by function.

    names:  {module index: label} overriding the auto-generated name. Print the returned
            frame once, then pin the ones that came out wrong -- auto-naming is a
            starting point, not an answer.
    colors: {module label: colour} overriding the functional-class colour.
    label:  genes to label; default is the `highlight` set plus the top `label_top`
            by degree in each module.

    Returns (ax, G, modules) where `modules` is a frame of module, class, colour, size
    and members.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("module_network needs networkx:  pip install networkx")

    e = network(genes, species=species, score=score, **kw)
    if e.empty:
        raise ValueError(f"no edges at score >= {score}")
    G = nx.Graph()
    for _, r in e.iterrows():
        G.add_edge(r.preferredName_A, r.preferredName_B, weight=float(r.score))

    # drop tiny fragments: a 2-node island carries no module structure and only adds
    # visual noise at the edge of the canvas
    keep = {n for c in nx.connected_components(G) if len(c) >= min_component for n in c}
    G = G.subgraph(keep).copy()

    comm = nx.community.louvain_communities(G, weight="weight", seed=seed,
                                            resolution=resolution)
    comm = sorted((set(c) for c in comm), key=len, reverse=True)
    big = [c for c in comm if len(c) >= min_module]
    rest = set().union(*[c for c in comm if len(c) < min_module]) if \
        any(len(c) < min_module for c in comm) else set()

    auto = name_modules(big, species=species)
    labels = [(names or {}).get(i, auto.get(i, f"module {i+1}")) for i in range(len(big))]
    cls_col = [classify_module(l) for l in labels]
    cols = [(colors or {}).get(l, c) for l, (_, c) in zip(labels, cls_col)]

    mod_of = {g: i for i, c in enumerate(big) for g in c}

    # every node must belong to some group for the layout, so the sub-threshold
    # communities are laid out as their own (unlabelled, grey) groups
    small = [c for c in comm if len(c) < min_module]
    pos = _module_layout(G, big + small, seed=seed)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    for u, v in G.edges():
        same = mod_of.get(u, -1) == mod_of.get(v, -2)
        ax.plot(*zip(pos[u], pos[v]), zorder=2,
                color=cols[mod_of[u]] if same else "#c8ccd0",
                lw=.9 if same else .5, alpha=.5 if same else .3,
                solid_capstyle="round")

    for i, c in enumerate(big):
        P = _blob([pos[g] for g in c], pad=hull_pad)
        ax.add_patch(Polygon(P, closed=True, facecolor=cols[i], alpha=.10,
                             edgecolor=cols[i], lw=1.3, zorder=1, joinstyle="round"))

    deg = dict(G.degree())
    for i, c in enumerate(big):
        pts = np.array([pos[g] for g in c])
        ax.scatter(pts[:, 0], pts[:, 1], s=[34 + 20 * deg[g] for g in c],
                   c=cols[i], zorder=3, edgecolor="white", linewidth=1.1)
    if rest:
        pts = np.array([pos[g] for g in rest])
        ax.scatter(pts[:, 0], pts[:, 1], s=[34 + 20 * deg[g] for g in rest],
                   c=_FALLBACK_COLOR, zorder=3, edgecolor="white", linewidth=1.1,
                   alpha=.75)

    hi = set(highlight) & set(G)
    if label is None:
        show = set(hi)
        for c in big:
            show |= set(sorted(c, key=lambda g: -deg[g])[:label_top])
    else:
        show = set(label) & set(G)

    # place labels, then separate them in display space so none overprint.
    # The decluttering happens in pixels, so the axes must already be at its final size
    # -- a subplots_adjust() applied by the caller afterwards would rescale every
    # position that was just resolved. Force a draw so transData is current, and set
    # the figure geometry BEFORE calling this function.
    show = sorted(show)
    ax.autoscale_view()
    ax.margins(.09)                  # must happen BEFORE labels are positioned:
    ax.set_xlim(*np.array(ax.get_xlim()) + [-.03, .03])   # changing the view
    ax.set_ylim(*np.array(ax.get_ylim()) + [-.03, .03])   # afterwards would rescale
    try:                                                   # every resolved position
        ax.figure.canvas.draw()
    except Exception:
        pass
    anchors = ax.transData.transform([pos[g] for g in show]) if show else np.empty((0, 2))
    if len(show):
        anchors[:, 1] += fontsize * 1.5                      # start just above the node
        adj = _declutter(anchors, show, ax, fontsize, anchors=anchors, pull=0.0)
        inv = ax.transData.inverted()
        for g, p in zip(show, inv.transform(adj)):
            node = pos[g]
            if np.hypot(*(p - node)) > .022:                 # drifted: draw a leader
                ax.plot([node[0], p[0]], [node[1], p[1]], color="#b9bec4", lw=.6,
                        zorder=4)
            ax.text(p[0], p[1], g, fontsize=fontsize, ha="center", va="center",
                    zorder=6, color="#111" if g in hi else "#333",
                    weight="bold" if g in hi else "normal",
                    bbox=dict(boxstyle="round,pad=.14", fc="white", ec="none",
                              alpha=.82 if g in hi else .68))

    if annotate_hulls == "numbers":
        # With many modules, names on the canvas cannot avoid drifting across
        # neighbouring outlines. A numbered badge on each blob plus a keyed legend has
        # no collisions at all, and gives long auto-generated names room to be read.
        for i, c in enumerate(big):
            P = _blob([pos[g] for g in c], pad=hull_pad)
            top = [P[np.argmax(P[:, 1]), 0], P[:, 1].max()]
            ax.scatter([top[0]], [top[1]], s=235, marker="o", c="white",
                       edgecolor=cols[i], linewidth=1.6, zorder=7)
            ax.text(top[0], top[1], str(i + 1), fontsize=8.6, ha="center",
                    va="center_baseline", color=cols[i], weight="bold", zorder=8)

    elif annotate_hulls:
        # anchor each module's name above its own outline, then separate the names from
        # each other in display space -- an un-deconflicted name drifts over a
        # neighbouring module and appears to label it
        anc, txt = [], []
        for i, c in enumerate(big):
            P = _blob([pos[g] for g in c], pad=hull_pad + .03)
            anc.append([P[np.argmax(P[:, 1]), 0], P[:, 1].max()])
            txt.append(_wrap(labels[i], 26))
        A = ax.transData.transform(np.array(anc))
        A[:, 1] += 15
        A = _declutter(A, [max(t.split("\n"), key=len) for t in txt], ax, 8.2,
                       iters=400, pad=3.0)
        for i, p in enumerate(ax.transData.inverted().transform(A)):
            if np.hypot(*(p - np.array(anc[i]))) > .05:
                ax.plot([anc[i][0], p[0]], [anc[i][1], p[1]], color=cols[i], lw=.7,
                        ls=":", zorder=4, alpha=.8)
            ax.text(p[0], p[1], txt[i], fontsize=8.2, ha="center", va="bottom",
                    color=cols[i], weight="bold", zorder=7,
                    bbox=dict(boxstyle="round,pad=.2", fc="white", ec=cols[i],
                              lw=.7, alpha=.92))

    ax.set_axis_off()
    if title:
        ax.set_title(title, loc="left", fontsize=12.5, weight="bold", pad=10)

    frame = pd.DataFrame(dict(
        module=labels, functional_class=[c for c, _ in cls_col], color=cols,
        n=[len(c) for c in big],
        members=[", ".join(sorted(c)) for c in big]))

    if legend and annotate_hulls == "numbers":
        handles = [plt.Line2D([], [], marker="o", ls="", color=c, markersize=9,
                              markeredgecolor="white") for c in cols]
        texts = [f"{i+1}. {l}  ({len(big[i])})" for i, l in enumerate(labels)]
        leg = ax.legend(handles, texts, loc="upper left", bbox_to_anchor=(1.005, 1.0),
                        frameon=False, fontsize=8.6, title="STRING modules",
                        title_fontsize=9.4, handletextpad=.6, labelspacing=.62)
        leg.get_title().set_fontweight("bold")
        leg._legend_box.align = "left"
        cls_seen, h2, t2 = [], [], []
        for cls, col in cls_col:
            if cls not in cls_seen:
                cls_seen.append(cls)
                h2.append(plt.Line2D([], [], marker="s", ls="", color=col,
                                     markersize=9))
                t2.append(cls)
        leg2 = ax.legend(h2, t2, loc="lower left", bbox_to_anchor=(1.005, 0),
                         frameon=False, fontsize=8.6, title="functional class",
                         title_fontsize=9.4)
        leg2.get_title().set_fontweight("bold")
        leg2._legend_box.align = "left"
        ax.add_artist(leg)
    elif legend:
        seen, handles, texts = set(), [], []
        for lab, (cls, _), col in zip(labels, cls_col, cols):
            if cls in seen:
                continue
            seen.add(cls)
            handles.append(plt.Line2D([], [], marker="o", ls="", color=col,
                                      markersize=9, markeredgecolor="white"))
            texts.append(cls)
        ax.legend(handles, texts, loc="upper left", frameon=False, fontsize=8.5,
                  title="functional class", title_fontsize=8.5, ncol=1,
                  bbox_to_anchor=(0, 1))
    ax.text(.995, .005, f"{G.number_of_nodes()} proteins · {G.number_of_edges()} edges · "
                        f"STRING score ≥ {score/1000:.2f} · {len(big)} modules",
            transform=ax.transAxes, fontsize=8, color="#666", ha="right")
    if verbose:
        miss = sorted(set(highlight) - set(G))
        if miss:
            print(f"no edge at score >= {score}: {', '.join(miss)}")
    return ax, G, frame
