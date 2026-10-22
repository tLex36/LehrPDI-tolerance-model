"""Gene Ontology and SynGO annotation, enrichment, and a sunburst view.

    from go import enrich, gene_terms, sunburst
    syn = enrich(hits, background=all_genes, source="SynGO")
    go_bp = enrich(hits, background=all_genes, aspect="P")

SynGO is not fetched from syngoportal.org (no stable bulk URL); the SynGO consortium's
annotations are deposited in the Gene Ontology itself and are selected here by the GAF's
"assigned by" field, which gives 1042 mouse genes across 213 synaptic terms.
"""
import gzip as _gzip
import os as _os

import pandas as pd

from netutil import pkg_path, urlopen_retry

OBO_URL = "https://current.geneontology.org/ontology/go-basic.obo"
GAF_URL = "https://current.geneontology.org/annotations/{organism}.gaf.gz"
SYNAPSE_ROOT = "GO:0045202"            # synapse -- the SynGO cellular-component root
ASPECTS = {"P": "biological_process", "C": "cellular_component", "F": "molecular_function"}


def _cache(name, url, binary=False):
    d = pkg_path("go_cache")
    _os.makedirs(d, exist_ok=True)
    path = _os.path.join(d, name)
    if not _os.path.exists(path):
        data = urlopen_retry(url, timeout=300)
        mode = "wb" if binary else "w"
        if not binary and isinstance(data, bytes):
            data = data.decode()
        open(path, mode).write(data)
    return path


# ---------------------------------------------------------------- ontology
_OBO = None

def load_obo(force=False) -> dict:
    """{GO id: {'name', 'namespace', 'parents', 'obsolete'}} from go-basic.obo."""
    global _OBO
    if _OBO is not None and not force:
        return _OBO
    path = _cache("go-basic.obo", OBO_URL)
    terms, cur = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line == "[Term]":
            cur = {"parents": [], "obsolete": False, "name": "", "namespace": "",
                   "subsets": set()}
            continue
        if line.startswith("["):                       # [Typedef] etc -- stop collecting
            if cur and cur.get("id"):
                terms[cur["id"]] = cur
            cur = None
            continue
        if cur is None or not line:
            if cur and cur.get("id") and not line:
                terms[cur["id"]] = cur
                cur = None
            continue
        k, _, v = line.partition(": ")
        if k == "id" and v.startswith("GO:"):
            cur["id"] = v
        elif k == "name":
            cur["name"] = v
        elif k == "namespace":
            cur["namespace"] = v
        elif k == "is_a":
            cur["parents"].append(v.split(" ! ")[0])
        elif k == "relationship" and v.startswith(("part_of", "regulates")):
            cur["parents"].append(v.split()[1])        # part_of behaves like is_a for propagation
        elif k == "subset":
            cur["subsets"].add(v)
        elif k == "is_obsolete" and v == "true":
            cur["obsolete"] = True
    if cur and cur.get("id"):
        terms[cur["id"]] = cur
    _OBO = {k: v for k, v in terms.items() if not v["obsolete"]}
    return _OBO


def ancestors(term, obo=None, _memo=None) -> set:
    """All ancestor terms of `term`, inclusive of itself."""
    obo = obo or load_obo()
    memo = _memo if _memo is not None else {}
    if term in memo:
        return memo[term]
    out = {term}
    for p in obo.get(term, {}).get("parents", []):
        if p in obo:
            out |= ancestors(p, obo, memo)
    memo[term] = out
    return out


def children_map(obo=None) -> dict:
    obo = obo or load_obo()
    kids = {}
    for t, d in obo.items():
        for p in d["parents"]:
            kids.setdefault(p, set()).add(t)
    return kids


# ---------------------------------------------------------------- annotations
_GAF_CACHE, _TERMS_CACHE = {}, {}


def load_gaf(organism="mgi") -> pd.DataFrame:
    """gene / go_id / evidence / aspect / assigned_by from a GO annotation file."""
    if organism in _GAF_CACHE:
        return _GAF_CACHE[organism]
    path = _cache(f"{organism}.gaf.gz", GAF_URL.format(organism=organism), binary=True)
    rows = []
    with _gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 15 or f[3].strip() == "NOT":    # drop explicit negative annotations
                continue
            rows.append((f[2], f[4], f[6], f[8], f[14]))
    _GAF_CACHE[organism] = pd.DataFrame(
        rows, columns=["gene", "go_id", "evidence", "aspect", "assigned_by"])
    return _GAF_CACHE[organism]


def gene_terms(organism="mgi", source=None, aspect=None, propagate=True,
               evidence_exclude=("IEA",)) -> dict:
    """{gene symbol: set of GO ids}.

    source:    'SynGO' restricts to the SynGO consortium's annotations; None = all of GO
    aspect:    'P' biological process, 'C' cellular component, 'F' molecular function
    propagate: also annotate every ancestor term, which is required for enrichment to be
               meaningful -- without it a gene annotated to a deep term does not count
               toward its parents.
    """
    key = (organism, source, aspect, propagate, tuple(evidence_exclude or ()))
    if key in _TERMS_CACHE:
        return _TERMS_CACHE[key]
    gaf = load_gaf(organism)
    if source:
        gaf = gaf[gaf.assigned_by.str.lower() == str(source).lower()]
    if aspect:
        gaf = gaf[gaf.aspect == aspect]
    if evidence_exclude:
        gaf = gaf[~gaf.evidence.isin(evidence_exclude)]
    obo = load_obo()
    memo = {}
    out = {}
    for g, t in zip(gaf.gene, gaf.go_id):
        if t not in obo:
            continue
        out.setdefault(g, set()).update(ancestors(t, obo, memo) if propagate else {t})
    _TERMS_CACHE[key] = out
    return out


def list_slims() -> dict:
    """{slim name: number of terms} for every GO slim declared in the ontology."""
    counts = {}
    for d in load_obo().values():
        for sub in d.get("subsets", ()):
            if sub.startswith("goslim"):
                counts[sub] = counts.get(sub, 0) + 1
    return dict(sorted(counts.items(), key=lambda kv: -kv[1]))


def term_genes(min_size=3, max_size=1000, slim=None, **kw) -> dict:
    """{GO id: set of gene symbols}, filtered by term size.

    slim: keep only terms in a GO slim subset (see list_slims()), e.g.
          'goslim_generic' (140 broad terms) or 'goslim_agr' (53, very broad).
          Annotations are still propagated from the full ontology first, so a gene
          annotated to a deep term counts toward its slim ancestors."""
    inv = {}
    for g, ts in gene_terms(**kw).items():
        for t in ts:
            inv.setdefault(t, set()).add(g)
    if slim:
        obo = load_obo()
        inv = {t: g for t, g in inv.items() if slim in obo.get(t, {}).get("subsets", ())}
    return {t: g for t, g in inv.items() if min_size <= len(g) <= max_size}


# ---------------------------------------------------------------- enrichment
def enrich(genes, background=None, organism="mgi", source=None, aspect=None,
           min_size=3, max_size=1000, method="fdr_bh", slim=None,
           exclude=(), exclude_descendants=True, **kw) -> pd.DataFrame:
    """Hypergeometric over-representation against GO (or the SynGO subset of GO).

    exclude: GO ids to drop from the result -- for terms you know are uninformative.
             GO:0043209 "myelin sheath" is the usual case: its annotations come from
             myelin biochemical purifications, so it carries mitochondria, chaperones and
             cytoskeleton and swamps any membrane-proteomics query.
    exclude_descendants: also drop terms beneath the excluded ones.
    """
    from scipy.stats import hypergeom
    obo = load_obo()
    sets = term_genes(min_size=min_size, max_size=max_size, organism=organism,
                      source=source, aspect=aspect, slim=slim, **kw)
    if exclude:
        drop = set(exclude if not isinstance(exclude, str) else [exclude])
        if exclude_descendants:
            memo = {}
            drop |= {t for t in sets if ancestors(t, obo, memo) & drop}
        sets = {t: g for t, g in sets.items() if t not in drop}
    annotated = set().union(*sets.values()) if sets else set()
    bg = set(background) & annotated if background is not None else annotated
    hits = set(genes) & bg
    M, N = len(bg), len(hits)
    rows = []
    for t, members in sets.items():
        inset = members & bg
        n, k = len(inset), len(inset & hits)
        if k == 0 or n == 0:
            continue
        rows.append({"go_id": t, "name": obo.get(t, {}).get("name", t),
                     "aspect": obo.get(t, {}).get("namespace", ""),
                     "n_term": n, "n_overlap": k,
                     "expected": round(n * N / M, 2) if M else float("nan"),
                     "fold_enrichment": round((k / N) / (n / M), 2) if M and N else float("nan"),
                     "p_value": hypergeom.sf(k - 1, M, n, N),
                     "genes": ",".join(sorted(inset & hits))})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    try:
        from statsmodels.stats.multitest import multipletests
        out["q_value"] = multipletests(out.p_value, method=method)[1]
    except ImportError:
        o = out.sort_values("p_value").reset_index(drop=True)
        out = o.assign(q_value=(o.p_value * len(o) / (o.index + 1)).cummin()[::-1].cummin()[::-1].clip(upper=1.0))
    out = out.sort_values("p_value").reset_index(drop=True)
    out.attrs["n_query"], out.attrs["n_background"] = N, M
    return out


def direction_table(d, organism="mgi", source="SynGO", aspect=None, min_size=3,
                    max_size=1000, method="fdr_bh", slim=None, **kw) -> pd.DataFrame:
    """Per-term fold-change direction: are a term's genes shifted up or down as a group?

    d: a prepared frame (gene / log2fc / pval). Mann-Whitney vs all other detected genes,
    the same test pathway_activity() uses for KEGG.
    """
    from scipy.stats import mannwhitneyu
    obo = load_obo()
    vals = d.set_index("gene").log2fc.groupby(level=0).mean()
    sets = term_genes(min_size=min_size, max_size=max_size, organism=organism,
                      source=source, aspect=aspect, slim=slim, **kw)
    rows = []
    for t, members in sets.items():
        inside_g = [g for g in members if g in vals.index]
        if len(inside_g) < min_size:
            continue
        inside = vals.loc[inside_g]
        outside = vals.drop(index=inside_g, errors="ignore")
        if len(outside) < min_size:
            continue
        try:
            _, p = mannwhitneyu(inside, outside, alternative="two-sided")
        except ValueError:
            continue
        rows.append({"go_id": t, "name": obo.get(t, {}).get("name", t),
                     "aspect": obo.get(t, {}).get("namespace", ""),
                     "n_genes": len(inside_g),
                     "median_log2fc": float(inside.median()),
                     "median_shift": float(inside.median() - outside.median()),
                     "n_up": int((inside > 0).sum()), "n_down": int((inside < 0).sum()),
                     "p_value": float(p),
                     "genes": ",".join(sorted(inside_g))})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    try:
        from statsmodels.stats.multitest import multipletests
        out["q_value"] = multipletests(out.p_value, method=method)[1]
    except ImportError:
        o = out.sort_values("p_value").reset_index(drop=True)
        out = o.assign(q_value=(o.p_value * len(o) / (o.index + 1)).cummin()[::-1].cummin()[::-1].clip(upper=1.0))
    return out.sort_values("median_shift", ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------- sunburst
def _tree(root, sets, obo, max_depth=4):
    """Induced hierarchy over `sets`: each term's parent is its most specific ancestor that
    is also in `sets`. Terms need not be contiguous in the DAG, so this works on filtered
    sets and GO slims where intermediate terms are absent. GO is a DAG, so a term with
    several eligible ancestors is attached to the smallest (most specific) one.
    """
    memo = {}
    anc = {t: ancestors(t, obo, memo) for t in sets}
    under = {t for t in sets if root in anc[t]}
    tree = {t: [] for t in under}
    for t in under:
        if t == root:
            continue
        cands = [a for a in anc[t] if a != t and a in under]
        if not cands:
            continue
        parent = min(cands, key=lambda a: (len(sets[a]), a))
        tree[parent].append(t)

    depth, keep = {root: 0}, {root}            # trim anything past max_depth
    frontier = [root]
    while frontier:
        nxt = []
        for n in frontier:
            if depth[n] >= max_depth:
                continue
            for c in tree.get(n, ()):
                if c not in keep:
                    depth[c] = depth[n] + 1
                    keep.add(c)
                    nxt.append(c)
        frontier = nxt
    return {n: [c for c in tree.get(n, ()) if c in keep] for n in keep}


def sunburst(d, source="SynGO", organism="mgi", root=SYNAPSE_ROOT, max_depth=4,
             min_size=3, max_size=1000, cmap="RdBu_r", vlim=None, center=None,
             ax=None, figsize=(9, 9), label_min_angle=11, label_size=7.5,
             title=None, return_data=False, **kw):
    """Sunburst of the SynGO/GO hierarchy, wedges coloured by median log2 fold change.

    Wedge angle is proportional to the number of detected genes annotated to that term,
    rings are ontology depth below `root`, and colour is the term's median log2fc relative
    to `center` (default: the median across all detected genes, which is the right
    reference for an IP where everything is enriched over control).
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import cm, colors as mcolors
    from matplotlib.patches import Wedge

    obo = load_obo()
    sets = term_genes(min_size=min_size, max_size=max_size, organism=organism,
                      source=source, slim=kw.pop("slim", None), **kw)
    vals = d.set_index("gene").log2fc.groupby(level=0).mean()
    sets = {t: (g & set(vals.index)) for t, g in sets.items()}
    sets = {t: g for t, g in sets.items() if len(g) >= min_size}
    if root not in sets:
        raise ValueError(f"root {root} has fewer than {min_size} detected genes")

    tree = _tree(root, sets, obo, max_depth=max_depth)
    base = float(vals.median()) if center is None else float(center)
    lim = vlim or max(abs(float(vals.loc[list(sets[root])].median() - base)), 0.3)
    for t, g in sets.items():
        lim = max(lim, abs(float(vals.loc[list(g)].median()) - base))
    norm = mcolors.TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=lim)
    cmo = plt.get_cmap(cmap)

    rows, wedges = [], []

    def place(node, depth, t0, t1):
        shift = float(vals.loc[list(sets[node])].median()) - base
        rows.append({"go_id": node, "name": obo[node]["name"], "depth": depth,
                     "n_genes": len(sets[node]), "median_log2fc": float(vals.loc[list(sets[node])].median()),
                     "median_shift": shift})
        wedges.append((depth, t0, t1, shift, obo[node]["name"], len(sets[node])))
        kids = [c for c in tree.get(node, []) if c in sets]
        if not kids or depth + 1 > max_depth:
            return
        total = sum(len(sets[c]) for c in kids)
        if not total:
            return
        a = t0
        for c in kids:                                    # children share the parent's arc
            w = (t1 - t0) * len(sets[c]) / total
            place(c, depth + 1, a, a + w)
            a += w

    place(root, 0, 0.0, 360.0)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={"aspect": "equal"})
    else:
        fig = ax.figure
    ring = 1.0
    max_depth_drawn = max(w[0] for w in wedges)
    for depth, t0, t1, shift, name, n in wedges:
        r0, r1 = depth * ring, (depth + 1) * ring
        ax.add_patch(Wedge((0, 0), r1, t0, t1, width=ring, facecolor=cmo(norm(shift)),
                           edgecolor="white", linewidth=0.8))
        span = t1 - t0
        if depth == 0:
            ax.text(0, 0, f"{name}\n({n})", ha="center", va="center",
                    fontsize=label_size + 1.5, fontweight="bold")
            continue
        # a wedge needs enough arc to hold text; outer rings are thinner in angle terms
        if span < label_min_angle:
            continue
        mid_deg = (t0 + t1) / 2
        mid = np.deg2rad(mid_deg)
        if depth == max_depth_drawn and span < label_min_angle * 2.2:
            # outermost thin wedges: label outside the ring with a leader line
            rr = r1 + 0.12
            ha = "left" if -90 <= (mid_deg % 360) <= 90 or (mid_deg % 360) >= 270 else "right"
            rot = mid_deg if ha == "left" else mid_deg - 180
            ax.text(np.cos(mid) * rr, np.sin(mid) * rr, f"{name} ({n})",
                    ha=ha, va="center", fontsize=label_size - 0.5, rotation=rot,
                    rotation_mode="anchor", color="#333")
            continue
        # inside the wedge: truncate to what the arc can hold
        budget = max(6, int(span / 3.2))
        label = name if len(name) <= budget else name[:budget - 1] + "\u2026"
        rot = mid_deg - 180 if 90 < (mid_deg % 360) <= 270 else mid_deg
        ax.text(np.cos(mid) * (r0 + r1) / 2, np.sin(mid) * (r0 + r1) / 2,
                f"{label}\n({n})", ha="center", va="center", fontsize=label_size,
                rotation=rot, rotation_mode="anchor")
    lim_r = (max_depth_drawn + 1) * ring * 1.38
    ax.set_xlim(-lim_r, lim_r); ax.set_ylim(-lim_r, lim_r)
    ax.axis("off")
    sm = cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label(f"median log$_2$FC vs all genes ({base:.2f})", fontsize=9)
    if title:
        ax.set_title(title, fontsize=12)
    out = pd.DataFrame(rows).sort_values(["depth", "n_genes"], ascending=[True, False])
    return (ax, out) if return_data else ax



def _prune(tree, root, passes):
    """Drop subtrees containing no passing term. A node survives if it passes itself or
    any descendant does, so an interesting deep term keeps its path back to the root."""
    kept = {}

    def visit(node):
        kids = [c for c in tree.get(node, []) if visit(c)] if tree.get(node) else []
        alive = passes(node) or bool(kids)
        if alive:
            kept[node] = kids
        return alive

    visit(root)
    kept.setdefault(root, [c for c in tree.get(root, []) if c in kept])
    return kept


# ---------------------------------------------------------------- hierarchy tree
def tree_plot(d=None, source="SynGO", organism="mgi", root=SYNAPSE_ROOT, max_depth=4,
              min_size=3, max_size=1000, aspect=None, slim=None, cmap="RdBu_r",
              vlim=None, center=None, ax=None, figsize=None, label_size=8,
              node_scale=9, branch_color="#999999", branch_lw=0.9, show_counts=True,
              prune_min_genes=None, prune_shift=None, prune_q=None, stats=None,
              return_data=False, **kw):
    """Rectangular cladogram of the GO hierarchy, in the style of a phylogenetic tree.

    One leaf per terminal term, depth on the x axis, elbow connectors between parent and
    child. With `d` (a prepared frame) each node is a dot sized by the number of detected
    genes and coloured by that term's median log2 fold change; without `d` the tree is
    drawn structurally, sized by annotated gene count.

    slim: restrict to a GO slim subset of broad terms (see list_slims()), e.g.
          'goslim_generic'. Use with root='GO:0008150' (biological process),
          'GO:0005575' (cellular component) or 'GO:0003674' (molecular function).
    prune_min_genes / prune_shift / prune_q: drop branches whose terms all fail these
          thresholds. A node survives if it passes or any descendant does, so a deep
          interesting term keeps its path back to the root. prune_q needs `stats`
          (a direction_table frame) for the q values.

    GO is a DAG rather than a tree -- a term can have several parents -- so each term is
    attached to the first parent that reaches it in a breadth-first walk from `root`,
    preferring the larger parent. The printed structure is therefore one valid rendering
    of the hierarchy, not the only one.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import cm, colors as mcolors

    obo = load_obo()
    sets = term_genes(min_size=min_size, max_size=max_size, organism=organism,
                      source=source, aspect=aspect, slim=slim, **kw)
    if root not in sets:        # aspect roots are in no slim, and exceed any max_size
        full = term_genes(min_size=1, max_size=10**9, organism=organism,
                          source=source, aspect=aspect, **kw)
        if root not in full:
            raise ValueError(f"root {root} has no annotations in this dataset")
        sets[root] = full[root]
    vals = None
    if d is not None:
        vals = d.set_index("gene").log2fc.groupby(level=0).mean()
        sets = {t: (g & set(vals.index)) for t, g in sets.items()}
        sets = {t: g for t, g in sets.items() if len(g) >= min_size}
    if root not in sets:
        raise ValueError(f"root {root} not present with >= {min_size} genes")

    tree = _tree(root, sets, obo, max_depth=max_depth)

    if prune_min_genes or prune_shift or prune_q:
        qmap = dict(zip(stats.go_id, stats.q_value)) if stats is not None else {}
        med = float(vals.median()) if vals is not None else 0.0

        def passes(t):
            if t == root:
                return True
            if prune_min_genes and len(sets.get(t, ())) < prune_min_genes:
                return False
            if prune_shift is not None and vals is not None:
                if abs(float(vals.loc[list(sets[t])].median()) - med) < prune_shift:
                    return False
            if prune_q is not None and qmap.get(t, 1.0) > prune_q:
                return False
            return True

        tree = _prune(tree, root, passes)
        sets = {t: g for t, g in sets.items() if t in tree}

    # ---- layout: leaves get consecutive rows, parents sit at the mean of their children
    ypos, order = {}, []
    counter = [0]

    def walk(node, depth):
        kids = sorted([c for c in tree.get(node, []) if c in sets],
                      key=lambda t: -len(sets[t]))
        order.append((node, depth))
        if not kids or depth >= max_depth:
            ypos[node] = counter[0]
            counter[0] += 1
            return ypos[node]
        ys = [walk(c, depth + 1) for c in kids]
        ypos[node] = (min(ys) + max(ys)) / 2
        return ypos[node]

    walk(root, 0)
    n_leaf = counter[0]
    depth_max = max(dp for _, dp in order)

    base = (float(vals.median()) if center is None else float(center)) if vals is not None else 0.0
    shifts = {}
    for t, g in sets.items():
        shifts[t] = (float(vals.loc[list(g)].median()) - base) if vals is not None else 0.0
    lim = vlim or (max((abs(shifts[t]) for t, _ in order), default=1.0) or 1.0)
    norm = mcolors.TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=lim)
    cmo = plt.get_cmap(cmap)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize or (12, max(4, 0.34 * n_leaf + 1)))
    else:
        fig = ax.figure

    for node, depth in order:                             # elbow connectors
        for c in [c for c in tree.get(node, []) if c in sets and c in ypos]:
            ax.plot([depth, depth, depth + 1], [ypos[node], ypos[c], ypos[c]],
                    color=branch_color, lw=branch_lw, zorder=1, solid_joinstyle="miter")

    for node, depth in order:
        n = len(sets[node])
        ax.scatter([depth], [ypos[node]], s=node_scale * np.sqrt(n) * 2.2,
                   c=[cmo(norm(shifts[node]))] if vals is not None else ["#cccccc"],
                   edgecolors="black", linewidths=0.5, zorder=3)
        name = obo[node]["name"]
        txt = f"{name} ({n})" if show_counts else name
        has_kids = bool([c for c in tree.get(node, []) if c in sets and c in ypos])
        if has_kids:            # internal node: label above the branch, clear of the children
            ax.text(depth + 0.04, ypos[node] - 0.22, txt, va="bottom", ha="left",
                    fontsize=label_size - 0.5, zorder=4, color="#333",
                    fontstyle="italic")
        else:                   # leaf: label to the right, as in a phylogram
            ax.text(depth + 0.07, ypos[node], " " + txt, va="center", ha="left",
                    fontsize=label_size, zorder=4)

    ax.set_xlim(-0.4, depth_max + 2.4)
    ax.set_ylim(-0.8, n_leaf - 0.2)
    ax.set_xticks(range(depth_max + 1))
    ax.set_xticklabels([f"depth {i}" for i in range(depth_max + 1)], fontsize=8)
    ax.set_yticks([])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.invert_yaxis()
    if vals is not None:
        sm = cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
        cb = fig.colorbar(sm, ax=ax, shrink=0.45, pad=0.01)
        cb.set_label(f"median log$_2$FC vs all genes ({base:.2f})", fontsize=9)
    fig.tight_layout()

    out = pd.DataFrame([{"go_id": t, "name": obo[t]["name"], "depth": dp,
                         "n_genes": len(sets[t]), "median_shift": shifts[t]}
                        for t, dp in order])
    return (ax, out) if return_data else ax


# ---------------------------------------------------------------- interactive HTML
def _node_stats(sets, vals, base, stats=None, top_genes=18):
    """Per-term hover payload: counts, medians, q, and the member genes with their values."""
    qmap = dict(zip(stats.go_id, stats.q_value)) if stats is not None else {}
    pmap = dict(zip(stats.go_id, stats.p_value)) if stats is not None else {}
    out = {}
    for t, genes in sets.items():
        g = sorted(genes)
        if vals is not None:
            v = vals.loc[g].sort_values(ascending=False)
            med = float(v.median())
            shown = list(v.items())
            if len(shown) > top_genes:
                half = top_genes // 2
                listing = ("<br>".join(f"&nbsp;&nbsp;{k} {x:+.2f}" for k, x in shown[:half])
                           + f"<br>&nbsp;&nbsp;… {len(shown) - top_genes} more …<br>"
                           + "<br>".join(f"&nbsp;&nbsp;{k} {x:+.2f}" for k, x in shown[-half:]))
            else:
                listing = "<br>".join(f"&nbsp;&nbsp;{k} {x:+.2f}" for k, x in shown)
            out[t] = {"n": len(g), "median": med, "shift": med - base,
                      "n_up": int((v > 0).sum()), "n_down": int((v < 0).sum()),
                      "genes": listing, "q": qmap.get(t), "p": pmap.get(t)}
        else:
            out[t] = {"n": len(g), "median": None, "shift": 0.0, "n_up": None, "n_down": None,
                      "genes": "<br>".join(f"&nbsp;&nbsp;{x}" for x in g[:top_genes]),
                      "q": qmap.get(t), "p": pmap.get(t)}
    return out


def _hover(name, go_id, st):
    rows = [f"<b>{name}</b>", f"{go_id} &nbsp;·&nbsp; {st['n']} genes"]
    if st["median"] is not None:
        rows.append(f"median log2FC {st['median']:+.2f} &nbsp; (shift {st['shift']:+.2f})")
        rows.append(f"{st['n_up']} up / {st['n_down']} down")
    if st.get("q") is not None:
        rows.append(f"p = {st['p']:.2g} &nbsp; q = {st['q']:.2g}")
    rows.append("<br>" + st["genes"])
    return "<br>".join(rows)


def tree_html(d=None, source="SynGO", organism="mgi", root=SYNAPSE_ROOT, max_depth=4,
              min_size=3, max_size=1000, aspect=None, slim=None, center=None,
              colorscale="RdBu_r", stats=None, prune_min_genes=None, prune_shift=None,
              prune_q=None, out=None, width=1100, height=None, title=None,
              top_genes=18, **kw):
    """Interactive cladogram: hover a node for its genes and their fold changes.

    Same arguments as tree_plot. `out` writes a standalone .html file; otherwise the
    plotly Figure is returned, which renders inline in a notebook via .show().
    """
    import numpy as np
    import plotly.graph_objects as pgo

    obo = load_obo()
    sets = term_genes(min_size=min_size, max_size=max_size, organism=organism,
                      source=source, aspect=aspect, slim=slim, **kw)
    if root not in sets:
        full = term_genes(min_size=1, max_size=10**9, organism=organism,
                          source=source, aspect=aspect, **kw)
        if root not in full:
            raise ValueError(f"root {root} has no annotations in this dataset")
        sets[root] = full[root]

    vals = None
    if d is not None:
        vals = d.set_index("gene").log2fc.groupby(level=0).mean()
        sets = {t: (g & set(vals.index)) for t, g in sets.items()}
        sets = {t: g for t, g in sets.items() if len(g) >= min_size or t == root}

    tree = _tree(root, sets, obo, max_depth=max_depth)
    base = (float(vals.median()) if center is None else float(center)) if vals is not None else 0.0

    if prune_min_genes or prune_shift or prune_q:
        qmap = dict(zip(stats.go_id, stats.q_value)) if stats is not None else {}

        def passes(t):
            if t == root:
                return True
            if prune_min_genes and len(sets.get(t, ())) < prune_min_genes:
                return False
            if prune_shift is not None and vals is not None:
                if abs(float(vals.loc[list(sets[t])].median()) - base) < prune_shift:
                    return False
            if prune_q is not None and qmap.get(t, 1.0) > prune_q:
                return False
            return True

        tree = _prune(tree, root, passes)
        sets = {t: g for t, g in sets.items() if t in tree}

    ypos, order, counter = {}, [], [0]

    def walk(node, depth):
        kids = sorted([c for c in tree.get(node, []) if c in sets], key=lambda t: -len(sets[t]))
        order.append((node, depth))
        if not kids or depth >= max_depth:
            ypos[node] = counter[0]; counter[0] += 1
            return ypos[node]
        ys = [walk(c, depth + 1) for c in kids]
        ypos[node] = (min(ys) + max(ys)) / 2
        return ypos[node]

    walk(root, 0)
    st = _node_stats(sets, vals, base, stats, top_genes)

    bx, by = [], []                                       # elbow connectors in one trace
    for node, depth in order:
        for c in [c for c in tree.get(node, []) if c in ypos]:
            bx += [depth, depth, depth + 1, None]
            by += [ypos[node], ypos[c], ypos[c], None]

    lim = max((abs(st[t]["shift"]) for t, _ in order), default=1.0) or 1.0
    fig = pgo.Figure()
    fig.add_trace(pgo.Scatter(x=bx, y=by, mode="lines", line=dict(color="#aaaaaa", width=1),
                              hoverinfo="skip", showlegend=False))
    fig.add_trace(pgo.Scatter(
        x=[dp for _, dp in order], y=[ypos[t] for t, _ in order],
        mode="markers+text",
        text=[f"  {obo[t]['name']} ({st[t]['n']})" for t, _ in order],
        textposition="middle right", textfont=dict(size=11),
        marker=dict(size=[8 + 2.2 * np.sqrt(st[t]["n"]) for t, _ in order],
                    color=[st[t]["shift"] for t, _ in order],
                    colorscale=colorscale, cmid=0, cmin=-lim, cmax=lim,
                    line=dict(color="black", width=0.6),
                    colorbar=dict(title=dict(text="median<br>log2FC<br>shift", side="right"),
                                  thickness=14, len=0.5)),
        customdata=[_hover(obo[t]["name"], t, st[t]) for t, _ in order],
        hovertemplate="%{customdata}<extra></extra>", showlegend=False))

    n_leaf, depth_max = counter[0], max(dp for _, dp in order)
    fig.update_layout(
        title=title or f"{source or 'GO'} hierarchy",
        width=width, height=height or max(420, 26 * n_leaf + 120),
        xaxis=dict(title="ontology depth", range=[-0.5, depth_max + 2.2],
                   tickmode="array", tickvals=list(range(depth_max + 1)),
                   showgrid=False, zeroline=False),
        yaxis=dict(visible=False, autorange="reversed"),
        plot_bgcolor="white", hoverlabel=dict(bgcolor="white", font_size=11,
                                             bordercolor="#888", align="left"),
        margin=dict(l=20, r=20, t=60, b=50))
    if out:
        fig.write_html(out, include_plotlyjs="cdn")
        print(f"wrote {out}")
    return fig


def sunburst_html(d, source="SynGO", organism="mgi", root=SYNAPSE_ROOT, max_depth=4,
                  min_size=3, max_size=1000, aspect=None, slim=None, center=None,
                  colorscale="RdBu_r", stats=None, out=None, width=820, height=820,
                  title=None, top_genes=18, **kw):
    """Interactive sunburst: click a wedge to zoom in, hover for the gene list."""
    import plotly.graph_objects as pgo

    obo = load_obo()
    sets = term_genes(min_size=min_size, max_size=max_size, organism=organism,
                      source=source, aspect=aspect, slim=slim, **kw)
    vals = d.set_index("gene").log2fc.groupby(level=0).mean()
    sets = {t: (g & set(vals.index)) for t, g in sets.items()}
    sets = {t: g for t, g in sets.items() if len(g) >= min_size}
    if root not in sets:
        raise ValueError(f"root {root} has fewer than {min_size} detected genes")

    tree = _tree(root, sets, obo, max_depth=max_depth)
    base = float(vals.median()) if center is None else float(center)
    st = _node_stats(sets, vals, base, stats, top_genes)

    parent_of = {c: p for p, cs in tree.items() for c in cs}
    ids, labels, parents, values, colors, hover = [], [], [], [], [], []
    for node, kids in tree.items():
        # sibling GO terms share genes, so a parent's gene count can be smaller than the
        # sum of its children's. Size each wedge by the genes NOT in any drawn child and
        # let plotly add the children on top ("remainder"), which can never overflow.
        own = sets[node] - set().union(*(sets[k] for k in kids)) if kids else sets[node]
        ids.append(node); labels.append(obo[node]["name"])
        parents.append(parent_of.get(node, ""))
        values.append(max(len(own), 1 if not kids else 0))
        colors.append(st[node]["shift"])
        hover.append(_hover(obo[node]["name"], node, st[node]))

    lim = max((abs(c) for c in colors), default=1.0) or 1.0
    fig = pgo.Figure(pgo.Sunburst(
        ids=ids, labels=labels, parents=parents, values=values, branchvalues="remainder",
        marker=dict(colors=colors, colorscale=colorscale, cmid=0, cmin=-lim, cmax=lim,
                    line=dict(color="white", width=1),
                    colorbar=dict(title=dict(text="median<br>log2FC<br>shift", side="right"),
                                  thickness=14, len=0.5)),
        customdata=hover, hovertemplate="%{customdata}<extra></extra>",
        insidetextorientation="radial"))
    fig.update_layout(title=title or f"{source or 'GO'} sunburst", width=width, height=height,
                      margin=dict(l=10, r=10, t=60, b=10),
                      hoverlabel=dict(bgcolor="white", font_size=11, bordercolor="#888",
                                      align="left"))
    if out:
        fig.write_html(out, include_plotlyjs="cdn")
        print(f"wrote {out}")
    return fig
