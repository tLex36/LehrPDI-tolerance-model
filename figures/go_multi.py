"""Interactive GO explorer across several datasets, with selectable contrasts.

    from go_multi import multi_explorer
    multi_explorer({
        "development": dict(raw=dev_proteins, denom="Ctrl",
                            order=["P6","P12","P18","P21","P35"]),
        "region":      dict(raw=reg_proteins, denom=None),
    }, source="SynGO")

Pick a dataset, then one sample (enrichment over its own control) or two (the difference
between them). Node colour, the gene table and the volcano all follow the selection.
"""
import html as _html
import json as _json
import math as _math

import numpy as np
import pandas as pd

import re as _re

import go as _go
import volcano as _vol

SYNGO_VIEWS = {"cellular localization": dict(aspect="C", root=_go.SYNAPSE_ROOT)}
FULL_VIEWS = {"cellular localization": dict(aspect="C", root="GO:0005575"),
              "biological process":    dict(aspect="P", root="GO:0008150")}

# annotation sources offered by the first dropdown. GO is far larger than SynGO, so it
# gets a higher min_size and a max_size that drops the uninformative top-level terms.
SOURCES = {
    "SynGO": dict(source="SynGO", views=SYNGO_VIEWS, min_size=3, max_size=1000),
    "GO (all)": dict(source=None, views=FULL_VIEWS, min_size=10, max_size=400),
}


def _collapse(d):
    """One row per gene symbol -- these sheets carry duplicates (isoforms, shared peptides)."""
    g = d.groupby("gene", as_index=True)
    return pd.DataFrame({"log2fc": g.log2fc.mean(), "pval": g.pval.min()})


def _contrast_frame(series, a, b=None):
    """One sample vs its control, or the difference between two samples."""
    A = _collapse(series[a])
    if b is None:
        return pd.DataFrame({"gene": A.index, "log2fc": A.log2fc.values,
                             "pval": A.pval.values}).reset_index(drop=True)
    B = _collapse(series[b])
    common = A.index.intersection(B.index)
    # both are ratios against a control, so the difference is the between-sample contrast.
    # there is no p for that difference; the weaker of the two inputs is carried instead.
    return pd.DataFrame({"gene": common,
                         "log2fc": (A.log2fc[common] - B.log2fc[common]).values,
                         "pval": np.maximum(A.pval[common], B.pval[common])}).reset_index(drop=True)


_UNSET = object()



def _p_col_for(df, ratio_col):
    """The p-value column paired with a given 'Abundance Ratio: (A) / (B)' column."""
    m = _re.match(r"abundance ratio:\s*\((?P<num>[^)]+)\)\s*/\s*\((?P<den>[^)]+)\)",
                  str(ratio_col).strip(), _re.I)
    if not m:
        return None
    num, den = m.group("num").strip(), m.group("den").strip()
    return next((c for c in df.columns if "p-value" in str(c).lower()
                 and f"({num})" in str(c) and f"({den})" in str(c)), None)


def _explicit_frames(spec, default_raw=None, **prep):
    """Build {label: frame} from a `contrasts` mapping.

    Each entry is dict(raw=<sheet>, col=<ratio column>, invert=False) -- `raw` may be
    omitted to use the dataset's own sheet. `invert` flips the ratio, so a stored
    '(DelESEV3) / (2D150)' column can be presented as '2D150 / DelESEV3'.
    """
    out = {}
    for label, cfg in spec.items():
        if isinstance(cfg, str):
            cfg = {"col": cfg}
        raw = cfg.get("raw", default_raw)
        if raw is None:
            raise ValueError(f"contrast {label!r}: no sheet given (pass raw= here or on the dataset)")
        col = cfg["col"]
        if col not in raw.columns:
            raise ValueError(f"contrast {label!r}: column {col!r} not in sheet; "
                             f"ratio columns are "
                             f"{[c for c in raw.columns if 'abundance ratio' in str(c).lower() and 'p-value' not in str(c).lower()]}")
        pc = cfg.get("p_col") or _p_col_for(raw, col)
        kw = {**prep, **{k: v for k, v in cfg.items()
                         if k in ("require_p", "missing_p", "contaminants")}}
        d = _vol.prepare(raw, ratio_col=col, p_col=pc, **kw)
        if cfg.get("invert"):
            d = d.assign(log2fc=-d.log2fc)
        out[label] = d
    return out


def build_multi(datasets, sources=None, source=_UNSET, views=None,
                min_size=None, max_size=None, organism="mgi", max_depth=12,
                pairs=True, verbose=True, **kw) -> dict:
    """Assemble the payload: gene matrix, a hierarchy per dataset/source/view, and node
    statistics for every contrast (each sample vs its control, and every pair).

    sources: {label: dict(source=..., views=..., min_size=..., max_size=...)} -- one entry
             per choice in the annotation dropdown. Defaults to SOURCES (SynGO + all GO).
    source/views/min_size/max_size: shorthand for a single annotation source, e.g.
             source="SynGO" builds only that one.
    """
    if sources is None:
        if source is not _UNSET:                     # single-source shorthand
            label = source or "GO (all)"
            base_cfg = dict(SOURCES.get(label, {"source": source}))
            base_cfg["source"] = source
            if views is not None:     base_cfg["views"] = views
            if min_size is not None:  base_cfg["min_size"] = min_size
            if max_size is not None:  base_cfg["max_size"] = max_size
            base_cfg.setdefault("views", SYNGO_VIEWS if source else FULL_VIEWS)
            base_cfg.setdefault("min_size", 3 if source else 10)
            base_cfg.setdefault("max_size", 1000 if source else 400)
            sources = {label: base_cfg}
        else:
            sources = SOURCES
    for bad in ("source", "views", "min_size", "max_size", "aspect", "root", "slim"):
        kw.pop(bad, None)                            # never let these reach term_genes
    obo = _go.load_obo()
    series_by, explicit = {}, set()
    for label, cfg in datasets.items():
        # resolved at call time, so reloading volcano alone is enough
        extra = {k: v for k, v in cfg.items()
                 if k in ("strip", "require_p", "missing_p", "contaminants")}
        if cfg.get("contrasts"):                     # explicitly named comparisons
            s = _explicit_frames(cfg["contrasts"], cfg.get("raw"),
                                 **{k: v for k, v in extra.items() if k != "strip"})
            series_by[label] = s
            explicit.add(label)
            continue
        s = _vol.ratio_series(cfg["raw"], denom=cfg.get("denom", "Ctrl"),
                              order=cfg.get("order"), **extra)
        if not s:
            cols = [str(c) for c in cfg["raw"].columns
                    if "abundance ratio" in str(c).lower() and "p-value" not in str(c).lower()]
            raise ValueError(
                f"dataset {label!r}: no ratio columns matched denom={cfg.get('denom', 'Ctrl')!r}.\n"
                f"  ratio columns present: {cols}\n"
                f"  pass denom=None to accept any denominator (each sample has its own control), "
                f"or denom='<exact control name>'.")
        series_by[label] = s

    genes = sorted(set().union(*(set(d.gene) for s in series_by.values() for d in s.values())))
    gidx = {g: i for i, g in enumerate(genes)}
    out = {"genes": genes, "datasets": {}}

    for label, series in series_by.items():
        samples = list(series)
        fc, pv = {}, {}
        for smp, d in series.items():
            v = d.set_index("gene")
            col_fc = [None] * len(genes)
            col_p = [None] * len(genes)
            for g, f, p in zip(v.index, v.log2fc, v.pval):
                col_fc[gidx[g]] = round(float(f), 3)
                col_p[gidx[g]] = (None if pd.isna(p) else float(p))
            fc[smp], pv[smp] = col_fc, col_p

        detected = set().union(*(set(d.gene) for d in series.values()))
        combos = [(s_, None) for s_ in samples]
        if pairs and label not in explicit:
            combos += [(samples[i], samples[j])
                       for i in range(len(samples)) for j in range(i + 1, len(samples))]
        frames = {(a if b is None else f"{a}|{b}"): _contrast_frame(series, a, b)
                  for a, b in combos}
        base = {k: round(float(f.log2fc.median()), 3) for k, f in frames.items()}

        src_out = {}
        for slabel, scfg in sources.items():
            views = scfg.get("views") or SYNGO_VIEWS
            mn, mx = scfg.get("min_size", 3), scfg.get("max_size", 1000)
            vmodels = {}
            for vlabel, over in views.items():
                cfg = dict(source=scfg.get("source"), aspect=None, root=None, slim=None)
                cfg.update(over)
                sets = _go.term_genes(min_size=mn, max_size=mx, organism=organism,
                                      source=cfg["source"], aspect=cfg["aspect"],
                                      slim=cfg["slim"], **kw)
                root = cfg["root"]
                if root not in sets:
                    full = _go.term_genes(min_size=1, max_size=10**9, organism=organism,
                                          source=cfg["source"], aspect=cfg["aspect"], **kw)
                    if root not in full:
                        continue
                    sets[root] = full[root]
                sets = {t: (g & detected) for t, g in sets.items()}
                sets = {t: g for t, g in sets.items() if len(g) >= mn or t == root}
                if root not in sets:
                    continue
                tree = _go._tree(root, sets, obo, max_depth=max_depth)
                vmodels[vlabel] = {"root": root,
                                   "nodes": {t: {"name": obo[t]["name"],
                                                 "kids": [k for k in kids if k in sets],
                                                 "g": sorted(gidx[g] for g in sets[t])}
                                             for t, kids in tree.items()}}
            stats = {}
            for k, d in frames.items():
                stats[k] = {}
                for vlabel in vmodels:
                    cfg = dict(views[vlabel])
                    t = _go.direction_table(d, source=scfg.get("source"), organism=organism,
                                            aspect=cfg.get("aspect"), slim=cfg.get("slim"),
                                            min_size=mn, max_size=mx, **kw)
                    stats[k][vlabel] = {r.go_id: [float(r.p_value), float(r.q_value)]
                                        for r in t.itertuples()
                                        if r.go_id in vmodels[vlabel]["nodes"]}
            src_out[slabel] = {"views": vmodels, "stats": stats}
            if verbose:
                print(f"  {label} / {slabel}: "
                      + ", ".join(f"{v} {len(m['nodes'])} nodes" for v, m in vmodels.items()))

        out["datasets"][label] = {"samples": samples, "fc": fc, "p": pv,
                                  "base": base, "sources": src_out,
                                  "pairs": label not in explicit}
    return out


_TPL = r"""<!doctype html>
<meta charset="utf-8">
<title>__TITLE__</title>
<style>
 body{font:13px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
      margin:0;padding:18px 22px;color:#222;background:#fff;
      padding-right:calc(__PANELW__px + 22px)}
 h1{font-size:16px;font-weight:600;margin:0 0 4px}
 .sub{color:#777;font-size:12px;margin-bottom:10px}
 .ctl{font-size:12px;color:#555;margin:6px 0 0}
 .ctl select{font-size:12px} .ctl b{color:#333}
 #crumbs{margin:10px 0 6px;font-size:12px}
 #crumbs span.c{color:#1a6fb5;cursor:pointer} #crumbs span.c:hover{text-decoration:underline}
 #crumbs span.sep{color:#bbb;margin:0 6px} #crumbs span.cur{font-weight:600}
 .link{fill:none;stroke:#bbb;stroke-width:1.1}
 .node circle{stroke:#333;stroke-width:.8;cursor:pointer}
 .node circle.leaf{stroke-dasharray:2 2;cursor:default}
 .node text{font-size:11.5px;cursor:pointer;dominant-baseline:middle}
 .node text.leaf{cursor:default;fill:#555}
 .node:hover circle{stroke-width:2}
 #grip{position:fixed;top:0;right:__PANELW__px;width:7px;height:100vh;z-index:60;
       cursor:col-resize}
 #grip:hover,#grip.drag{background:#cfe3f5}
 body.resizing{user-select:none;cursor:col-resize}
 #tip{position:fixed;top:0;right:0;width:__PANELW__px;height:100vh;z-index:50;overflow:auto;
      background:#fcfcfc;border-left:1px solid #ddd;padding:16px 16px 28px;box-sizing:border-box;
      font-size:11.5px}
 #tip h3{margin:0 0 2px;font-size:13px} #tip .meta{color:#666;margin-bottom:9px;line-height:1.5}
 #tip .hint{color:#aaa;font-style:italic;text-align:center;padding:30px 0}
 #tip table{border-collapse:collapse;width:100%}
 #tip th{text-align:left;font-weight:600;border-bottom:1px solid #ddd;padding:2px 8px 3px 0;
         position:sticky;top:0;background:#fcfcfc}
 #tip td{padding:1px 8px 1px 0;font-variant-numeric:tabular-nums}
 #tip td.v{text-align:right}
 .up{color:#b2182b}.dn{color:#2166ac}
 #vol{display:block;margin:0 0 10px}
 #vol circle.bg,#vol circle.hl{cursor:pointer}
 #vol circle.near{stroke:#111;stroke-width:1.4}
 #gtip{position:absolute;display:none;z-index:70;pointer-events:none;
       background:#222;color:#fff;border-radius:4px;padding:3px 7px;
       font-size:11px;white-space:nowrap;box-shadow:0 2px 8px rgba(0,0,0,.3)}
 #vol .bg{fill:#d6d6d6;fill-opacity:.55} #vol .hl{stroke:#000;stroke-width:.5}
 #vol .ax{stroke:#bbb;stroke-width:.7}
 #vol text.lab{font-size:8.5px;fill:#222} #vol text.axl{font-size:8.5px;fill:#888}
 #volbar{text-align:right;margin:-4px 0 10px}
 #volbar button,#bighead button{font-size:11px;padding:2px 9px;cursor:pointer;
   border:1px solid #bbb;border-radius:4px;background:#fff;color:#333}
 #volbar button:hover,#bighead button:hover{background:#f0f6fc;border-color:#7aa7d4}
 #big{display:none;position:fixed;inset:0;z-index:90;background:#fff}
 #bighead{display:flex;align-items:center;gap:18px;padding:10px 18px;
   border-bottom:1px solid #e3e3e3;font-size:12px;color:#444}
 #bighead b{font-size:14px;color:#111}
 #bighead .sp{color:#666} #bighead .hintx{color:#aaa;font-style:italic}
 #bighead input{width:62px;font-size:12px}
 #bighead button{margin-left:auto}
 #bigvol{display:block}
 #bigvol .bg{fill:#c9c9c9;fill-opacity:.5}
 #bigvol .sig{stroke:#000;stroke-width:.4;stroke-opacity:.45}
 #bigvol .thr{stroke:#999;stroke-width:1;stroke-dasharray:5 4}
 #bigvol .ax{stroke:#333;stroke-width:.9}
 #bigvol text.lab{font-size:11px;fill:#111;paint-order:stroke;
   stroke:#fff;stroke-width:3px;stroke-linejoin:round}
 #bigvol text.axl{font-size:12px;fill:#555}
 #legend{margin-top:14px;font-size:11px;color:#666}
 #bar{display:inline-block;width:150px;height:10px;vertical-align:middle;
      border:1px solid #ccc;margin:0 6px}
</style>
<h1>__TITLE__</h1>
<div class="sub">__SUBTITLE__</div>
<div class="ctl">
  <b>annotation</b> <select id="src"></select>
  &nbsp; <b>dataset</b> <select id="ds"></select>
  &nbsp; <b id="lblA">sample</b> <select id="sa"></select>
  <span id="vsWrap">&nbsp; <b>vs</b> <select id="sb"></select></span>
  &nbsp;|&nbsp; organise by <select id="view"></select>
</div>
<div class="ctl">
  rank children by <select id="rank">
    <option value="q">significance (q)</option><option value="shift">|shift|</option>
    <option value="size">gene count</option></select>
  &nbsp; show <select id="topk">
    <option>2</option><option>3</option><option selected>4</option><option>6</option>
    <option>8</option><option>12</option><option value="999">all</option></select> per node
  &nbsp; <label><input type="checkbox" id="dim" checked> dim q &gt;
    <select id="alpha"><option>0.01</option><option selected>0.05</option>
      <option>0.1</option><option>0.25</option></select></label>
</div>
<div id="crumbs"></div>
<svg id="plot"></svg>
<div id="legend">median log<sub>2</sub> shift vs all genes (<span id="base"></span>):
  <span class="dn">down</span><span id="bar"></span><span class="up">up</span>
  &nbsp;·&nbsp; circle area &prop; genes &nbsp;·&nbsp; dashed = no children &nbsp;·&nbsp; faded = n.s.</div>
<div id="grip" title="drag to resize"></div>
<div id="tip"><svg id="vol"></svg>
  <div id="volbar"><button id="expand" title="expand the volcano">&#10530; expand</button></div>
  <div id="tipbody"></div></div>
<div id="big"><div id="bighead">
    <b id="bigtitle"></b>
    <span class="sp">p &lt; <input id="bp" type="number" step="0.01" min="0" max="1" value="0.05"></span>
    <span class="sp">|log<sub>2</sub>FC| &gt; <input id="bfc" type="number" step="0.5" min="0" value="1"></span>
    <span class="sp" id="bigcount"></span>
    <span class="sp hintx">click a point to toggle its label</span>
    <button id="bigclose">close &#10005;</button>
  </div><svg id="bigvol"></svg></div>
<div id="gtip"></div>
<script>
const DATA = __DATA__, GENES = DATA.genes, DS = DATA.datasets;
let ds = Object.keys(DS)[0], src = null, view = null, A = null, B = null;
const paths = {};                                   // breadcrumbs per dataset+view

const NS = 'http://www.w3.org/2000/svg';
const el = (t, at) => { const e = document.createElementNS(NS, t);
                        for (const k in at) e.setAttribute(k, at[k]); return e; };
const opt = (v, t) => Object.assign(document.createElement('option'), {value: v, textContent: t});
const $ = id => document.getElementById(id);

function key(){ return B ? `${A}|${B}` : A; }
function D(){ return DS[ds]; }
function S(){ return D().sources[src]; }
function nodes(){ return S().views[view].nodes; }
function statsOf(){                                  // q lookup, handling reversed pairs
  const st = S().stats, k = key();
  if (st[k]) return [st[k][view] || {}, 1];
  const rev = B ? `${B}|${A}` : null;
  if (rev && st[rev]) return [st[rev][view] || {}, -1];
  return [{}, 1];
}
function baseOf(){
  const d = D(), k = key();
  if (d.base[k] !== undefined) return d.base[k];
  const rev = B ? `${B}|${A}` : null;
  return rev && d.base[rev] !== undefined ? -d.base[rev] : 0;
}
function val(i){                                     // this contrast's value for gene index i
  const d = D(), a = d.fc[A][i];
  if (a === null || a === undefined) return null;
  if (!B) return a;
  const b = d.fc[B][i];
  return (b === null || b === undefined) ? null : a - b;
}
function pval(i){
  const d = D(), a = d.p[A][i];
  if (!B) return a;
  const b = d.p[B][i];
  return (a === null || b === null) ? null : Math.max(a, b);
}
function med(xs){ if(!xs.length) return NaN; const s=[...xs].sort((p,q)=>p-q), m=s.length>>1;
                  return s.length%2 ? s[m] : (s[m-1]+s[m])/2; }

const GIDX = {};              // '<ds>||<src>||<view>' -> {gene index: [node ids, smallest first]}
function geneIndex(){
  const k = ds + '||' + src + '||' + view;
  if (GIDX[k]) return GIDX[k];
  const N = nodes(), idx = {};
  for (const id in N) for (const i of N[id].g) (idx[i] = idx[i] || []).push(id);
  for (const i in idx) idx[i].sort((a, b) => N[a].g.length - N[b].g.length);
  return (GIDX[k] = idx);
}
function pathTo(target){                             // root -> target, via the kids map
  const N = nodes(), root = S().views[view].root, parent = {};
  for (const id in N) for (const c of N[id].kids) parent[c] = id;
  const out = [target];
  let cur = target, guard = 0;
  while (cur !== root && parent[cur] && guard++ < 50){ cur = parent[cur]; out.unshift(cur); }
  return out[0] === root ? out : [root, target];
}

let NSTAT = {};                                      // per-node summary for the current contrast
function recompute(){
  const base = baseOf(), [st] = statsOf(), N = nodes();
  NSTAT = {};
  for (const id in N){
    const vs = N[id].g.map(val).filter(v => v !== null);
    const m = med(vs), s = st[id];
    NSTAT[id] = {n: vs.length, median: m, shift: m - base,
                 up: vs.filter(v => v > 0).length, dn: vs.filter(v => v < 0).length,
                 p: s ? s[0] : null, q: s ? s[1] : null};
  }
}
let LIM = 1, NMAX = 1;
const radius = n => 5 + 10 * Math.sqrt(n / Math.max(NMAX, 1));
function color(v){
  const t = Math.max(-1, Math.min(1, v / LIM));
  const a=[33,102,172], b=[247,247,247], c=[178,24,43];
  const [x,y,f] = t < 0 ? [a,b,t+1] : [b,c,t];
  return `rgb(${x.map((v0,i)=>Math.round(v0+(y[i]-v0)*f)).join(',')})`;
}

function rank(ids, by){
  const k = id => { const n = NSTAT[id];
    if (by === 'q') return [n.q === null ? 1 : n.q, -Math.abs(n.shift)];
    if (by === 'shift') return [-Math.abs(n.shift)];
    return [-n.n]; };
  return ids.slice().sort((p,q) => { const a=k(p), b=k(q);
    for (let i=0;i<a.length;i++) if (a[i]!==b[i]) return a[i]-b[i];
    return nodes()[p].name.localeCompare(nodes()[q].name); });
}

function render(){
  recompute();
  const N = nodes();
  LIM = Math.max(0.3, ...Object.values(NSTAT).map(n => Math.abs(n.shift)).filter(isFinite));
  NMAX = Math.max(...Object.values(NSTAT).map(n => n.n));
  $('base').textContent = (baseOf() >= 0 ? '+' : '') + baseOf().toFixed(2);
  $('bar').style.background = `linear-gradient(to right, ${color(-LIM)}, ${color(0)}, ${color(LIM)})`;

  const pk = ds + '||' + src + '||' + view;
  if (!paths[pk]) paths[pk] = [S().views[view].root];
  const path = paths[pk], rootId = path[path.length - 1];
  const by = $('rank').value, K = +$('topk').value,
        DIM = $('dim').checked, ALPHA = +$('alpha').value;

  const l1 = rank(N[rootId].kids, by).slice(0, K), l2 = {};
  l1.forEach(id => l2[id] = rank(N[id].kids, by).slice(0, K));
  const rowsOf = id => Math.max(1, l2[id].length);
  const total = l1.reduce((s,id)=>s+rowsOf(id), 0) || 1;
  const RMAX = 15, ROW = 2*RMAX + 12, PAD = 34, COLW = [320, 320], LABEL = 300;
  const svg = $('plot');
  svg.setAttribute('width', 70 + COLW[0] + COLW[1] + LABEL);
  svg.setAttribute('height', Math.max(150, total*ROW + PAD*2));
  svg.innerHTML = '';
  const X = [70, 70+COLW[0], 70+COLW[0]+COLW[1]], pos = {};
  let row = 0;
  l1.forEach(id => { const kids = l2[id], span = rowsOf(id);
    kids.forEach((k,i) => pos[k] = [X[2], PAD + (row+i+0.5)*ROW]);
    pos[id] = [X[1], PAD + (row + span/2)*ROW]; row += span; });
  pos[rootId] = [X[0], PAD + total*ROW/2];
  const link = (a,b) => { const [x1,y1]=pos[a], [x2,y2]=pos[b], mx=(x1+x2)/2;
    svg.appendChild(el('path',{class:'link', d:`M${x1},${y1} H${mx} V${y2} H${x2}`})); };
  l1.forEach(id => { link(rootId,id); l2[id].forEach(k => link(id,k)); });

  const draw = (id, isRoot) => {
    const n = NSTAT[id], [x,y] = pos[id], leaf = N[id].kids.length === 0;
    const ns = DIM && n.q !== null && n.q > ALPHA, r = radius(n.n);
    const g = el('g', {class:'node'});
    const c = el('circle', {cx:x, cy:y, r:r, fill: color(n.shift), class: leaf?'leaf':'',
                            'fill-opacity': ns?0.35:1, stroke: ns?'#aaa':'#333'});
    const lab = `${N[id].name} (${n.n})`;
    const t = el('text', {x: x+r+7, y: y, class: leaf?'leaf':''});
    t.textContent = (lab.length>40 ? lab.slice(0,39)+'…' : lab) + (ns?'  n.s.':'');
    if (ns) t.setAttribute('fill', '#999');
    t.appendChild(el('title', {})).textContent = lab;
    g.appendChild(c); g.appendChild(t);
    if (!leaf && !isRoot) c.onclick = t.onclick = () => { path.push(id); render(); };
    if (isRoot && path.length > 1) c.onclick = t.onclick = () => { path.pop(); render(); };
    g.onmouseenter = () => tip(id);
    svg.appendChild(g);
  };
  draw(rootId, true);
  l1.forEach(id => { draw(id,false); l2[id].forEach(k => draw(k,false)); });

  const cr = $('crumbs'); cr.innerHTML = '';
  path.forEach((id,i) => {
    if (i) cr.appendChild(Object.assign(document.createElement('span'),
                                        {className:'sep', textContent:'›'}));
    const s = document.createElement('span'); s.textContent = N[id].name;
    if (i === path.length-1) s.className = 'cur';
    else { s.className='c'; s.onclick = () => { paths[pk] = path.slice(0,i+1); render(); }; }
    cr.appendChild(s);
  });
  drawVolcano();
}

// ---------------- volcano ----------------
const VH = 210, VM = {l:30, r:8, t:10, b:24};
const panelW = () => $('tip').clientWidth - 32;
let CUR = {hl:null, n:0};
function drawVolcano(hl, nlabel){
  if (hl !== undefined) CUR = {hl:hl, n:nlabel};
  hl = CUR.hl; nlabel = CUR.n;
  const VW = Math.max(180, panelW());
  const pts = [];
  for (let i = 0; i < GENES.length; i++){
    const v = val(i); if (v === null) continue;
    const p = pval(i);
    pts.push([i, v, p === null ? 0 : -Math.log10(Math.max(p, 1e-300))]);
  }
  if (!pts.length) return;
  const xmin = Math.min(...pts.map(p=>p[1])), xmax = Math.max(...pts.map(p=>p[1]));
  const ymax = Math.max(...pts.map(p=>p[2]), 1);
  const vx = v => VM.l + (v - xmin)/(xmax - xmin || 1)*(VW - VM.l - VM.r);
  const vy = v => VH - VM.b - v/ymax*(VH - VM.t - VM.b);
  const svg = $('vol'); svg.setAttribute('width', VW); svg.setAttribute('height', VH);
  svg.innerHTML = '';
  svg.appendChild(el('line',{class:'ax', x1:vx(0), x2:vx(0), y1:VM.t, y2:VH-VM.b}));
  svg.appendChild(el('line',{class:'ax', x1:VM.l, x2:VW-VM.r, y1:VH-VM.b, y2:VH-VM.b}));
  const marks = [];
  VPTS = [];
  pts.forEach(([i,v,y]) => {
    const X = vx(v), Y = vy(y);
    VPTS.push({x: X, y: Y, i: i});
    if (hl && hl.has(i)) marks.push([X, Y, GENES[i], v, y]);
    else svg.appendChild(el('circle',{cx:X, cy:Y, r:1.5, class:'bg'}));
  });
  marks.forEach(([x,y,g,v]) => svg.appendChild(
    el('circle',{cx:x, cy:y, r:2.8, class:'hl', fill: v>0?'#b2182b':'#2166ac'})));
  const placed = [];
  marks.sort((a,b)=>b[4]-a[4]).slice(0, nlabel||8).forEach(([x,y,g]) => {
    const w = g.length*5.2 + 6;
    let lx = Math.min(x+4, VW-w-2), ly = y-4, dir = -1, step = 0;
    const hits = (a,b) => placed.some(q => Math.abs(b-q.y) < 8 && a < q.x+q.w && q.x < a+w);
    while (hits(lx, ly) && step < 12){ step++; dir = -dir; ly = y-4 + dir*Math.ceil(step/2)*8; }
    placed.push({x:lx, y:ly, w:w});
    if (Math.abs(ly - (y-4)) > 5)
      svg.appendChild(el('line',{x1:x, y1:y, x2:lx-2, y2:ly-3,
                                 stroke:'#9ca3af','stroke-width':.6}));
    const t = el('text',{x: lx, y: ly, class:'lab'});
    t.textContent = g; svg.appendChild(t);
  });
  [[xmin,0],[xmax,1]].forEach(([v,i]) => {
    const t = el('text',{x: i?VW-VM.r:VM.l, y: VH-12, class:'axl','text-anchor': i?'end':'start'});
    t.textContent = v.toFixed(1); svg.appendChild(t); });
  let t = el('text',{x:2, y:VM.t+8, class:'axl'}); t.textContent='-log10 p'; svg.appendChild(t);
  t = el('text',{x:(VM.l+VW)/2, y:VH-2, class:'axl','text-anchor':'middle'});
  t.textContent = B ? `log2 ${A} - ${B}`
                    : (D().pairs === false ? `log2 ${A}` : `log2 ${A} / control`);
  svg.appendChild(t);
}

let VPTS = [];                                       // drawn volcano points, screen coords
function nearestPoint(ev){
  const r = $('vol').getBoundingClientRect();
  const mx = ev.clientX - r.left, my = ev.clientY - r.top;
  let best = null, bd = 49;                          // within 7 px
  for (const p of VPTS){
    const d = (p.x - mx) ** 2 + (p.y - my) ** 2;
    if (d < bd){ bd = d; best = p; }
  }
  return best;
}
(function(){
  const svg = $('vol'), lab = $('gtip');
  svg.addEventListener('mousemove', ev => {
    const p = nearestPoint(ev);
    if (!p){ lab.style.display = 'none'; return; }
    const v = val(p.i), q = pval(p.i);
    lab.textContent = `${GENES[p.i]}  ${fmt(v)}` + (q === null ? '' : `  p=${q.toExponential(1)}`);
    lab.style.display = 'block';
    lab.style.left = (ev.clientX + 12) + 'px';
    lab.style.top  = (ev.clientY - 26) + 'px';
  });
  svg.addEventListener('mouseleave', () => { lab.style.display = 'none'; });
  svg.addEventListener('click', ev => {
    const p = nearestPoint(ev); if (!p) return;
    const ids = geneIndex()[p.i];
    if (!ids || !ids.length){
      lab.textContent = `${GENES[p.i]} — not annotated in this view`;
      return;
    }
    // the most specific node is often a 3-gene leaf; the next one up gives context
    const target = ids.length > 1 ? ids[1] : ids[0];
    paths[ds + '||' + src + '||' + view] = pathTo(target);
    render();
    tip(target);
  });
})();

// ---------------- full-screen volcano ----------------
let BPTS = [], BLAB = new Set(), BOPEN = false;

function bigPoints(){
  const out = [];
  for (let i = 0; i < GENES.length; i++){
    const v = val(i); if (v === null) continue;
    const p = pval(i);
    out.push({i: i, v: v, y: p === null ? 0 : -Math.log10(Math.max(p, 1e-300)), p: p});
  }
  return out;
}
function openBig(){
  BOPEN = true;
  const pts = bigPoints();
  const pc = +$('bp').value, fc = +$('bfc').value;
  const sig = pts.filter(d => d.p !== null && d.p < pc && Math.abs(d.v) > fc);
  BLAB = new Set(sig.slice().sort((a,b) => b.y - a.y).slice(0, 12).map(d => d.i));
  $('big').style.display = 'block';
  $('bigtitle').textContent = B ? `${A} − ${B}` : A;
  drawBig();
}
function closeBig(){ BOPEN = false; $('big').style.display = 'none'; $('gtip').style.display = 'none'; }

function drawBig(){
  const svg = $('bigvol'), W = innerWidth, H = innerHeight - $('bighead').offsetHeight - 4;
  const M = {l: 72, r: 40, t: 28, b: 56};
  svg.setAttribute('width', W); svg.setAttribute('height', H);
  svg.innerHTML = '';
  const pts = bigPoints();
  if (!pts.length) return;
  const pc = +$('bp').value, fc = +$('bfc').value;
  const xlo = Math.min(...pts.map(d => d.v)), xhi = Math.max(...pts.map(d => d.v));
  const yhi = Math.max(...pts.map(d => d.y), 1);
  const pad = (xhi - xlo) * 0.04 || 1;
  const X = v => M.l + (v - xlo + pad) / (xhi - xlo + 2*pad) * (W - M.l - M.r);
  const Y = v => H - M.b - v / (yhi * 1.06) * (H - M.t - M.b);

  // threshold guides
  if (pc > 0 && pc < 1){
    const yl = Y(-Math.log10(pc));
    svg.appendChild(el('line', {class:'thr', x1:M.l, x2:W-M.r, y1:yl, y2:yl}));
    const t = el('text', {x: W-M.r, y: yl-5, class:'axl', 'text-anchor':'end'});
    t.textContent = `p = ${pc}`; svg.appendChild(t);
  }
  [-fc, fc].forEach(v => { if (!fc) return;
    svg.appendChild(el('line', {class:'thr', x1:X(v), x2:X(v), y1:M.t, y2:H-M.b})); });
  svg.appendChild(el('line', {class:'ax', x1:M.l, x2:W-M.r, y1:H-M.b, y2:H-M.b}));
  svg.appendChild(el('line', {class:'ax', x1:M.l, x2:M.l, y1:M.t, y2:H-M.b}));

  BPTS = [];
  let nUp = 0, nDn = 0;
  pts.forEach(d => {
    const x = X(d.v), y = Y(d.y);
    BPTS.push({x: x, y: y, i: d.i});
    const isSig = d.p !== null && d.p < pc && Math.abs(d.v) > fc;
    if (isSig) (d.v > 0 ? nUp++ : nDn++);
    svg.appendChild(el('circle', {cx:x, cy:y, r: isSig ? 3.4 : 2,
      fill: isSig ? (d.v > 0 ? '#b2182b' : '#2166ac') : '#c9c9c9',
      class: isSig ? 'sig' : 'bg'}));
  });
  $('bigcount').textContent = `${nUp} up / ${nDn} down of ${pts.length}`;

  const placed = [];
  pts.filter(d => BLAB.has(d.i)).sort((a,b) => b.y - a.y).forEach(d => {
    const px = X(d.v), py = Y(d.y), g = GENES[d.i], w = g.length*6.2 + 8;
    let lx = px + 7, ly = py - 6, dir = -1, step = 0;
    const hits = (a,b) => placed.some(q => Math.abs(b-q.y) < 12 && a < q.x+q.w && q.x < a+w);
    while (hits(lx, ly) && step < 14){ step++; dir = -dir; ly = py-6 + dir*Math.ceil(step/2)*12; }
    if (hits(lx, ly)){ lx = px - w - 7; ly = py - 6; }
    placed.push({x:lx, y:ly, w:w});
    if (Math.abs(ly - (py-6)) > 6 || lx < px)
      svg.appendChild(el('line', {x1:px, y1:py, x2: lx < px ? lx+w+2 : lx-2, y2: ly-3,
                                  stroke:'#9ca3af','stroke-width':.7}));
    const t = el('text', {x: lx, y: ly, class:'lab'});
    t.textContent = g; svg.appendChild(t);
  });

  let t = el('text', {x: (M.l + W - M.r)/2, y: H - 16, class:'axl', 'text-anchor':'middle'});
  t.textContent = B ? `log2 ${A} − ${B}` : `log2 ${A}`; svg.appendChild(t);
  t = el('text', {x: 18, y: (M.t + H - M.b)/2, class:'axl', 'text-anchor':'middle',
                  transform: `rotate(-90 18 ${(M.t + H - M.b)/2})`});
  t.textContent = '-log10 p'; svg.appendChild(t);
  [xlo, 0, xhi].forEach(v => {
    const q = el('text', {x: X(v), y: H - M.b + 18, class:'axl', 'text-anchor':'middle'});
    q.textContent = v.toFixed(1); svg.appendChild(q); });
  for (let k = 0; k <= Math.floor(yhi); k += Math.max(1, Math.ceil(yhi/8))){
    const q = el('text', {x: M.l - 10, y: Y(k) + 4, class:'axl', 'text-anchor':'end'});
    q.textContent = k; svg.appendChild(q); }
}
function nearestBig(ev){
  let best = null, bd = 121;                         // within 11 px
  for (const p of BPTS){
    const d = (p.x - ev.clientX) ** 2 + (p.y - (ev.clientY - $('bighead').offsetHeight)) ** 2;
    if (d < bd){ bd = d; best = p; }
  }
  return best;
}
$('expand').onclick = openBig;
$('bigclose').onclick = closeBig;
$('bp').oninput = $('bfc').oninput = () => { if (BOPEN) drawBig(); };
addEventListener('keydown', e => { if (e.key === 'Escape' && BOPEN) closeBig(); });
addEventListener('resize', () => { if (BOPEN) drawBig(); });
$('bigvol').addEventListener('mousemove', ev => {
  const p = nearestBig(ev), lab = $('gtip');
  if (!p){ lab.style.display = 'none'; return; }
  const v = val(p.i), q = pval(p.i);
  lab.textContent = `${GENES[p.i]}  ${fmt(v)}` + (q === null ? '' : `  p=${q.toExponential(1)}`);
  lab.style.display = 'block';
  lab.style.left = (ev.clientX + 12) + 'px'; lab.style.top = (ev.clientY - 26) + 'px';
});
$('bigvol').addEventListener('mouseleave', () => { $('gtip').style.display = 'none'; });
$('bigvol').addEventListener('click', ev => {
  const p = nearestBig(ev); if (!p) return;
  BLAB.has(p.i) ? BLAB.delete(p.i) : BLAB.add(p.i);
  drawBig();
});

function fmt(v){ return (v>=0?'+':'') + v.toFixed(2); }
function tip(id){
  const n = NSTAT[id], N = nodes(), box = $('tipbody');
  drawVolcano(new Set(N[id].g), 8);
  const rows = N[id].g.map(i => [i, val(i)]).filter(r => r[1] !== null)
    .sort((a,b) => b[1]-a[1]).map(([i,v]) => {
      const p = pval(i);
      return `<tr><td>${GENES[i]}</td><td class="v ${v>0?'up':'dn'}">${fmt(v)}</td>` +
             `<td class="v">${p===null?'':p.toExponential(1)}</td></tr>`; }).join('');
  const stat = n.q === null ? '' :
    ` &nbsp;·&nbsp; p = ${n.p.toExponential(1)} &nbsp; q = ${n.q.toExponential(1)}`;
  box.innerHTML = `<h3>${N[id].name}</h3><div class="meta">${id} &nbsp;·&nbsp; ${n.n} genes` +
    ` &nbsp;·&nbsp; median ${fmt(n.median)} (shift ${fmt(n.shift)}) &nbsp;·&nbsp; ` +
    `${n.up} up / ${n.dn} down${stat}<br><i>${B ? A+' − '+B : A+' vs control'}</i></div>` +
    `<table><thead><tr><th>gene</th><th class="v">log2</th><th class="v">p</th></tr></thead>` +
    `<tbody>${rows}</tbody></table>`;
  $('tip').scrollTop = 0;
}
function placeholder(){ drawVolcano(null, 0);
  $('tipbody').innerHTML = '<div class="hint">hover a node to highlight its genes</div>'; }

// ---------------- selectors ----------------
function fillSamples(){
  const smp = D().samples, canPair = D().pairs !== false;
  $('sa').innerHTML = ''; smp.forEach(s => $('sa').appendChild(opt(s, s)));
  if (!smp.includes(A)) A = smp[0];
  $('sa').value = A;
  $('lblA').textContent = canPair ? 'sample' : 'comparison';
  $('vsWrap').style.display = canPair ? '' : 'none';
  if (!canPair){ B = null; return; }
  $('sb').innerHTML = ''; $('sb').appendChild(opt('', '— none (vs control) —'));
  smp.filter(s => s !== A).forEach(s => $('sb').appendChild(opt(s, s)));
  if (B === A || !smp.includes(B)) B = null;
  $('sb').value = B || '';
}
function fillSources(){
  const ss = Object.keys(D().sources);
  $('src').innerHTML = ''; ss.forEach(x => $('src').appendChild(opt(x, x)));
  if (!ss.includes(src)) src = ss[0];
  $('src').value = src;
}
function fillViews(){
  const vs = Object.keys(S().views);
  $('view').innerHTML = ''; vs.forEach(v => $('view').appendChild(opt(v, v)));
  if (!vs.includes(view)) view = vs[0];
  $('view').value = view;
}
Object.keys(DS).forEach(k => $('ds').appendChild(opt(k, k)));
$('ds').onchange = () => { ds = $('ds').value; fillSources(); fillViews(); fillSamples();
                           placeholder(); render(); };
$('src').onchange = () => { src = $('src').value; fillViews(); placeholder(); render(); };
$('sa').onchange = () => { A = $('sa').value; fillSamples(); placeholder(); render(); };
$('sb').onchange = () => { B = $('sb').value || null; placeholder(); render(); };
$('view').onchange = () => { view = $('view').value; placeholder(); render(); };
['rank','topk','dim','alpha'].forEach(id => $(id).onchange = render);

(function(){
  const grip = $('grip'), tip = $('tip');
  const MINW = 240, MAXW = () => Math.min(900, innerWidth - 320);
  function setW(w){ w = Math.max(MINW, Math.min(MAXW(), w));
    tip.style.width = w+'px'; grip.style.right = w+'px';
    document.body.style.paddingRight = (w+22)+'px'; drawVolcano(); }
  grip.addEventListener('mousedown', e => { e.preventDefault();
    grip.classList.add('drag'); document.body.classList.add('resizing');
    const mv = ev => setW(innerWidth - ev.clientX);
    const up = () => { grip.classList.remove('drag');
      document.body.classList.remove('resizing');
      removeEventListener('mousemove', mv); removeEventListener('mouseup', up); };
    addEventListener('mousemove', mv); addEventListener('mouseup', up); });
  grip.addEventListener('dblclick', () => setW(__PANELW__));
  addEventListener('resize', () => drawVolcano());
})();

fillSources(); fillViews(); fillSamples(); placeholder(); render();
</script>
"""


class MultiExplorer:
    def __init__(self, html, payload, path=None, height=800, width="100%"):
        self.html, self.payload, self.path = html, payload, path
        self.height, self.width = height, width

    def _repr_html_(self):
        doc = _html.escape(self.html, quote=True)
        return (f'<iframe srcdoc="{doc}" width="{self.width}" height="{self.height}" '
                f'style="border:1px solid #ddd;border-radius:6px;resize:vertical" '
                f'sandbox="allow-scripts"></iframe>')

    def save(self, path):
        open(path, "w").write(self.html)
        self.path = path
        return self

    def __repr__(self):
        bits = ", ".join(f"{k}: {len(v['samples'])} samples, {len(v['base'])} contrasts, "
                         f"{len(v['sources'])} sources"
                         for k, v in self.payload["datasets"].items())
        return f"<MultiExplorer [{bits}]{' -> ' + self.path if self.path else ''}>"


def multi_explorer(datasets, out=None, title=None, subtitle=None, height=800,
                   width="100%", panel_width=380, **kw) -> MultiExplorer:
    """Build the multi-dataset explorer. `datasets` maps a label to
    dict(raw=<Proteins sheet>, denom=<control name or None>, order=[...])."""
    payload = build_multi(datasets, **kw)
    title = title or "GO explorer"
    n = sum(len(d["base"]) for d in payload["datasets"].values())
    srcs = sorted({s for d in payload["datasets"].values() for s in d["sources"]})
    sub = subtitle or (f"{len(payload['datasets'])} datasets · {n} contrasts · "
                       f"{' / '.join(srcs)} · {len(payload['genes'])} genes · "
                       f"pick a sample, or two to compare")
    html = (_TPL.replace("__DATA__", _json.dumps(payload, separators=(",", ":")))
                .replace("__TITLE__", title).replace("__SUBTITLE__", sub)
                .replace("__PANELW__", str(int(panel_width))))
    ex = MultiExplorer(html, payload, height=height, width=width)
    if out:
        ex.save(out)
        print(f"wrote {out}  ({len(html)//1024} KB)")
    return ex
