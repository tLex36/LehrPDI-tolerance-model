"""Self-contained interactive GO/SynGO explorer: drill down two layers at a time.

    from go_explorer import explorer_html
    explorer_html(contrast_df, source="SynGO", out="syngo_explorer.html")

Renders the root term, its most enriched children, and their most enriched children.
Click any node to re-root on it; click a breadcrumb to go back. Hover shows the full
gene table. No CDN or external library -- the file opens offline.
"""
import html as _html
import json as _json
import os as _os

import pandas as pd

import go as _go


def _rank_key(stats_row, by):
    """Sort key for 'most enriched': lower is better."""
    if by == "q":
        return (stats_row.get("q") if stats_row.get("q") is not None else 1.0,
                -abs(stats_row.get("shift") or 0.0))
    if by == "shift":
        return (-abs(stats_row.get("shift") or 0.0),)
    if by == "size":
        return (-stats_row.get("n", 0),)
    raise ValueError(f"rank_by must be 'q', 'shift' or 'size'; got {by!r}")


def build_model(d, source="SynGO", organism="mgi", root=None, aspect=None, slim=None,
                min_size=3, max_size=1000, max_depth=12, stats=None, center=None,
                **kw) -> dict:
    """{nodes, genes, root} -- the whole hierarchy as plain JSON-able data."""
    obo = _go.load_obo()
    root = root or _go.SYNAPSE_ROOT
    sets = _go.term_genes(min_size=min_size, max_size=max_size, organism=organism,
                          source=source, aspect=aspect, slim=slim, **kw)
    if root not in sets:
        full = _go.term_genes(min_size=1, max_size=10**9, organism=organism,
                              source=source, aspect=aspect, **kw)
        if root not in full:
            raise ValueError(f"root {root} has no annotations here")
        sets[root] = full[root]

    vals = d.set_index("gene").log2fc.groupby(level=0).mean()
    pv = (d.set_index("gene").pval.groupby(level=0).min()
          if "pval" in d.columns else None)
    sets = {t: (g & set(vals.index)) for t, g in sets.items()}
    sets = {t: g for t, g in sets.items() if len(g) >= min_size or t == root}

    tree = _go._tree(root, sets, obo, max_depth=max_depth)
    base = float(vals.median()) if center is None else float(center)
    qmap = dict(zip(stats.go_id, stats.q_value)) if stats is not None else {}
    pmap = dict(zip(stats.go_id, stats.p_value)) if stats is not None else {}

    nodes = {}
    for t, kids in tree.items():
        g = sorted(sets[t], key=lambda s: -vals[s])
        v = vals.loc[g]
        nodes[t] = {
            "id": t, "name": obo[t]["name"], "kids": [k for k in kids if k in sets],
            "n": len(g), "median": round(float(v.median()), 3),
            "shift": round(float(v.median()) - base, 3),
            "up": int((v > 0).sum()), "down": int((v < 0).sum()),
            "q": (None if t not in qmap else float(qmap[t])),
            "p": (None if t not in pmap else float(pmap[t])),
            "g": g,
        }
    genes = {s: [round(float(vals[s]), 3),
                 (None if pv is None or s not in pv.index else float(pv[s]))]
             for s in sorted(set().union(*(set(n["g"]) for n in nodes.values())))}
    return {"root": root, "nodes": nodes, "genes": genes, "base": round(base, 3)}


_TEMPLATE = r"""<!doctype html>
<meta charset="utf-8">
<title>__TITLE__</title>
<style>
 body{font:13px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
      margin:0;padding:18px 22px;color:#222;background:#fff}
 h1{font-size:16px;font-weight:600;margin:0 0 4px}
 .sub{color:#777;font-size:12px;margin-bottom:12px}
 #crumbs{margin:10px 0 6px;font-size:12px}
 #crumbs span.c{color:#1a6fb5;cursor:pointer;text-decoration:none}
 #crumbs span.c:hover{text-decoration:underline}
 #crumbs span.sep{color:#bbb;margin:0 6px}
 #crumbs span.cur{font-weight:600;color:#222}
 svg{display:block}
 .link{fill:none;stroke:#bbb;stroke-width:1.1}
 .node circle{stroke:#333;stroke-width:.8;cursor:pointer}
 .node circle.leaf{stroke-dasharray:2 2;cursor:default}
 .node text{font-size:11.5px;cursor:pointer;dominant-baseline:middle}
 .node text.leaf{cursor:default;fill:#555}
 .node:hover circle{stroke-width:2}
 #tip{position:fixed;top:0;right:0;width:__PANELW__px;height:100vh;z-index:50;
      overflow:auto;background:#fcfcfc;border-left:1px solid #ddd;padding:16px 16px 28px;
      font-size:11.5px;box-sizing:border-box}
 #tip h3{margin:0 0 2px;font-size:13px}
 #tip .meta{color:#666;margin-bottom:9px;line-height:1.5}
 #tip .hint{color:#aaa;font-style:italic;text-align:center;padding:30px 0}
 #grip{position:fixed;top:0;right:__PANELW__px;width:7px;height:100vh;z-index:60;
       cursor:col-resize;background:transparent}
 #grip:hover,#grip.drag{background:#cfe3f5}
 body.resizing{user-select:none;cursor:col-resize}
 #vol{display:block;margin:0 0 10px}
 #vol .bg{fill:#d6d6d6;fill-opacity:.55}
 #vol .hl{stroke:#000;stroke-width:.5}
 #vol .ax{stroke:#bbb;stroke-width:.7}
 #vol text.lab{font-size:8.5px;fill:#222}
 #vol text.axl{font-size:8.5px;fill:#888}
 body{padding-right:calc(__PANELW__px + 22px)}
 #tip table{border-collapse:collapse;width:100%}
 #tip th{text-align:left;font-weight:600;border-bottom:1px solid #ddd;padding:2px 8px 3px 0;
         position:sticky;top:0;background:#fff}
 #tip td{padding:1px 8px 1px 0;font-variant-numeric:tabular-nums}
 #tip td.v{text-align:right}
 .up{color:#b2182b}.dn{color:#2166ac}
 #legend{margin-top:14px;font-size:11px;color:#666}
 #bar{display:inline-block;width:150px;height:10px;vertical-align:middle;
      border:1px solid #ccc;margin:0 6px}
 .ctl{font-size:12px;color:#555;margin:6px 0 0}
 .ctl select{font-size:12px}
</style>
<h1>__TITLE__</h1>
<div class="sub">__SUBTITLE__</div>
<div class="ctl">organise by
  <select id="view"></select>
  &nbsp;|&nbsp; rank children by
  <select id="rank">
    <option value="q">significance (q)</option>
    <option value="shift">|shift|</option>
    <option value="size">gene count</option>
  </select>
  &nbsp; show
  <select id="topk">
    <option>2</option><option>3</option><option selected>4</option><option>6</option>
    <option>8</option><option>12</option><option value="999">all</option>
  </select>
  children per node
  &nbsp; <label><input type="checkbox" id="dim" checked> dim q &gt;
    <select id="alpha"><option>0.01</option><option selected>0.05</option>
      <option>0.1</option><option>0.25</option></select></label>
</div>
<div id="crumbs"></div>
<svg id="plot"></svg>
<div id="legend">median log<sub>2</sub>FC shift vs all genes (<span id="base">__BASE__</span>):
  <span class="dn">down</span><span id="bar"></span><span class="up">up</span>
  &nbsp;·&nbsp; circle area ∝ genes &nbsp;·&nbsp; dashed outline = no children &nbsp;·&nbsp; faded = n.s.
</div>
<div id="grip" title="drag to resize"></div>
<div id="tip"><svg id="vol"></svg><div id="tipbody"></div></div>
<script>
const DATA = __DATA__;
const MODELS = DATA.views;                     // {view label: model}
const PTS = DATA.points;                       // [[gene, log2fc, p], ...] all detected
const VIEWS = Object.keys(MODELS);
let view = VIEWS[0];
const paths = {};                              // breadcrumb trail per view
VIEWS.forEach(v => paths[v] = [MODELS[v].root]);

let MODEL = MODELS[view], N = MODEL.nodes, G = MODEL.genes, path = paths[view];

const sel = document.getElementById('view');
VIEWS.forEach(v => sel.appendChild(Object.assign(document.createElement('option'),
                                                 {value: v, textContent: v})));
sel.onchange = () => {
  view = sel.value;
  MODEL = MODELS[view]; N = MODEL.nodes; G = MODEL.genes; path = paths[view];
  document.getElementById('base').textContent = (MODEL.base >= 0 ? '+' : '') +
                                                MODEL.base.toFixed(2);
  placeholder();
  render();
};

const LIM = Math.max(0.3, ...VIEWS.flatMap(v =>
              Object.values(MODELS[v].nodes).map(n => Math.abs(n.shift))));
const NMAX = Math.max(...VIEWS.flatMap(v =>
              Object.values(MODELS[v].nodes).map(n => n.n)));
const RMIN = 5, RMAX = 15;
const radius = n => RMIN + (RMAX - RMIN) * Math.sqrt(n / NMAX);
function color(v){                             // blue -> white -> red, centred on 0
  const t = Math.max(-1, Math.min(1, v / LIM));
  const a = [33,102,172], b = [247,247,247], c = [178,24,43];
  const [x,y,f] = t < 0 ? [a,b,t+1] : [b,c,t];
  return `rgb(${x.map((v0,i)=>Math.round(v0+(y[i]-v0)*f)).join(',')})`;
}
document.getElementById('bar').style.background =
  `linear-gradient(to right, ${color(-LIM)}, ${color(0)}, ${color(LIM)})`;

function rank(ids, by){
  const key = id => {
    const n = N[id];
    if (by === 'q')     return [n.q === null ? 1 : n.q, -Math.abs(n.shift)];
    if (by === 'shift') return [-Math.abs(n.shift)];
    return [-n.n];
  };
  return ids.slice().sort((p,q) => {
    const a = key(p), b = key(q);
    for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return a[i] - b[i];
    return N[p].name.localeCompare(N[q].name);
  });
}

function render(){
  const by = document.getElementById('rank').value;
  const K  = +document.getElementById('topk').value;
  const DIM = document.getElementById('dim').checked;
  const ALPHA = +document.getElementById('alpha').value;
  path = paths[view];
  const rootId = path[path.length - 1];

  // two layers below the current root
  const l1 = rank(N[rootId].kids, by).slice(0, K);
  const l2 = {};
  l1.forEach(id => l2[id] = rank(N[id].kids, by).slice(0, K));

  const rowsOf = id => Math.max(1, l2[id].length);
  const total  = l1.reduce((s,id) => s + rowsOf(id), 0) || 1;
  const ROW = 2 * RMAX + 12, PAD = 34, COLW = [320, 320], LABEL = 300;
  const H = Math.max(150, total * ROW + PAD * 2);
  const W = 70 + COLW[0] + COLW[1] + LABEL;

  const svg = document.getElementById('plot');
  svg.setAttribute('width', W); svg.setAttribute('height', H);
  svg.innerHTML = '';
  const NS = 'http://www.w3.org/2000/svg';
  const el = (t, at) => { const e = document.createElementNS(NS, t);
                          for (const k in at) e.setAttribute(k, at[k]); return e; };

  const X = [70, 70 + COLW[0], 70 + COLW[0] + COLW[1]];
  const pos = {};
  let row = 0;
  l1.forEach(id => {
    const kids = l2[id], span = rowsOf(id);
    kids.forEach((k, i) => pos[k] = [X[2], PAD + (row + i + 0.5) * ROW]);
    pos[id] = [X[1], PAD + (row + span / 2) * ROW];
    row += span;
  });
  pos[rootId] = [X[0], PAD + total * ROW / 2];

  const link = (a, b) => {
    const [x1,y1] = pos[a], [x2,y2] = pos[b], mx = (x1 + x2) / 2;
    svg.appendChild(el('path', {class:'link', d:`M${x1},${y1} H${mx} V${y2} H${x2}`}));
  };
  l1.forEach(id => { link(rootId, id); l2[id].forEach(k => link(id, k)); });

  const draw = (id, isRoot) => {
    const n = N[id], [x,y] = pos[id];
    const g = el('g', {class:'node'});
    const leaf = n.kids.length === 0;
    const r = radius(n.n);
    const ns = DIM && n.q !== null && n.q > ALPHA;       // not significant at ALPHA
    const c = el('circle', {cx:x, cy:y, r: r, fill: color(n.shift),
                            class: leaf ? 'leaf' : '',
                            'fill-opacity': ns ? 0.35 : 1,
                            stroke: ns ? '#aaa' : '#333'});
    const label = `${n.name} (${n.n})`;
    const t = el('text', {x: x + r + 7, y: y,
                          class: leaf ? 'leaf' : ''});
    t.textContent = (label.length > 40 ? label.slice(0, 39) + '…' : label)
                    + (ns ? '  n.s.' : '');
    if (ns) t.setAttribute('fill', '#999');
    t.appendChild(el('title', {})).textContent = label;
    g.appendChild(c); g.appendChild(t);
    if (!leaf && !isRoot) { c.onclick = t.onclick = () => { paths[view].push(id);
                                                        path = paths[view]; render(); }; }
    if (isRoot && path.length > 1) { c.onclick = t.onclick = () => { paths[view].pop();
                                                                 path = paths[view]; render(); }; }
    g.onmouseenter = e => tip(e, id);
    g.onmouseleave = hide;
    svg.appendChild(g);
  };
  draw(rootId, true);
  l1.forEach(id => { draw(id, false); l2[id].forEach(k => draw(k, false)); });

  const cr = document.getElementById('crumbs');
  cr.innerHTML = '';
  path.forEach((id, i) => {
    if (i) cr.appendChild(Object.assign(document.createElement('span'),
                                        {className:'sep', textContent:'›'}));
    const s = document.createElement('span');
    s.textContent = N[id].name;
    if (i === path.length - 1) s.className = 'cur';
    else { s.className = 'c'; s.onclick = () => { paths[view] = path.slice(0, i + 1);
                                                  path = paths[view]; render(); }; }
    cr.appendChild(s);
  });
}

const VH = 210, VM = {l:30, r:8, t:10, b:24};
const panelW = () => document.getElementById('tip').clientWidth - 32;
const NSV = 'http://www.w3.org/2000/svg';
const sel2 = (t, at) => { const e = document.createElementNS(NSV, t);
                          for (const k in at) e.setAttribute(k, at[k]); return e; };
const XS = PTS.map(p => p[1]), YS = PTS.map(p => p[2] === null ? 0 : -Math.log10(Math.max(p[2], 1e-300)));
const XMIN = Math.min(...XS), XMAX = Math.max(...XS), YMAX = Math.max(...YS, 1);
const vy = v => VH - VM.b - v / YMAX * (VH - VM.t - VM.b);
let CUR = {hl: null, n: 0};                    // last highlight, so a resize can redraw

function drawVolcano(hl, nlabel){
  if (hl !== undefined) CUR = {hl: hl, n: nlabel};
  hl = CUR.hl; nlabel = CUR.n;
  const VW = Math.max(180, panelW());
  const vx = v => VM.l + (v - XMIN) / (XMAX - XMIN || 1) * (VW - VM.l - VM.r);
  const svg = document.getElementById('vol');
  svg.setAttribute('width', VW); svg.setAttribute('height', VH);
  svg.innerHTML = '';
  svg.appendChild(sel2('line', {class:'ax', x1:vx(0), x2:vx(0), y1:VM.t, y2:VH - VM.b}));
  svg.appendChild(sel2('line', {class:'ax', x1:VM.l, x2:VW - VM.r,
                                y1:VH - VM.b, y2:VH - VM.b}));
  const bg = [], marks = [];
  PTS.forEach((pt, i) => {
    const [g] = pt, inSet = hl && hl.has(g);
    (inSet ? marks : bg).push([vx(XS[i]), vy(YS[i]), g, XS[i], YS[i]]);
  });
  bg.forEach(([x,y]) => svg.appendChild(sel2('circle', {cx:x, cy:y, r:1.5, class:'bg'})));
  marks.forEach(([x,y,g,fc]) =>
    svg.appendChild(sel2('circle', {cx:x, cy:y, r:2.8, class:'hl',
                                    fill: fc > 0 ? '#b2182b' : '#2166ac'})));
  // label the most significant highlighted genes, nudging apart vertical collisions
  const used = [];
  marks.sort((a,b) => b[4] - a[4]).slice(0, nlabel || 8).forEach(([x,y,g]) => {
    let ly = y - 4;
    while (used.some(u => Math.abs(u - ly) < 8)) ly -= 8;
    used.push(ly);
    const t = sel2('text', {x: Math.min(x + 4, VW - 34), y: ly, class:'lab'});
    t.textContent = g;
    svg.appendChild(t);
  });
  [[XMIN,'left'],[XMAX,'end']].forEach(([v, anc], i) => {
    const t = sel2('text', {x: i ? VW - VM.r : VM.l, y: VH - 12, class:'axl',
                            'text-anchor': i ? 'end' : 'start'});
    t.textContent = v.toFixed(1); svg.appendChild(t);
  });
  const yl = sel2('text', {x: 2, y: VM.t + 8, class:'axl'});
  yl.textContent = '-log10 p'; svg.appendChild(yl);
  const xl = sel2('text', {x: (VM.l + VW) / 2, y: VH - 2, class:'axl',
                           'text-anchor':'middle'});
  xl.textContent = 'log2FC'; svg.appendChild(xl);
}

function fmt(v){ return (v >= 0 ? '+' : '') + v.toFixed(2); }
function tip(e, id){
  const n = N[id], box = document.getElementById('tipbody');
  const stat = (n.q === null) ? '' :
      ` &nbsp;·&nbsp; p = ${n.p.toExponential(1)} &nbsp; q = ${n.q.toExponential(1)}`;
  const rows = n.g.map(s => {
    const [fc, p] = G[s];
    return `<tr><td>${s}</td><td class="v ${fc>0?'up':'dn'}">${fmt(fc)}</td>` +
           `<td class="v">${p === null ? '' : p.toExponential(1)}</td></tr>`;
  }).join('');
  drawVolcano(new Set(n.g), 8);
  box.innerHTML =
    `<h3>${n.name}</h3>` +
    `<div class="meta">${n.id} &nbsp;·&nbsp; ${n.n} genes &nbsp;·&nbsp; ` +
    `median ${fmt(n.median)} (shift ${fmt(n.shift)}) &nbsp;·&nbsp; ` +
    `${n.up} up / ${n.down} down${stat}</div>` +
    `<table><thead><tr><th>gene</th><th class="v">log2FC</th><th class="v">p</th></tr>` +
    `</thead><tbody>${rows}</tbody></table>`;
  document.getElementById('tip').scrollTop = 0;
}
function hide(){ /* panel persists until another node is hovered */ }
function placeholder(){
  drawVolcano(null, 0);
  document.getElementById('tipbody').innerHTML =
    '<div class="hint">hover a node to highlight its genes</div>';
}

['rank','topk','dim','alpha'].forEach(id =>
  document.getElementById(id).onchange = render);

(function(){                                   // drag the seam to resize the panel
  const grip = document.getElementById('grip'), tip = document.getElementById('tip');
  const MINW = 240, MAXW = () => Math.min(900, innerWidth - 320);
  function setW(w){
    w = Math.max(MINW, Math.min(MAXW(), w));
    tip.style.width = w + 'px';
    grip.style.right = w + 'px';
    document.body.style.paddingRight = (w + 22) + 'px';
    drawVolcano();                             // no args -> redraw the current highlight
  }
  grip.addEventListener('mousedown', e => {
    e.preventDefault();
    grip.classList.add('drag'); document.body.classList.add('resizing');
    const move = ev => setW(innerWidth - ev.clientX);
    const up = () => {
      grip.classList.remove('drag'); document.body.classList.remove('resizing');
      removeEventListener('mousemove', move); removeEventListener('mouseup', up);
    };
    addEventListener('mousemove', move); addEventListener('mouseup', up);
  });
  grip.addEventListener('dblclick', () => setW(__PANELW__));   // double-click to reset
  addEventListener('resize', () => drawVolcano());
})();
placeholder();
render();
</script>
"""


SYNGO_VIEWS = {
    "cellular localization": dict(aspect="C", root=_go.SYNAPSE_ROOT),
    "biological process":    dict(aspect="P", root="GO:0008150"),
    "molecular function":    dict(aspect="F", root="GO:0003674"),
}
GO_VIEWS = {
    "cellular localization": dict(aspect="C", root="GO:0005575"),
    "biological process":    dict(aspect="P", root="GO:0008150"),
    "molecular function":    dict(aspect="F", root="GO:0003674"),
}


class Explorer:
    """The built explorer. Displays itself in a notebook; `.models` holds the data.

    attributes: models ({view: dict}), html (str), path (file written, or None)
    """

    def __init__(self, html, models, path=None, height=760, width="100%"):
        self.html, self.models, self.path = html, models, path
        self.height, self.width = height, width

    @property
    def model(self):
        """The first view's model, for the single-view case."""
        return next(iter(self.models.values()))

    def _repr_html_(self):
        # srcdoc keeps the page's CSS, ids and fixed-position tooltip inside the frame,
        # so nothing leaks into the notebook and it can be displayed more than once
        doc = _html.escape(self.html, quote=True)
        return (f'<iframe srcdoc="{doc}" width="{self.width}" height="{self.height}" '
                f'style="border:1px solid #ddd;border-radius:6px;resize:vertical" '
                f'sandbox="allow-scripts"></iframe>')

    def save(self, path):
        with open(path, "w") as fh:
            fh.write(self.html)
        self.path = path
        return self

    def __repr__(self):
        bits = ", ".join(f"{v}: {len(m['nodes'])} terms" for v, m in self.models.items())
        where = f" -> {self.path}" if self.path else ""
        return f"<Explorer [{bits}]{where}>"


def explorer_html(d, out=None, source="SynGO", organism="mgi", root=None, aspect=None,
                  slim=None, views=None, min_size=3, max_size=1000, max_depth=12,
                  stats=None, center=None, title=None, subtitle=None,
                  height=760, width="100%", panel_width=360, **kw) -> Explorer:
    """Build the interactive explorer.

    views: {label: kwargs} -- builds one hierarchy per entry and adds an "organise by"
           selector to the page. Each entry's kwargs override the call's (typically
           `aspect` and `root`). Use SYNGO_VIEWS or GO_VIEWS for the standard
           localization / process / function split, or None for a single hierarchy.
           Each view gets its own breadcrumb trail and its own direction_table, computed
           automatically unless the entry supplies `stats`.

    Returns an Explorer, which renders inline in a notebook -- leave it as the last
    expression in the cell. `out` additionally writes a standalone .html file.
    """
    views = views or {"": dict(aspect=aspect, root=root)}
    models = {}
    for label, over in views.items():
        cfg = dict(source=source, organism=organism, root=root, aspect=aspect, slim=slim,
                   min_size=min_size, max_size=max_size, max_depth=max_depth,
                   center=center, **kw)
        cfg.update({k: v for k, v in over.items() if k != "stats"})
        st = over.get("stats", stats)
        if st is None:
            st = _go.direction_table(d, source=cfg["source"], organism=cfg["organism"],
                                     aspect=cfg["aspect"], slim=cfg["slim"],
                                     min_size=min_size, max_size=max_size)
        try:
            models[label or (cfg["aspect"] or "all")] = build_model(d, stats=st, **cfg)
        except ValueError as e:
            print(f"  skipped view {label!r}: {e}")
    if not models:
        raise ValueError("no view produced a hierarchy")

    title = title or f"{source or 'GO'} explorer"
    first = next(iter(models.values()))
    n_terms = sum(len(m["nodes"]) for m in models.values())
    sub = subtitle or (f"{n_terms} terms across {len(models)} view(s) · "
                       f"{len(first['genes'])} genes · click a node to drill in, "
                       f"a breadcrumb to go back")
    pts = [[g, round(float(fc), 3), (None if pv is None else float(pv))]
           for g, fc, pv in zip(d.gene, d.log2fc,
                                d.pval if "pval" in d.columns else [None] * len(d))]
    payload = {"views": models, "points": pts}
    html = (_TEMPLATE.replace("__DATA__", _json.dumps(payload, separators=(",", ":")))
                     .replace("__TITLE__", title).replace("__SUBTITLE__", sub)
                     .replace("__BASE__", f"{first['base']:+.2f}")
                     .replace("__PANELW__", str(int(panel_width))))
    ex = Explorer(html, models, height=height, width=width)
    if out:
        ex.save(out)
        print(f"wrote {out}  ({len(html) // 1024} KB, {len(models)} views)")
    return ex
