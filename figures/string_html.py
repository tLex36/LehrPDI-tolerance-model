"""Interactive STRING module network as a self-contained HTML page.

    import string_html as sh
    page = sh.module_network_html(genes, species="rat", names=NAMES,
                                  out="cnksr2_network.html")

Pan, zoom, drag proteins, lasso a cluster and ask what it is enriched for. The
enrichment runs IN THE BROWSER: a local file cannot call STRING's API (no backend, and
the browser blocks the cross-origin request), so the term-to-gene table for the
displayed proteins is embedded and the hypergeometric test is done in JavaScript. That
is faster than a round trip anyway -- a lasso returns terms immediately -- and it keeps
the page working with no network at all.

Two universes are offered for that test, because they answer different questions:
  * "displayed network" -- is this term over-represented in the selection relative to
    the other proteins on the canvas? This is what you want when comparing clusters.
  * "rat proteome"      -- the question STRING answers by default. Uses STRING's own
    background counts, so numbers are comparable to string_db.enrichment().
"""
import html as _html
import json as _json

import numpy as np
import pandas as pd

import string_db as _sdb

DEFAULT_BACKGROUND = 20000          # rat proteins with any STRING annotation, approx.
EXCLUDE_CATEGORIES = ("PMID", "NetworkNeighborAL")


def _build(genes, species="rat", score=700, resolution=1.1, min_module=4,
           min_component=3, seed=7, names=None, colors=None, **kw):
    """Graph, communities, names and layout -- the same pipeline as module_network."""
    import networkx as nx

    e = _sdb.network(genes, species=species, score=score, **kw)
    if e.empty:
        raise ValueError(f"no edges at score >= {score}")
    G = nx.Graph()
    for _, r in e.iterrows():
        G.add_edge(r.preferredName_A, r.preferredName_B, weight=float(r.score))
    keep = {n for c in nx.connected_components(G) if len(c) >= min_component for n in c}
    G = G.subgraph(keep).copy()

    comm = nx.community.louvain_communities(G, weight="weight", seed=seed,
                                            resolution=resolution)
    comm = sorted((set(c) for c in comm), key=len, reverse=True)
    big = [c for c in comm if len(c) >= min_module]
    small = [c for c in comm if len(c) < min_module]

    auto = _sdb.name_modules(big, species=species)
    labels = [(names or {}).get(i, auto.get(i, f"module {i+1}")) for i in range(len(big))]
    cls_col = [_sdb.classify_module(l) for l in labels]
    cols = [(colors or {}).get(l, c) for l, (_, c) in zip(labels, cls_col)]
    pos = _sdb._module_layout(G, big + small, seed=seed)
    return G, pos, big, small, labels, [c for c, _ in cls_col], cols


def module_network_html(genes, species="rat", score=700, out=None, title="STRING network",
                        subtitle="", names=None, colors=None, highlight=(),
                        default_labels=None, label_top=3, resolution=1.1,
                        min_module=4, min_component=3, seed=7, background=None,
                        width="100%", height=820, **kw):
    """Build the interactive page.

    default_labels: genes labelled on load. Defaults to `highlight` plus the top
            `label_top` by degree in each module; every other label is one click away.
    background: size of the proteome universe for the "rat proteome" test.
    """
    G, pos, big, small, labels, classes, cols = _build(
        genes, species=species, score=score, resolution=resolution,
        min_module=min_module, min_component=min_component, seed=seed,
        names=names, colors=colors, **kw)

    nodes = sorted(G.nodes())
    idx = {g: i for i, g in enumerate(nodes)}
    deg = dict(G.degree())
    mod_of = {g: i for i, c in enumerate(big) for g in c}

    P = np.array([pos[g] for g in nodes], float)
    P -= P.min(0)
    P /= max(P.max(), 1e-9)                       # normalise to the unit square

    hi = set(highlight) & set(G)
    if default_labels is None:
        default_labels = set(hi)
        for c in big:
            default_labels |= set(sorted(c, key=lambda g: -deg[g])[:label_top])
    default_labels = set(default_labels) & set(G)

    node_js = [dict(n=g, x=round(float(P[i, 0]), 5), y=round(float(P[i, 1]), 5),
                    m=mod_of.get(g, -1), d=deg[g],
                    h=1 if g in hi else 0, l=1 if g in default_labels else 0)
               for i, g in enumerate(nodes)]
    edge_js = [dict(s=idx[u], t=idx[v], w=round(dd["weight"], 3),
                    m=mod_of.get(u, -1) if mod_of.get(u, -1) == mod_of.get(v, -2) else -1)
               for u, v, dd in G.edges(data=True)]
    mod_js = [dict(name=labels[i], cls=classes[i], col=cols[i],
                   g=sorted(idx[x] for x in c)) for i, c in enumerate(big)]

    # term catalogue for the in-browser test: every STRING term touching >= 2 of the
    # displayed proteins, with its proteome-wide size for the second universe
    enr = _sdb.enrichment(nodes, species=species, fdr=None)
    terms = []
    if not enr.empty:
        enr = enr[~enr.category.isin(EXCLUDE_CATEGORIES)]
        for _, r in enr.iterrows():
            gs = [idx[x.strip()] for x in str(r.preferredNames).split(",")
                  if x.strip() in idx]
            if len(gs) >= 2:
                terms.append(dict(c=r.category, d=r.description,
                                  K=int(r.number_of_genes_in_background), g=sorted(gs)))

    data = dict(nodes=node_js, edges=edge_js, modules=mod_js, terms=terms,
                N=int(background or DEFAULT_BACKGROUND), score=score,
                catLabel=_sdb.CAT_LABEL)

    page = (_TPL.replace("__DATA__", _json.dumps(data, separators=(",", ":")))
                .replace("__TITLE__", _html.escape(title))
                .replace("__SUBTITLE__", subtitle or
                         f"{G.number_of_nodes()} proteins &middot; {G.number_of_edges()} "
                         f"edges &middot; STRING score &ge; {score/1000:.2f} &middot; "
                         f"{len(big)} modules"))
    if out:
        open(out, "w").write(page)
        print(f"wrote {out}  ({len(page)//1024} KB)")
    return Page(page, width, height, out)


class Page:
    def __init__(self, html, width, height, path=None):
        self.html, self.width, self.height, self.path = html, width, height, path

    def _repr_html_(self):
        doc = _html.escape(self.html, quote=True)
        return (f'<iframe srcdoc="{doc}" width="{self.width}" height="{self.height}" '
                f'style="border:1px solid #ddd;border-radius:6px;resize:vertical" '
                f'sandbox="allow-scripts"></iframe>')

    def save(self, path):
        open(path, "w").write(self.html); self.path = path; return self

    def __repr__(self):
        return f"<Page{' -> ' + self.path if self.path else ''}>"


_TPL = r"""<!doctype html>
<meta charset="utf-8">
<title>__TITLE__</title>
<style>
 :root{
   --ink:#1a1d21; --ink2:#4a5056; --mut:#878d94; --line:#e3e6e9;
   --bg:#ffffff; --panel:#f7f8f9; --accent:#0072B2;
 }
 *{box-sizing:border-box}
 body{margin:0;font:13px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
      color:var(--ink);background:var(--bg)}
 h1{font-size:17px;margin:0 0 2px}
 .sub{color:var(--mut);font-size:12px;margin-bottom:10px}
 #wrap{display:flex;gap:12px;padding:14px 16px;height:100vh}
 #left{flex:1 1 auto;position:relative;border:1px solid var(--line);border-radius:8px;
       overflow:hidden;background:#fff;min-width:0}
 #right{flex:0 0 310px;display:flex;flex-direction:column;gap:10px;overflow:hidden}
 svg{display:block;width:100%;height:100%;cursor:grab}
 svg.dragging{cursor:grabbing}
 svg.selmode{cursor:crosshair}
 #tools{position:absolute;top:10px;right:10px;display:flex;flex-direction:column;gap:6px;
        align-items:flex-end;z-index:5}
 .grp{display:flex;background:#fff;border:1px solid var(--line);border-radius:7px;
      box-shadow:0 1px 3px rgba(0,0,0,.07);overflow:hidden}
 .grp button{border:0;background:#fff;padding:6px 9px;font:12px inherit;cursor:pointer;
             color:var(--ink2);border-right:1px solid var(--line)}
 .grp button:last-child{border-right:0}
 .grp button:hover{background:#eef2f6}
 .grp button.on{background:var(--accent);color:#fff}
 .card{border:1px solid var(--line);border-radius:8px;background:var(--panel);padding:10px 11px}
 .card h2{font-size:12px;margin:0 0 7px;text-transform:uppercase;letter-spacing:.04em;color:var(--ink2)}
 #legend{overflow:auto;flex:0 0 38%;min-height:110px}
 .mrow{display:flex;gap:7px;align-items:flex-start;padding:3px 4px;border-radius:5px;cursor:pointer;font-size:12px}
 .mrow:hover{background:#eceff2}
 .mrow.off{opacity:.38}
 .sw{flex:0 0 11px;height:11px;border-radius:3px;margin-top:3px}
 #res{flex:1 1 auto;overflow:auto;min-height:0}
 table{border-collapse:collapse;width:100%;font-size:11.5px}
 th,td{text-align:left;padding:4px 5px;border-bottom:1px solid var(--line);vertical-align:top}
 th{color:var(--ink2);font-weight:600;position:sticky;top:0;background:var(--panel)}
 td.n{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
 .pill{display:inline-block;padding:1px 5px;border-radius:9px;font-size:10px;background:#e7ebef;color:var(--ink2)}
 .hint{color:var(--mut);font-size:11.5px;line-height:1.5}
 #menu{position:fixed;display:none;background:#fff;border:1px solid var(--line);border-radius:7px;
       box-shadow:0 4px 16px rgba(0,0,0,.14);z-index:50;min-width:196px;padding:4px;font-size:12.5px}
 #menu div{padding:6px 9px;border-radius:5px;cursor:pointer}
 #menu div:hover{background:#eef2f6}
 #menu .sep{height:1px;background:var(--line);margin:4px 2px;padding:0;cursor:default}
 #menu .sep:hover{background:var(--line)}
 #tip{position:fixed;display:none;background:#1a1d21;color:#fff;padding:4px 7px;border-radius:5px;
      font-size:11.5px;pointer-events:none;z-index:60}
 label.ck{display:flex;gap:6px;align-items:center;font-size:12px;color:var(--ink2);cursor:pointer}
 select{font:12px inherit;padding:3px 5px;border:1px solid var(--line);border-radius:5px;background:#fff}
 .count{color:var(--mut);font-size:11.5px}
</style>
<div id="wrap">
  <div style="display:flex;flex-direction:column;flex:1 1 auto;min-width:0">
    <h1>__TITLE__</h1>
    <div class="sub">__SUBTITLE__</div>
    <div id="left">
      <svg id="svg"></svg>
      <div id="tools">
        <div class="grp">
          <button id="t-pan" class="on" title="Pan (drag background)">✋ pan</button>
          <button id="t-box" title="Drag a rectangle to select">▭ box</button>
          <button id="t-cir" title="Drag from the centre to select a circle">◯ circle</button>
        </div>
        <div class="grp">
          <button id="zin">＋</button><button id="zout">−</button><button id="zfit">fit</button>
        </div>
        <div class="grp">
          <button id="lab">labels: default</button>
          <button id="relayout" title="Undo any proteins you have dragged">reset layout</button>
        </div>
      </div>
    </div>
  </div>
  <div id="right">
    <div class="card" id="legend"><h2>Modules</h2><div id="mods"></div></div>
    <div class="card" id="res">
      <h2>Selection</h2>
      <div id="selinfo" class="hint">
        Pick <b>box</b> or <b>circle</b>, drag across a cluster, then <b>right-click</b>
        the selection to test what it is enriched for.<br><br>
        Drag any protein to move it. Scroll to zoom. Click a protein to pin its label.
      </div>
      <div id="out"></div>
    </div>
  </div>
</div>
<div id="menu"></div>
<div id="tip"></div>
<script>
const D = __DATA__;
const SVG = document.getElementById("svg"), NS = "http://www.w3.org/2000/svg";
const el = (t,a={}) => { const e=document.createElementNS(NS,t);
  for(const k in a) e.setAttribute(k,a[k]); return e; };
const $ = id => document.getElementById(id);

/* ---------- state ---------- */
const N = D.nodes.length;
const home = D.nodes.map(n => [n.x, n.y]);          // layout as built, for reset
let   P    = home.map(p => p.slice());              // live positions, 0..1
let   view = {k:1, x:0, y:0};
let   tool = "pan";
let   labelMode = "default";                        // default | all | none
const pinned = new Set(D.nodes.map((n,i)=>n.l?i:-1).filter(i=>i>=0));
let   sel = new Set();
const modOff = new Set();

/* ---------- geometry ---------- */
let W = 100, H = 100, PAD = 34;
function fit(){
  const r = SVG.getBoundingClientRect();
  W = Math.max(r.width, 50); H = Math.max(r.height, 50);
}
const sx = u => PAD + u * (W - 2*PAD);
const sy = u => PAD + u * (H - 2*PAD);
const ux = p => (p - PAD) / (W - 2*PAD);
const uy = p => (p - PAD) / (H - 2*PAD);
/* screen position of node i, including the zoom/pan transform */
const SX = i => sx(P[i][0]) * view.k + view.x;
const SY = i => sy(P[i][1]) * view.k + view.y;

/* ---------- layers ---------- */
const gRoot = el("g"), gHull = el("g"), gEdge = el("g"), gNode = el("g"),
      gLead = el("g"), gLab = el("g"), gSel = el("g");
gRoot.appendChild(gHull); gRoot.appendChild(gEdge); gRoot.appendChild(gNode);
gRoot.appendChild(gLead); gRoot.appendChild(gLab);
SVG.appendChild(gRoot); SVG.appendChild(gSel);

/* ---------- convex hull, padded + smoothed (mirrors the matplotlib figure) ---- */
function hull(pts){
  if(pts.length < 3) return pts.slice();
  const p = pts.slice().sort((a,b)=> a[0]-b[0] || a[1]-b[1]);
  const cross=(o,a,b)=>(a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0]);
  const lo=[]; for(const q of p){ while(lo.length>=2 && cross(lo[lo.length-2],lo[lo.length-1],q)<=0) lo.pop(); lo.push(q); }
  const up=[]; for(let i=p.length-1;i>=0;i--){ const q=p[i];
    while(up.length>=2 && cross(up[up.length-2],up[up.length-1],q)<=0) up.pop(); up.push(q); }
  lo.pop(); up.pop(); return lo.concat(up);
}
function blob(pts, pad){
  if(!pts.length) return "";
  const cx = pts.reduce((s,p)=>s+p[0],0)/pts.length,
        cy = pts.reduce((s,p)=>s+p[1],0)/pts.length;
  let V;
  if(pts.length <= 5){                       // small modules: a circle, not a triangle
    let r = 0; for(const p of pts) r = Math.max(r, Math.hypot(p[0]-cx,p[1]-cy));
    r += pad; V = [];
    for(let a=0;a<12;a++) V.push([cx+r*Math.cos(a/12*2*Math.PI), cy+r*Math.sin(a/12*2*Math.PI)]);
  } else {
    V = hull(pts).map(p=>{
      const d = Math.hypot(p[0]-cx,p[1]-cy) || 1;
      return [p[0]+(p[0]-cx)/d*pad, p[1]+(p[1]-cy)/d*pad];
    });
  }
  if(V.length < 3) return "";
  let s = `M ${V[0][0].toFixed(1)} ${V[0][1].toFixed(1)}`;   // Catmull-Rom -> bezier
  const n = V.length;
  for(let i=0;i<n;i++){
    const p0=V[(i-1+n)%n], p1=V[i], p2=V[(i+1)%n], p3=V[(i+2)%n];
    const c1=[p1[0]+(p2[0]-p0[0])/6, p1[1]+(p2[1]-p0[1])/6];
    const c2=[p2[0]-(p3[0]-p1[0])/6, p2[1]-(p3[1]-p1[1])/6];
    s += ` C ${c1[0].toFixed(1)} ${c1[1].toFixed(1)}, ${c2[0].toFixed(1)} ${c2[1].toFixed(1)}, ${p2[0].toFixed(1)} ${p2[1].toFixed(1)}`;
  }
  return s + " Z";
}

/* ---------- draw ---------- */
const nodeEls = [], labEls = [], hullEls = [], edgeEls = [];
function visible(i){ const m = D.nodes[i].m; return m < 0 || !modOff.has(m); }
const rOf = i => (3.1 + 0.46*Math.min(D.nodes[i].d, 16)) * Math.sqrt(view.k);

function build(){
  gHull.textContent=""; gEdge.textContent=""; gNode.textContent="";
  gLab.textContent=""; gLead.textContent="";
  hullEls.length=0; edgeEls.length=0; nodeEls.length=0; labEls.length=0;

  D.modules.forEach((m,mi)=>{
    const p = el("path", {fill:m.col, "fill-opacity":.10, stroke:m.col,
                          "stroke-width":1.3, "stroke-linejoin":"round"});
    gHull.appendChild(p); hullEls.push(p);
  });
  D.edges.forEach(e=>{
    const col = e.m >= 0 ? D.modules[e.m].col : "#c8ccd0";
    const l = el("line",{stroke:col, "stroke-opacity": e.m>=0?.5:.3,
                         "stroke-width": e.m>=0?1.0:.6, "stroke-linecap":"round"});
    gEdge.appendChild(l); edgeEls.push(l);
  });
  D.nodes.forEach((n,i)=>{
    const col = n.m>=0 ? D.modules[n.m].col : "#9aa0a6";
    const c = el("circle",{fill:col, stroke:"#fff", "stroke-width":1.1, cursor:"pointer"});
    c.dataset.i = i; gNode.appendChild(c); nodeEls.push(c);
    const t = el("text",{"font-size":10.5,"text-anchor":"middle","pointer-events":"none",
                         fill: n.h? "#111":"#333", "font-weight": n.h? 700:400,
                         stroke:"#fff","stroke-width":2.6,"paint-order":"stroke"});
    t.textContent = n.n; gLab.appendChild(t); labEls.push(t);
  });
  draw();
}

function labelled(i){
  if(labelMode === "all")  return true;
  if(labelMode === "none") return sel.has(i);
  return pinned.has(i) || sel.has(i);
}

function draw(){
  D.modules.forEach((m,mi)=>{
    const pts = m.g.filter(visible).map(i=>[SX(i),SY(i)]);
    const d = modOff.has(mi) ? "" : blob(pts, 17*Math.sqrt(view.k));
    hullEls[mi].setAttribute("d", d);
    hullEls[mi].style.display = d ? "" : "none";
  });
  D.edges.forEach((e,k)=>{
    const L = edgeEls[k];
    if(!visible(e.s) || !visible(e.t)){ L.style.display="none"; return; }
    L.style.display="";
    L.setAttribute("x1",SX(e.s)); L.setAttribute("y1",SY(e.s));
    L.setAttribute("x2",SX(e.t)); L.setAttribute("y2",SY(e.t));
  });
  const fs = 10.5*Math.sqrt(view.k), shown = [];
  D.nodes.forEach((n,i)=>{
    const c = nodeEls[i];
    if(!visible(i)){ c.style.display="none"; labEls[i].style.display="none"; return; }
    c.style.display="";
    c.setAttribute("cx",SX(i)); c.setAttribute("cy",SY(i)); c.setAttribute("r",rOf(i));
    const on = sel.has(i);
    c.setAttribute("stroke", on ? "#111" : "#fff");
    c.setAttribute("stroke-width", on ? 2.2 : 1.1);
    const t = labEls[i];
    if(labelled(i)){
      t.style.display="";
      t.setAttribute("font-size", fs.toFixed(1));
      shown.push(i);
    } else t.style.display="none";
  });
  placeLabels(shown, fs);
}

/* Labels are separated in screen space, exactly as the static figure does it: a
   dense module holds more names than fit above its nodes, so they have to spread and
   take a leader line with them. Skipped mid-drag to keep dragging smooth. */
function placeLabels(ids, fs){
  gLead.textContent = "";
  if(!ids.length) return;
  const pos = ids.map(i=>[SX(i), SY(i)-rOf(i)-5]);
  if(!drag){
    const w = ids.map(i=>D.nodes[i].n.length*fs*0.62 + 5),
          h = fs*1.35 + 4;
    for(let it=0; it<90; it++){
      let moved = false;
      for(let a=0;a<ids.length;a++) for(let b=a+1;b<ids.length;b++){
        const dx = pos[b][0]-pos[a][0], dy = pos[b][1]-pos[a][1];
        const ox = (w[a]+w[b])/2 - Math.abs(dx), oy = h - Math.abs(dy);
        if(ox>0 && oy>0){
          moved = true;
          if(oy*2.2 < ox){ const s=(oy/2+.5)*(dy>=0?1:-1); pos[a][1]-=s; pos[b][1]+=s; }
          else           { const s=(ox/2+.5)*(dx>=0?1:-1); pos[a][0]-=s; pos[b][0]+=s; }
        }
      }
      if(!moved) break;
    }
  }
  ids.forEach((i,k)=>{
    const t = labEls[i];
    t.setAttribute("x", pos[k][0].toFixed(1));
    t.setAttribute("y", pos[k][1].toFixed(1));
    const nx = SX(i), ny = SY(i);
    if(Math.hypot(pos[k][0]-nx, pos[k][1]-ny) > rOf(i)+14){
      const L = el("line",{x1:nx, y1:ny, x2:pos[k][0], y2:pos[k][1]+fs*0.32,
                           stroke:"#b9bec4", "stroke-width":.7});
      gLead.appendChild(L);
    }
  });
}

/* ---------- zoom & pan ---------- */
function zoomAt(cx, cy, f){
  const k2 = Math.min(9, Math.max(.25, view.k*f));
  f = k2/view.k;
  view.x = cx - (cx-view.x)*f; view.y = cy - (cy-view.y)*f; view.k = k2;
  draw();
}
SVG.addEventListener("wheel", ev=>{
  ev.preventDefault();
  const r = SVG.getBoundingClientRect();
  zoomAt(ev.clientX-r.left, ev.clientY-r.top, ev.deltaY<0 ? 1.13 : 1/1.13);
},{passive:false});
$("zin").onclick  = ()=> zoomAt(W/2,H/2,1.25);
$("zout").onclick = ()=> zoomAt(W/2,H/2,1/1.25);
$("zfit").onclick = ()=>{ view={k:1,x:0,y:0}; draw(); };
$("relayout").onclick = ()=>{ P = home.map(p=>p.slice()); draw(); };

/* ---------- tools ---------- */
function setTool(t){
  tool=t;
  for(const [id,v] of [["t-pan","pan"],["t-box","box"],["t-cir","circle"]])
    $(id).classList.toggle("on", v===t);
  SVG.classList.toggle("selmode", t!=="pan");
}
$("t-pan").onclick=()=>setTool("pan");
$("t-box").onclick=()=>setTool("box");
$("t-cir").onclick=()=>setTool("circle");
$("lab").onclick = ()=>{
  labelMode = labelMode==="default" ? "all" : labelMode==="all" ? "none" : "default";
  $("lab").textContent = "labels: "+labelMode;
  draw();
};

/* ---------- interaction: drag a node, pan, or draw a selection ---------- */
let drag = null, shape = null;
SVG.addEventListener("mousedown", ev=>{
  if(ev.button === 2) return;                       // right-click handled separately
  hideMenu();
  const r = SVG.getBoundingClientRect(), mx = ev.clientX-r.left, my = ev.clientY-r.top;
  const i = ev.target.dataset && ev.target.dataset.i !== undefined
          ? +ev.target.dataset.i : -1;
  if(i >= 0){                                        // moving a protein
    drag = {i, moved:false, ox:mx-SX(i), oy:my-SY(i)};
    return;
  }
  if(tool === "pan"){ drag = {pan:true, x:mx-view.x, y:my-view.y}; SVG.classList.add("dragging"); }
  else {
    shape = {x0:mx, y0:my, x1:mx, y1:my};
    gSel.textContent="";
    const s = tool==="box" ? el("rect",{}) : el("circle",{});
    s.setAttribute("fill","rgba(0,114,178,.10)");
    s.setAttribute("stroke","#0072B2"); s.setAttribute("stroke-width",1.2);
    s.setAttribute("stroke-dasharray","4 3");
    gSel.appendChild(s); shape.el = s;
  }
});
window.addEventListener("mousemove", ev=>{
  const r = SVG.getBoundingClientRect(), mx = ev.clientX-r.left, my = ev.clientY-r.top;
  if(drag && drag.pan){ view.x = mx-drag.x; view.y = my-drag.y; draw(); return; }
  if(drag){
    drag.moved = true;
    P[drag.i] = [ux((mx-drag.ox-view.x)/view.k), uy((my-drag.oy-view.y)/view.k)];
    draw(); return;
  }
  if(shape){
    shape.x1=mx; shape.y1=my;
    if(tool==="box"){
      shape.el.setAttribute("x",Math.min(shape.x0,mx)); shape.el.setAttribute("y",Math.min(shape.y0,my));
      shape.el.setAttribute("width",Math.abs(mx-shape.x0)); shape.el.setAttribute("height",Math.abs(my-shape.y0));
    } else {
      shape.el.setAttribute("cx",shape.x0); shape.el.setAttribute("cy",shape.y0);
      shape.el.setAttribute("r",Math.hypot(mx-shape.x0,my-shape.y0));
    }
  }
});
window.addEventListener("mouseup", ev=>{
  SVG.classList.remove("dragging");
  if(drag && !drag.pan && !drag.moved){              // a click, not a drag: pin label
    const i = drag.i;
    if(pinned.has(i)) pinned.delete(i); else pinned.add(i);
    draw();
  }
  drag = null;
  if(shape){
    const picked = new Set();
    for(let i=0;i<N;i++){
      if(!visible(i)) continue;
      const x=SX(i), y=SY(i);
      if(tool==="box"){
        if(x>=Math.min(shape.x0,shape.x1) && x<=Math.max(shape.x0,shape.x1) &&
           y>=Math.min(shape.y0,shape.y1) && y<=Math.max(shape.y0,shape.y1)) picked.add(i);
      } else {
        if(Math.hypot(x-shape.x0,y-shape.y0) <= Math.hypot(shape.x1-shape.x0,shape.y1-shape.y0))
          picked.add(i);
      }
    }
    if(ev.shiftKey) picked.forEach(i=>sel.add(i)); else sel = picked;
    gSel.textContent=""; shape=null;
    showSelection(); draw();
  }
});

/* ---------- hover ---------- */
SVG.addEventListener("mousemove", ev=>{
  const i = ev.target.dataset && ev.target.dataset.i !== undefined ? +ev.target.dataset.i : -1;
  const tip = $("tip");
  if(i < 0 || drag || shape){ tip.style.display="none"; return; }
  const n = D.nodes[i];
  tip.innerHTML = `<b>${n.n}</b> &middot; ${n.d} edges` +
                  (n.m>=0 ? `<br>${D.modules[n.m].name}` : "");
  tip.style.display="block";
  tip.style.left = (ev.clientX+13)+"px"; tip.style.top = (ev.clientY+13)+"px";
});
SVG.addEventListener("mouseleave", ()=> $("tip").style.display="none");

/* ---------- context menu ---------- */
const menu = $("menu");
function hideMenu(){ menu.style.display="none"; }
window.addEventListener("click", hideMenu);
SVG.addEventListener("contextmenu", ev=>{
  ev.preventDefault();
  const i = ev.target.dataset && ev.target.dataset.i !== undefined ? +ev.target.dataset.i : -1;
  if(i>=0 && !sel.has(i)){ sel = new Set([i]); showSelection(); draw(); }
  const items = [];
  if(sel.size){
    items.push(["Find enriched STRING terms", ()=>runEnrich()]);
    items.push(["Select whole module of selection", ()=>{
      const ms = new Set([...sel].map(j=>D.nodes[j].m).filter(m=>m>=0));
      ms.forEach(m => D.modules[m].g.forEach(j=>sel.add(j)));
      showSelection(); draw();
    }]);
    items.push(["Grow selection by one edge", ()=>{
      const add=new Set();
      D.edges.forEach(e=>{ if(sel.has(e.s)) add.add(e.t); if(sel.has(e.t)) add.add(e.s); });
      add.forEach(j=>sel.add(j)); showSelection(); draw();
    }]);
    items.push(["sep"]);
    items.push(["Pin labels for selection", ()=>{ sel.forEach(j=>pinned.add(j)); draw(); }]);
    items.push(["Unpin labels for selection", ()=>{ sel.forEach(j=>pinned.delete(j)); draw(); }]);
    items.push(["sep"]);
    items.push(["Clear selection", ()=>{ sel=new Set(); showSelection(); draw(); }]);
  } else {
    items.push(["Select all", ()=>{ sel=new Set(D.nodes.map((_,j)=>j).filter(visible)); showSelection(); draw(); }]);
  }
  menu.textContent="";
  for(const it of items){
    if(it[0]==="sep"){ const s=document.createElement("div"); s.className="sep"; menu.appendChild(s); continue; }
    const d=document.createElement("div"); d.textContent=it[0];
    d.onclick=e=>{ e.stopPropagation(); hideMenu(); it[1](); };
    menu.appendChild(d);
  }
  menu.style.display="block";
  menu.style.left=Math.min(ev.clientX, innerWidth-215)+"px";
  menu.style.top =Math.min(ev.clientY, innerHeight-230)+"px";
});

/* ---------- hypergeometric enrichment, in the browser ---------- */
const LG=[0]; (function(){ let s=0; for(let i=1;i<=60000;i++){ s+=Math.log(i); LG.push(s); } })();
const lnFact = n => n<LG.length ? LG[n] : (n+.5)*Math.log(n)-n+0.9189385332046727+1/(12*n);
const lnC = (n,k) => (k<0||k>n) ? -Infinity : lnFact(n)-lnFact(k)-lnFact(n-k);
function hyperSF(k,N,K,n){            /* P(X >= k) */
  if(k<=0) return 1;
  const hi=Math.min(n,K); if(k>hi) return 0;
  const den=lnC(N,n); let s=0;
  for(let i=k;i<=hi;i++){ const v=lnC(K,i)+lnC(N-K,n-i)-den; if(v>-745) s+=Math.exp(v); }
  return Math.min(1,s);
}
function bh(ps){
  const m=ps.length, o=ps.map((p,i)=>[p,i]).sort((a,b)=>a[0]-b[0]);
  const q=new Array(m); let prev=1;
  for(let r=m-1;r>=0;r--){ prev=Math.min(prev, o[r][0]*m/(r+1)); q[o[r][1]]=prev; }
  return q;
}
function runEnrich(){
  const universe = $("uni") ? $("uni").value : "net";
  const S = [...sel];
  const dest = $("out2") || $("out");
  if(S.length < 2){ dest.innerHTML='<div class="hint">Select at least two proteins.</div>'; return; }
  const inSel = new Array(N).fill(false); S.forEach(i=>inSel[i]=true);
  const n = S.length, Nbg = universe==="net" ? N : D.N;
  const rows=[];
  for(const t of D.terms){
    let k=0; for(const g of t.g) if(inSel[g]) k++;
    if(k<2) continue;
    const K = universe==="net" ? t.g.length : t.K;
    if(K<k) continue;
    const p = hyperSF(k,Nbg,K,n);
    rows.push({c:t.c,d:t.d,k,K,p, fold:(k/n)/(K/Nbg),
               genes:t.g.filter(g=>inSel[g]).map(g=>D.nodes[g].n)});
  }
  if(!rows.length){ dest.innerHTML='<div class="hint">No term covers 2+ of the selected proteins.</div>'; return; }
  const q = bh(rows.map(r=>r.p));
  rows.forEach((r,i)=>r.q=q[i]);
  rows.sort((a,b)=>a.p-b.p || b.k-a.k);
  const sig = rows.filter(r=>r.q<0.05);
  const show = (sig.length?sig:rows).slice(0,40);
  const fmt = v => v<1e-4 ? v.toExponential(1) : v.toFixed(4);
  let h = `<div class="count" style="margin-bottom:6px">${sig.length} terms at FDR&lt;0.05`
        + ` of ${rows.length} tested &middot; universe: `
        + `${universe==="net" ? "displayed network ("+N+")" : "rat proteome ("+D.N+")"}`
        + (sig.length? "" : " &mdash; showing top by p") + `</div>`;
  h += "<table><thead><tr><th>term</th><th class='n'>hits</th><th class='n'>fold</th>"
     + "<th class='n'>q</th></tr></thead><tbody>";
  for(const r of show){
    const cat = D.catLabel[r.c] || r.c;
    h += `<tr><td><b>${esc(r.d)}</b><br><span class="pill">${esc(cat)}</span> `
       + `<span class="count">${esc(r.genes.slice(0,9).join(", "))}`
       + `${r.genes.length>9?" …":""}</span></td>`
       + `<td class="n">${r.k}/${r.K}</td><td class="n">${r.fold.toFixed(1)}×</td>`
       + `<td class="n">${fmt(r.q)}</td></tr>`;
  }
  ($("out2") || $("out")).innerHTML = h + "</tbody></table>";
}
const esc = s => String(s).replace(/[&<>"]/g, c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));

function showSelection(){
  const S=[...sel];
  if(!S.length){
    $("selinfo").style.display="";
    $("out").innerHTML=""; return;
  }
  $("selinfo").style.display="none";
  const names=S.map(i=>D.nodes[i].n).sort();
  $("out").innerHTML =
    `<div style="margin-bottom:7px">
       <b>${S.length} selected</b>
       <div class="count" style="margin:4px 0 8px">${esc(names.slice(0,28).join(", "))}${names.length>28?" …":""}</div>
       <label class="ck" style="margin-bottom:7px">universe
         <select id="uni"><option value="net">displayed network</option>
                          <option value="prot">rat proteome</option></select></label>
       <button id="go" style="font:12px inherit;padding:5px 10px;border:1px solid var(--line);
         border-radius:6px;background:#fff;cursor:pointer">Find enriched terms</button>
     </div><div id="out2"></div>`;
  $("go").onclick = runEnrich;
  $("uni").onchange = runEnrich;
}

/* ---------- legend ---------- */
function legend(){
  const box=$("mods"); box.textContent="";
  D.modules.forEach((m,i)=>{
    const row=document.createElement("div"); row.className="mrow";
    row.innerHTML=`<span class="sw" style="background:${m.col}"></span>
      <span><b>${i+1}.</b> ${esc(m.name)} <span class="count">(${m.g.length})</span></span>`;
    row.onclick=()=>{ if(modOff.has(i)) modOff.delete(i); else modOff.add(i);
                      row.classList.toggle("off"); draw(); };
    row.oncontextmenu=e=>{ e.preventDefault(); e.stopPropagation();
                           sel=new Set(m.g); showSelection(); draw(); };
    box.appendChild(row);
  });
}

/* ---------- go ---------- */
function boot(){ fit(); build(); legend(); }
window.addEventListener("resize", ()=>{ fit(); draw(); });
boot();
</script>
"""


# =============================================================== enrichment dot plot

def enrichment_html(sets, out=None, title="STRING enrichment", subtitle="",
                    categories=("KEGG", "Process", "Component", "RCTM", "Function"),
                    top=8, max_fdr=0.05, max_term_size=600, x="fdr",
                    width="100%", height=900):
    """Interactive dot plot of STRING enrichment: hover a dot for its genes.

    sets: {name: enrichment frame} -- one tab per entry, so two baits can be compared
          without scrolling past one another.
    top:  terms kept per category, ranked by FDR.
    x:    "fdr" plots -log10(FDR); "fold" plots fold enrichment.

    A dot plot rather than bars: the gene count becomes the dot's area instead of
    needing a number printed at the end of every bar, which leaves the hover layer free
    to carry the thing a static figure cannot show at all -- which genes are in the term.
    """
    if isinstance(sets, pd.DataFrame):
        sets = {"": sets}

    packed = []
    for name, e in sets.items():
        if e is None or e.empty:
            packed.append(dict(name=name, terms=[])); continue
        x_ = e.copy()
        if max_term_size is not None:
            x_ = x_[x_.number_of_genes_in_background <= max_term_size]
        x_ = x_[x_.fdr < max_fdr]
        x_ = x_[x_.category.isin(categories)]
        keep = []
        for c in categories:
            sub = x_[x_.category == c].nsmallest(top, "fdr")
            for _, r in sub.iterrows():
                genes = [g.strip() for g in str(r.preferredNames).split(",") if g.strip()]
                keep.append(dict(
                    d=str(r.description), c=str(r.category),
                    k=int(r.number_of_genes), K=int(r.number_of_genes_in_background),
                    f=round(float(r.fold_enrichment), 2),
                    q=float(r.fdr), g=genes))
        packed.append(dict(name=str(name), terms=keep))

    data = dict(sets=packed, catLabel=_sdb.CAT_LABEL, catColor=_sdb.CAT_COLORS,
                order=list(categories), x=x)
    page = (_ETPL.replace("__DATA__", _json.dumps(data, separators=(",", ":")))
                 .replace("__TITLE__", _html.escape(title))
                 .replace("__SUBTITLE__", subtitle))
    if out:
        open(out, "w").write(page)
        print(f"wrote {out}  ({len(page)//1024} KB)")
    return Page(page, width, height, out)


_ETPL = r"""<!doctype html>
<meta charset="utf-8">
<title>__TITLE__</title>
<style>
 :root{--ink:#1a1d21;--ink2:#4a5056;--mut:#878d94;--line:#e3e6e9;--accent:#0072B2}
 *{box-sizing:border-box}
 body{margin:0;padding:18px 22px 40px;color:var(--ink);background:#fff;
      font:13px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif}
 h1{font-size:18px;margin:0 0 2px}
 .sub{color:var(--mut);font-size:12px;margin-bottom:12px}
 .bar{display:flex;gap:14px;align-items:center;flex-wrap:wrap;margin:0 0 6px}
 .tabs{display:flex;gap:3px;border-bottom:2px solid var(--line);flex:1 1 auto}
 .tab{border:0;background:none;padding:8px 16px;font:600 13px inherit;color:var(--mut);
      cursor:pointer;border-bottom:2px solid transparent;margin-bottom:-2px}
 .tab:hover{color:var(--ink2)} .tab.on{color:var(--accent);border-bottom-color:var(--accent)}
 .grp{display:flex;border:1px solid var(--line);border-radius:7px;overflow:hidden}
 .grp button{border:0;background:#fff;padding:5px 10px;font:12px inherit;cursor:pointer;
             color:var(--ink2);border-right:1px solid var(--line)}
 .grp button:last-child{border-right:0}
 .grp button:hover{background:#eef2f6} .grp button.on{background:var(--accent);color:#fff}
 .hint{color:var(--mut);font-size:11.5px}
 #tip{position:fixed;display:none;max-width:430px;background:#1a1d21;color:#fff;
      padding:8px 10px;border-radius:7px;font-size:11.5px;line-height:1.5;z-index:60;
      pointer-events:none;box-shadow:0 4px 18px rgba(0,0,0,.25)}
 #tip b{font-size:12.5px} #tip .gl{color:#cdd3d9;margin-top:5px;word-break:break-word}
 #tip .meta{color:#9aa4ad;margin-top:3px}
 svg{display:block;width:100%;height:auto}
 .lab{font-size:11px;fill:#2b2f34} .cat{font-size:10.5px;font-weight:700}
 .ax{stroke:#c9ccd1;stroke-width:1} .gl{stroke:#eef0f2;stroke-width:1}
 .axl{font-size:10.5px;fill:#878d94}
 .pin{font-size:10px;fill:#4a5056}
</style>
<h1>__TITLE__</h1>
<div class="sub">__SUBTITLE__</div>
<div class="bar">
  <div class="tabs" id="tabs"></div>
  <div class="grp">
    <button id="x-fdr" class="on">−log10 FDR</button><button id="x-fold">fold</button>
  </div>
  <div class="grp"><button id="labs">labels: none</button></div>
</div>
<div class="hint">hover a dot for its genes &middot; click a dot to pin its gene list &middot;
  dot area = genes in term</div>
<div id="plot"></div>
<div id="tip"></div>
<script>
const D = __DATA__, NS="http://www.w3.org/2000/svg";
const el=(t,a={})=>{const e=document.createElementNS(NS,t);for(const k in a)e.setAttribute(k,a[k]);return e;};
const $=id=>document.getElementById(id);
const esc=s=>String(s).replace(/[&<>"]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
let cur=0, xmode=D.x==="fold"?"fold":"fdr", labMode="none";
const pinned = D.sets.map(()=>new Set());

const nl = q => -Math.log10(Math.max(q,1e-300));
const val = t => xmode==="fold" ? t.f : nl(t.q);

function tabs(){
  const b=$("tabs"); b.textContent="";
  if(D.sets.length<2 && !D.sets[0].name) return;
  D.sets.forEach((s,i)=>{
    const t=document.createElement("button");
    t.className="tab"+(i===cur?" on":""); t.textContent=s.name||("set "+(i+1));
    t.onclick=()=>{cur=i;tabs();draw();}; b.appendChild(t);
  });
}

function draw(){
  const S=D.sets[cur], T=S.terms;
  const host=$("plot"); host.textContent="";
  if(!T.length){ host.innerHTML='<div class="hint" style="padding:22px 0">nothing enriched at this threshold</div>'; return; }

  // rows grouped by category, in the order the caller asked for
  const rows=[]; let y=0;
  for(const c of D.order){
    const sub=T.filter(t=>t.c===c);
    if(!sub.length) continue;
    sub.sort((a,b)=>val(b)-val(a));
    rows.push({head:c, y:y}); y+=1;
    for(const t of sub){ rows.push({t, y:y}); y+=1; }
    y+=.5;
  }
  const RH=19, L=330, R=230, TOP=16, BOT=42;
  const W=1180, H=TOP+BOT+y*RH;
  const svg=el("svg",{viewBox:`0 0 ${W} ${H}`}); host.appendChild(svg);

  const vals=T.map(val), vmax=Math.max(...vals), vmin=0;
  const X=v=>L+(v-vmin)/((vmax-vmin)||1)*(W-L-R);
  const Y=r=>TOP+r*RH+RH/2;
  const kmax=Math.max(...T.map(t=>t.k));
  const rad=k=>3.2+7.5*Math.sqrt(k/kmax);

  // gridlines + axis
  const ticks=4;
  for(let i=0;i<=ticks;i++){
    const v=vmin+(vmax-vmin)*i/ticks, px=X(v);
    svg.appendChild(el("line",{class:"gl",x1:px,x2:px,y1:TOP,y2:H-BOT}));
    const tx=el("text",{x:px,y:H-BOT+15,class:"axl","text-anchor":"middle"});
    tx.textContent = xmode==="fold" ? v.toFixed(1)+"×" : v.toFixed(0);
    svg.appendChild(tx);
  }
  svg.appendChild(el("line",{class:"ax",x1:L,x2:W-R,y1:H-BOT,y2:H-BOT}));
  const xl=el("text",{x:(L+W-R)/2,y:H-BOT+32,class:"axl","text-anchor":"middle"});
  xl.textContent = xmode==="fold" ? "fold enrichment" : "−log₁₀ FDR";
  svg.appendChild(xl);

  for(const r of rows){
    if(r.head!==undefined){
      const c=el("text",{x:8,y:Y(r.y)+3,class:"cat",fill:D.catColor[r.head]||"#666"});
      c.textContent=(D.catLabel[r.head]||r.head).toUpperCase(); svg.appendChild(c);
      continue;
    }
    const t=r.t, yy=Y(r.y), col=D.catColor[t.c]||"#666";
    const lab=el("text",{x:L-12,y:yy+3.5,class:"lab","text-anchor":"end"});
    lab.textContent = t.d.length>52 ? t.d.slice(0,51)+"…" : t.d;
    svg.appendChild(lab);
    svg.appendChild(el("line",{x1:L,x2:X(val(t)),y1:yy,y2:yy,stroke:col,
                               "stroke-opacity":.28,"stroke-width":1.4}));
    const c=el("circle",{cx:X(val(t)),cy:yy,r:rad(t.k),fill:col,
                         stroke:"#fff","stroke-width":1.2,cursor:"pointer"});
    c.addEventListener("mousemove",ev=>{
      const tip=$("tip");
      tip.innerHTML=`<b>${esc(t.d)}</b>`
        +`<div class="meta">${esc(D.catLabel[t.c]||t.c)} · ${t.k} of ${t.K} proteins · `
        +`${t.f}× · FDR ${t.q.toExponential(1)}</div>`
        +`<div class="gl">${esc(t.g.join(", "))}</div>`;
      tip.style.display="block";
      const w=tip.offsetWidth, h=tip.offsetHeight;
      tip.style.left=Math.min(ev.clientX+14, innerWidth-w-10)+"px";
      tip.style.top =Math.min(ev.clientY+14, innerHeight-h-10)+"px";
    });
    c.addEventListener("mouseleave",()=>$("tip").style.display="none");
    c.addEventListener("click",()=>{
      const P=pinned[cur]; if(P.has(t.d)) P.delete(t.d); else P.add(t.d);
      $("tip").style.display="none"; draw();
    });
    svg.appendChild(c);

    const show = labMode==="genes" || (labMode==="pinned" && pinned[cur].has(t.d))
                 || pinned[cur].has(t.d);
    if(labMode==="counts" && !pinned[cur].has(t.d)){
      const n=el("text",{x:X(val(t))+rad(t.k)+5,y:yy+3.5,class:"pin"});
      n.textContent=`${t.k}/${t.K}`; svg.appendChild(n);
    } else if(show){
      const avail=W-R-X(val(t))-rad(t.k)-10, per=5.4;
      let s=t.g.join(", ");
      const max=Math.max(8, Math.floor((avail+R-18)/per));
      if(s.length>max) s=s.slice(0,max-1)+"…";
      const n=el("text",{x:X(val(t))+rad(t.k)+6,y:yy+3.5,class:"pin"});
      n.textContent=s; svg.appendChild(n);
    }
  }
}

$("x-fdr").onclick =()=>{xmode="fdr"; $("x-fdr").classList.add("on"); $("x-fold").classList.remove("on"); draw();};
$("x-fold").onclick=()=>{xmode="fold";$("x-fold").classList.add("on");$("x-fdr").classList.remove("on"); draw();};
$("labs").onclick=()=>{
  labMode = labMode==="none" ? "counts" : labMode==="counts" ? "genes" : "none";
  $("labs").textContent="labels: "+labMode; draw();
};
tabs(); draw();
</script>
"""
