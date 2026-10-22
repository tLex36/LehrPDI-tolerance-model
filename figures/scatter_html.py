"""Interactive gene scatter: one dot per gene, hover for its numbers, click to label.

    import scatter_html as sc
    sc.gene_scatter_html(d, x="2D150_log2", y="DelESEV3_log2", gene="gene",
                         cls="cls", out="scatter.html")

Built for two-condition comparisons where both axes are enrichments against a shared
control, so the diagonal is the null and distance from it is the effect. Pan, zoom,
hover any dot for the fields you name, click to pin its label, search for a gene by
name, and click the legend to mute a class.
"""
import html as _html
import json as _json

import numpy as np
import pandas as pd

CLASS_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7",
                "#E69F00", "#56B4E9", "#8c6bb1", "#B07A3C"]
MUTED = "#c9ccd1"


def gene_scatter_html(d, x, y, gene="gene", cls=None, out=None,
                      title="Gene scatter", subtitle="",
                      xlabel=None, ylabel=None, hover=(), hover_lower=None,
                      colors=None,
                      default_labels=(), rings=(), cap=None, cap_note="",
                      diagonal=True, order=None, muted=(), xlim=None, ylim=None,
                      lower_bound=None,
                      width="100%", height=860):
    """d: frame with one row per gene.

    hover:  columns to show in the tooltip, as names or (label, column) pairs.
    hover_lower: a different field list for lower-bound points. Their control-relative
            ratios are a software ceiling rather than a measurement, so showing the raw
            per-bait abundances instead tells the reader something true.
    cls:    column of class labels; each gets a colour and a legend entry you can click
            to mute. `muted` names classes drawn greyed out from the start.
    cap:    axis value where points were clamped (PD's ratio ceiling). Draws the guide
            lines and jitters the pile-up so the stack is countable.
    rings:  genes to outline, for a set you want findable without labelling everything.
    lower_bound: boolean column. True marks a point whose position is a LOWER BOUND
            rather than a measurement -- drawn hollow, so a reader cannot mistake an
            estimate for a number. Use it wherever a ratio hit a software ceiling.
    xlim/ylim: default view, as (lo, hi). Points outside are still drawn and reachable
            by zooming out -- this only decides where the view starts, so a handful of
            extreme values cannot leave most of the canvas empty.
    """
    f = d.copy()
    f = f[np.isfinite(f[x]) & np.isfinite(f[y])]

    norm = [(h, h) if isinstance(h, str) else tuple(h) for h in hover]
    lowr = norm if hover_lower is None else [
        (h, h) if isinstance(h, str) else tuple(h) for h in hover_lower]
    norm = [(l, c) for l, c in norm if c in f.columns]
    lowr = [(l, c) for l, c in lowr if c in f.columns]
    # one value array per point covering the union, plus the row indices each kind of
    # point should display -- cheaper than carrying two tooltips for every protein
    hv = list(dict.fromkeys(norm + lowr))
    pos = {lc: i for i, lc in enumerate(hv)}
    show_norm = [pos[lc] for lc in norm]
    show_lowr = [pos[lc] for lc in lowr]
    classes = (list(order) if order else
               (list(pd.unique(f[cls].dropna())) if cls else ["all"]))
    cmap = dict(colors or {})
    for i, c in enumerate(classes):
        cmap.setdefault(c, MUTED if c in muted else CLASS_COLORS[i % len(CLASS_COLORS)])

    rng = np.random.default_rng(3)
    jit = rng.normal(0, .07, len(f))
    xs, ys = f[x].to_numpy(float).copy(), f[y].to_numpy(float).copy()
    capped = np.zeros(len(f), bool)
    if cap is not None:
        cx, cy = xs >= cap - .01, ys >= cap - .01
        capped = cx | cy
        xs[cx] = cap + jit[cx]
        ys[cy] = cap + jit[cy]

    def cell(v):
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            return ""
        if isinstance(v, (float, np.floating)):
            return f"{v:.3g}"
        return str(v)

    pts = []
    for i, (_, r) in enumerate(f.iterrows()):
        pts.append(dict(
            g=str(r[gene]), x=round(float(xs[i]), 4), y=round(float(ys[i]), 4),
            c=classes.index(r[cls]) if cls and r[cls] in classes else 0,
            k=1 if capped[i] else 0,
            r=1 if str(r[gene]) in set(rings) else 0,
            lb=1 if (lower_bound and bool(r.get(lower_bound))) else 0,
            v=[cell(r[col]) for _, col in hv]))

    data = dict(
        pts=pts, classes=[str(c) for c in classes],
        colors=[cmap[c] for c in classes],
        muted=[classes.index(m) for m in muted if m in classes],
        fields=[lab for lab, col in hv],
        showNorm=show_norm, showLower=show_lowr,
        labels=sorted(set(default_labels) & set(f[gene].astype(str))),
        xlabel=xlabel or x, ylabel=ylabel or y,
        cap=None if cap is None else float(cap), capNote=cap_note,
        diag=bool(diagonal),
        xlim=list(map(float, xlim)) if xlim else None,
        ylim=list(map(float, ylim)) if ylim else None)

    page = (_TPL.replace("__DATA__", _json.dumps(data, separators=(",", ":")))
                .replace("__TITLE__", _html.escape(title))
                .replace("__SUBTITLE__", subtitle))
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
 :root{--ink:#1a1d21;--ink2:#4a5056;--mut:#878d94;--line:#e3e6e9;--accent:#0072B2}
 *{box-sizing:border-box}
 body{margin:0;padding:16px 20px;color:var(--ink);background:#fff;
      font:13px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif}
 h1{font-size:18px;margin:0 0 2px}
 .sub{color:var(--mut);font-size:12px;margin-bottom:10px}
 #wrap{display:flex;gap:12px;align-items:flex-start}
 #left{flex:1 1 auto;min-width:0;position:relative;border:1px solid var(--line);
       border-radius:8px;overflow:hidden}
 #side{flex:0 0 212px;display:flex;flex-direction:column;gap:9px}
 svg{display:block;width:100%;height:72vh;min-height:470px;cursor:grab}
 svg.drag{cursor:grabbing}
 .card{border:1px solid var(--line);border-radius:8px;padding:9px 10px;background:#f7f8f9}
 .card h2{font-size:11px;margin:0 0 6px;text-transform:uppercase;letter-spacing:.04em;
          color:var(--ink2)}
 .row{display:flex;gap:7px;align-items:center;padding:3px 4px;border-radius:5px;
      cursor:pointer;font-size:12px}
 .row:hover{background:#eceff2} .row.off{opacity:.36}
 .sw{flex:0 0 10px;height:10px;border-radius:3px}
 .n{margin-left:auto;color:var(--mut);font-size:11px;font-variant-numeric:tabular-nums}
 #bar{display:flex;justify-content:flex-end;margin:0 224px 7px 0}
 #tools{display:flex;gap:7px;align-items:center}
 #bar #note{margin-left:2px}
 .grp{display:flex;background:#fff;border:1px solid var(--line);border-radius:7px;
      overflow:hidden;box-shadow:0 1px 3px rgba(0,0,0,.07)}
 .grp button{border:0;background:#fff;padding:5px 9px;font:12px inherit;cursor:pointer;
             color:var(--ink2);border-right:1px solid var(--line)}
 .grp button:last-child{border-right:0}
 .grp button:hover{background:#eef2f6} .grp button.on{background:var(--accent);color:#fff}
 input[type=search]{width:100%;font:12px inherit;padding:5px 7px;border:1px solid var(--line);
                    border-radius:6px}
 #hits{max-height:168px;overflow:auto;margin-top:6px}
 #hits div{padding:3px 5px;border-radius:4px;cursor:pointer;font-size:12px}
 #hits div:hover{background:#eceff2}
 #tip{position:fixed;display:none;background:#1a1d21;color:#fff;padding:7px 9px;
      border-radius:7px;font-size:11.5px;line-height:1.55;z-index:60;pointer-events:none;
      box-shadow:0 4px 18px rgba(0,0,0,.25)}
 #tip b{font-size:12.5px} #tip table{border-collapse:collapse;margin-top:4px}
 #tip td{padding:0 0 0 8px;color:#dfe4e8} #tip td:first-child{padding:0;color:#9aa4ad}
 .hint{color:var(--mut);font-size:11.5px;line-height:1.5}
 .axl{font-size:11px;fill:#5c6166} .tik{font-size:10.5px;fill:#878d94}
 .lab{font-size:10.5px;paint-order:stroke;stroke:#fff;stroke-width:2.6}
</style>
<h1>__TITLE__</h1>
<div class="sub">__SUBTITLE__</div>
<div id="bar">
  <div id="tools">
    <div class="grp"><button id="zin">＋</button><button id="zout">−</button>
      <button id="zfit">fit</button></div>
    <div class="grp"><button id="labs">labels: default</button></div>
    <span id="note" class="hint"></span>
  </div>
</div>
<div id="wrap">
  <div id="left"><svg id="svg"></svg></div>
  <div id="side">
    <div class="card"><h2>Classes</h2><div id="legend"></div></div>
    <div class="card"><h2>Evidence</h2><div id="legendb"></div></div>
    <div class="card"><h2>Find a gene</h2>
      <input type="search" id="q" placeholder="type a gene symbol…" autocomplete="off">
      <div id="hits"></div></div>
    <div class="card"><div class="hint">
      scroll to zoom &middot; drag to pan<br>hover a dot for its values<br>
      click a dot to pin its label</div></div>
  </div>
</div>
<div id="tip"></div>
<script>
const D=__DATA__, NS="http://www.w3.org/2000/svg";
const el=(t,a={})=>{const e=document.createElementNS(NS,t);for(const k in a)e.setAttribute(k,a[k]);return e;};
const $=id=>document.getElementById(id);
const esc=s=>String(s).replace(/[&<>"]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));

const P=D.pts, N=P.length;
const off=new Set(D.muted);
const offB=new Set();              // muted evidence kinds: 0 measured, 1 lower bound
const vis=i=>!off.has(P[i].c) && !offB.has(P[i].lb);
const pinned=new Set(D.labels.map(g=>P.findIndex(p=>p.g===g)).filter(i=>i>=0));
let labMode="default", view={k:1,x:0,y:0}, W=0,H=0;
const M={l:58,r:16,t:14,b:46};

const xs=P.map(p=>p.x), ys=P.map(p=>p.y);
const pad=.5;
const DOM=D.xlim&&D.ylim
  ? {x0:D.xlim[0],x1:D.xlim[1],y0:D.ylim[0],y1:D.ylim[1]}
  : {x0:Math.min(...xs)-pad, x1:Math.max(...xs)+pad,
     y0:Math.min(...ys)-pad, y1:Math.max(...ys)+pad};
// one scale for both axes: with the same quantity on each, an unequal aspect would
// tilt the diagonal and make equal binding look like an effect
const span=Math.max(DOM.x1-DOM.x0, DOM.y1-DOM.y0);
DOM.x1=DOM.x0+span; DOM.y1=DOM.y0+span;

const SVG=$("svg");
const gGrid=el("g"), gPts=el("g"), gLead=el("g"), gLab=el("g"), gAx=el("g");
SVG.appendChild(gGrid); SVG.appendChild(gPts); SVG.appendChild(gLead);
SVG.appendChild(gLab); SVG.appendChild(gAx);

function fit(){ const r=SVG.getBoundingClientRect(); W=r.width; H=r.height; }
const px=v=>M.l+(v-DOM.x0)/(DOM.x1-DOM.x0)*(W-M.l-M.r);
const py=v=>H-M.b-(v-DOM.y0)/(DOM.y1-DOM.y0)*(H-M.t-M.b);
const SX=i=>px(P[i].x)*view.k+view.x, SY=i=>py(P[i].y)*view.k+view.y;
const inv=(sx)=>DOM.x0+((sx-view.x)/view.k-M.l)/(W-M.l-M.r)*(DOM.x1-DOM.x0);
const invY=(sy)=>DOM.y0+(H-M.b-(sy-view.y)/view.k)/(H-M.t-M.b)*(DOM.y1-DOM.y0);

const els=[];
function build(){
  gPts.textContent="";els.length=0;
  P.forEach((p,i)=>{
    // hollow = the position is a lower bound, not a measurement
    const c=el("circle",{r:p.r?4.2:3.1,
      fill:p.lb?"none":D.colors[p.c],
      stroke:p.lb?D.colors[p.c]:(p.r?"#111":"#fff"),
      "stroke-width":p.lb?1.3:(p.r?1.3:.6),cursor:"pointer"});
    c.dataset.i=i; gPts.appendChild(c); els.push(c);
  });
}

const MAX_LABELS=260;      // beyond this the declutter is O(n^2) on a useless pile
function inView(i){
  const x=SX(i), y=SY(i);
  return x>=M.l-4 && x<=W-M.r+4 && y>=M.t-4 && y<=H-M.b+4;
}
function labelled(i){
  if(!vis(i)) return false;
  // "all" means all currently in view: labelling 1591 genes at once is never
  // readable, and zooming into a region is how you actually want to use it
  if(labMode==="all") return inView(i);
  if(labMode==="none") return false;
  return pinned.has(i);
}

function draw(){
  // grid + axes
  gGrid.textContent=""; gAx.textContent="";
  const step=Math.max(1,Math.round((DOM.x1-DOM.x0)/view.k/8));
  for(let v=Math.ceil(DOM.x0/step)*step; v<=DOM.x1; v+=step){
    const X=px(v)*view.k+view.x;
    if(X<M.l-1||X>W-M.r+1) continue;
    gGrid.appendChild(el("line",{x1:X,x2:X,y1:M.t,y2:H-M.b,stroke:"#eef0f2"}));
    const t=el("text",{x:X,y:H-M.b+15,class:"tik","text-anchor":"middle"});
    t.textContent=v; gAx.appendChild(t);
  }
  for(let v=Math.ceil(DOM.y0/step)*step; v<=DOM.y1; v+=step){
    const Y=py(v)*view.k+view.y;
    if(Y<M.t-1||Y>H-M.b+1) continue;
    gGrid.appendChild(el("line",{x1:M.l,x2:W-M.r,y1:Y,y2:Y,stroke:"#eef0f2"}));
    const t=el("text",{x:M.l-7,y:Y+3.5,class:"tik","text-anchor":"end"});
    t.textContent=v; gAx.appendChild(t);
  }
  if(D.diag){
    const a=Math.max(DOM.x0,DOM.y0), b=Math.min(DOM.x1,DOM.y1);
    gGrid.appendChild(el("line",{x1:px(a)*view.k+view.x,y1:py(a)*view.k+view.y,
      x2:px(b)*view.k+view.x,y2:py(b)*view.k+view.y,
      stroke:"#9aa0a6","stroke-width":1,"stroke-dasharray":"5 4"}));
  }
  if(D.cap!==null){
    const CX=px(D.cap)*view.k+view.x, CY=py(D.cap)*view.k+view.y;
    gGrid.appendChild(el("line",{x1:CX,x2:CX,y1:M.t,y2:H-M.b,stroke:"#cfd4d9","stroke-dasharray":"2 4"}));
    gGrid.appendChild(el("line",{x1:M.l,x2:W-M.r,y1:CY,y2:CY,stroke:"#cfd4d9","stroke-dasharray":"2 4"}));
    if(D.capNote){
      const t=el("text",{x:CX-6,y:H-M.b-7,class:"tik","text-anchor":"end"});
      t.textContent=D.capNote; gAx.appendChild(t);
    }
  }
  // axis titles
  const xt=el("text",{x:(M.l+W-M.r)/2,y:H-6,class:"axl","text-anchor":"middle"});
  xt.textContent=D.xlabel; gAx.appendChild(xt);
  const yt=el("text",{x:14,y:(M.t+H-M.b)/2,class:"axl","text-anchor":"middle",
                      transform:`rotate(-90 14 ${(M.t+H-M.b)/2})`});
  yt.textContent=D.ylabel; gAx.appendChild(yt);

  // points
  P.forEach((p,i)=>{
    const c=els[i];
    if(!vis(i)){ c.style.display="none"; return; }
    c.style.display="";
    c.setAttribute("cx",SX(i)); c.setAttribute("cy",SY(i));
  });
  placeLabels();
}

function placeLabels(){
  gLab.textContent=""; gLead.textContent="";
  const ids=[]; for(let i=0;i<N;i++) if(labelled(i)) ids.push(i);
  if(!ids.length) return;
  const fs=10.5, pos=ids.map(i=>[SX(i),SY(i)-9]);
  const w=ids.map(i=>P[i].g.length*fs*0.62+6), h=fs*1.5+4;
  const solve=ids.length<=MAX_LABELS;     // too many to separate; draw them as they lie
  // the zoom/label buttons float over the top-right of the plot; a label placed there
  // is simply invisible, so push any that land under them clear first
  const tb=$("tools").getBoundingClientRect(), sv=SVG.getBoundingClientRect();
  const KO={x0:tb.left-sv.left-8, y0:0, x1:W, y1:tb.bottom-sv.top+8};
  for(let k=0;k<ids.length;k++){
    if(pos[k][0]+w[k]/2>KO.x0 && pos[k][1]-h/2<KO.y1)
      pos[k][1]=KO.y1+h/2+2;
  }
  for(let it=0; solve && it<400; it++){
    let moved=false;
    for(let a=0;a<ids.length;a++) for(let b=a+1;b<ids.length;b++){
      const dx=pos[b][0]-pos[a][0], dy=pos[b][1]-pos[a][1];
      const ox=(w[a]+w[b])/2-Math.abs(dx), oy=h-Math.abs(dy);
      if(ox>0&&oy>0){ moved=true;
        if(oy*2.2<ox){const s=(oy/2+.5)*(dy>=0?1:-1);pos[a][1]-=s;pos[b][1]+=s;}
        else{const s=(ox/2+.5)*(dx>=0?1:-1);pos[a][0]-=s;pos[b][0]+=s;}
      }
    }
    if(!moved) break;
  }
  ids.forEach((i,k)=>{
    const [lx,ly]=pos[k];
    if(Math.hypot(lx-SX(i),ly-SY(i))>13)
      gLead.appendChild(el("line",{x1:SX(i),y1:SY(i),x2:lx,y2:ly+3,
                                   stroke:"#b9bec4","stroke-width":.7}));
    const t=el("text",{x:lx,y:ly,class:"lab","text-anchor":"middle",
                       fill:P[i].r?"#111":"#333","font-weight":P[i].r?700:400});
    t.textContent=P[i].g; gLab.appendChild(t);
  });
}

/* hover */
SVG.addEventListener("mousemove",ev=>{
  const i=ev.target.dataset&&ev.target.dataset.i!==undefined?+ev.target.dataset.i:-1;
  const tip=$("tip");
  if(i<0||drag){tip.style.display="none";return;}
  const p=P[i];
  let rows="";
  (p.lb?D.showLower:D.showNorm).forEach(k=>{
    if(p.v[k]!=="") rows+=`<tr><td>${esc(D.fields[k])}</td><td>${esc(p.v[k])}</td></tr>`; });
  tip.innerHTML=`<b>${esc(p.g)}</b><div style="color:#9aa4ad">${esc(D.classes[p.c])}`
                +(p.lb?" · position is a lower bound":"")+`</div>`
                +(rows?`<table>${rows}</table>`:"");
  tip.style.display="block";
  const w=tip.offsetWidth,h=tip.offsetHeight;
  tip.style.left=Math.min(ev.clientX+14,innerWidth-w-10)+"px";
  tip.style.top =Math.min(ev.clientY+14,innerHeight-h-10)+"px";
});
SVG.addEventListener("mouseleave",()=>$("tip").style.display="none");

/* click to pin, drag to pan */
let drag=null;
SVG.addEventListener("mousedown",ev=>{
  const i=ev.target.dataset&&ev.target.dataset.i!==undefined?+ev.target.dataset.i:-1;
  drag={x:ev.clientX-view.x,y:ev.clientY-view.y,i,moved:false};
  SVG.classList.add("drag");
});
window.addEventListener("mousemove",ev=>{
  if(!drag) return;
  if(Math.abs(ev.clientX-view.x-drag.x)>2||Math.abs(ev.clientY-view.y-drag.y)>2) drag.moved=true;
  if(drag.moved){ view.x=ev.clientX-drag.x; view.y=ev.clientY-drag.y; draw(); }
});
window.addEventListener("mouseup",()=>{
  SVG.classList.remove("drag");
  if(drag&&!drag.moved&&drag.i>=0){
    if(pinned.has(drag.i)) pinned.delete(drag.i); else pinned.add(drag.i);
    if(labMode==="none"){labMode="default";$("labs").textContent="labels: default";}
    placeLabels();
  }
  drag=null;
});
SVG.addEventListener("wheel",ev=>{
  ev.preventDefault();
  const r=SVG.getBoundingClientRect(), cx=ev.clientX-r.left, cy=ev.clientY-r.top;
  const f=ev.deltaY<0?1.14:1/1.14, k2=Math.min(40,Math.max(.6,view.k*f)), g=k2/view.k;
  view.x=cx-(cx-view.x)*g; view.y=cy-(cy-view.y)*g; view.k=k2; draw();
},{passive:false});
$("zin").onclick=()=>{const f=1.3,k2=Math.min(40,view.k*f),g=k2/view.k;
  view.x=W/2-(W/2-view.x)*g;view.y=H/2-(H/2-view.y)*g;view.k=k2;draw();};
$("zout").onclick=()=>{const f=1/1.3,k2=Math.max(.6,view.k*f),g=k2/view.k;
  view.x=W/2-(W/2-view.x)*g;view.y=H/2-(H/2-view.y)*g;view.k=k2;draw();};
$("zfit").onclick=()=>{view={k:1,x:0,y:0};draw();};
$("labs").onclick=()=>{
  labMode=labMode==="default"?"all":labMode==="all"?"none":"default";
  $("labs").textContent="labels: "+labMode; placeLabels();
  const n=document.querySelectorAll("text.lab").length;
  $("note").textContent = (labMode==="all"&&n>MAX_LABELS)
    ? `${n} labels in view — zoom in to separate them` : "";
};

/* legend */
function legend(){
  const b=$("legend"); b.textContent="";
  D.classes.forEach((c,i)=>{
    const n=P.filter(p=>p.c===i).length;
    const row=document.createElement("div");
    row.className="row"+(off.has(i)?" off":"");
    row.innerHTML=`<span class="sw" style="background:${D.colors[i]}"></span>
                   <span>${esc(c)}</span><span class="n">${n}</span>`;
    row.onclick=()=>{ if(off.has(i)) off.delete(i); else off.add(i);
                      row.classList.toggle("off"); draw(); };
    b.appendChild(row);
  });
}

function legendB(){
  const b=$("legendb"); b.textContent="";
  [[0,"measured","filled"],[1,"lower bound","hollow"]].forEach(([k,lab,shape])=>{
    const n=P.filter(p=>p.lb===k).length;
    if(!n) return;
    const row=document.createElement("div");
    row.className="row"+(offB.has(k)?" off":"");
    const sw = k ? `<span class="sw" style="background:#fff;border:1.6px solid #6b7177"></span>`
                 : `<span class="sw" style="background:#6b7177"></span>`;
    row.innerHTML=`${sw}<span>${lab} <span class="count">(${shape})</span></span>`
                 +`<span class="n">${n}</span>`;
    row.onclick=()=>{ if(offB.has(k)) offB.delete(k); else offB.add(k);
                      row.classList.toggle("off"); draw(); };
    b.appendChild(row);
  });
}

/* search */
$("q").addEventListener("input",()=>{
  const v=$("q").value.trim().toLowerCase(), box=$("hits"); box.textContent="";
  if(v.length<1) return;
  const hit=[]; for(let i=0;i<N;i++) if(P[i].g.toLowerCase().includes(v)) hit.push(i);
  hit.slice(0,40).forEach(i=>{
    const r=document.createElement("div");
    r.textContent=P[i].g;
    r.onclick=()=>{ pinned.add(i); off.delete(P[i].c); offB.delete(P[i].lb);
      // centre the view on the hit so a pinned gene is never off-screen
      view.x=W/2-px(P[i].x)*view.k; view.y=H/2-py(P[i].y)*view.k;
      legend(); legendB(); draw(); };
    box.appendChild(r);
  });
});

function boot(){ fit(); build(); legend(); legendB(); draw(); }
window.addEventListener("resize",()=>{fit();draw();});
boot();
</script>
"""
