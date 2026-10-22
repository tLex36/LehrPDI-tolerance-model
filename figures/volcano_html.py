"""Self-contained interactive volcano page with tables underneath.

    from volcano_html import volcano_report
    volcano_report(d, out="report.html", title="...", tables=[...])

The volcano is SVG with hover tooltips and click-to-toggle labels; tables are plain HTML
with per-row highlight rules. No CDN -- the file opens offline.
"""
import html as _html
import json as _json

import numpy as np
import pandas as pd

_TPL = r"""<!doctype html>
<meta charset="utf-8">
<title>__TITLE__</title>
<style>
 :root{--ink:#1f2328;--muted:#6b7280;--line:#e5e7eb;--surface:#ffffff}
 body{font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
      margin:0;padding:26px 30px 60px;color:var(--ink);background:var(--surface);max-width:1180px}
 h1{font-size:19px;margin:0 0 4px} h2{font-size:15px;margin:30px 0 8px}
 .sub{color:var(--muted);font-size:13px;margin-bottom:14px}
 .ctl{font-size:12.5px;color:var(--muted);margin:8px 0 2px}
 .ctl input{width:66px;font-size:12.5px} .ctl b{color:var(--ink)}
 #grid{display:flex;flex-wrap:wrap;gap:18px 26px;margin-top:6px}
 .panel{flex:0 0 auto}
 .ptitle{font-size:13px;font-weight:600;margin:0 0 2px}
 .pmeta{font-size:11.5px;color:var(--muted);margin:0 0 2px}
 svg{display:block}
 .bg{fill:#c9ccd1;fill-opacity:.55} .sig{stroke:#1f2328;stroke-width:.4;stroke-opacity:.4}
 .thr{stroke:#9ca3af;stroke-width:1;stroke-dasharray:5 4}
 .ax{stroke:#374151;stroke-width:.9}
 text.lab{font-size:11px;fill:var(--ink);paint-order:stroke;stroke:#fff;stroke-width:3px;
          stroke-linejoin:round}
 text.axl{font-size:11.5px;fill:var(--muted)}
 #gtip{position:absolute;display:none;z-index:70;pointer-events:none;background:#1f2328;
       color:#fff;border-radius:4px;padding:3px 8px;font-size:11.5px;white-space:nowrap}
 table{border-collapse:collapse;font-size:12.5px;margin-top:4px}
 th{text-align:left;font-weight:600;border-bottom:1.5px solid var(--line);padding:5px 14px 5px 0;
    color:var(--muted);white-space:nowrap}
 td{padding:3px 14px 3px 0;border-bottom:1px solid #f3f4f6;font-variant-numeric:tabular-nums}
 td.n{text-align:right} .cap{color:var(--muted);font-size:12px;margin:2px 0 8px}
 .key{display:inline-block;width:11px;height:11px;border-radius:2px;vertical-align:-1px;
      margin-right:5px;border:1px solid rgba(0,0,0,.12)}
 .legend{font-size:12px;color:var(--muted);margin:6px 0 10px}
 .legend span{margin-right:16px}
</style>
<h1>__TITLE__</h1>
<div class="sub">__SUBTITLE__</div>
<div class="ctl">
  <b>p &lt;</b> <input id="bp" type="number" step="0.01" min="0" max="1" value="__PTHR__">
  &nbsp; <b>__FCLABEL__</b>
  <input id="bfc" type="number" step="0.25" min="0" value="__FCTHR__">
  <span __CENHIDE__>&nbsp; <label title="centre the null on the median instead of zero">
    <input type="checkbox" id="cen" __CENCHECKED__> <b>centre on median</b></label></span>
  &nbsp; <i>click a point to toggle its label &mdash; labels apply to every panel</i>
</div>
<div id="grid"></div>
<div id="gtip"></div>
__TABLES__
<script>
const PANELS = __DATA__;
const $ = id => document.getElementById(id);
const NS = "http://www.w3.org/2000/svg";
const el = (t, a) => { const e = document.createElementNS(NS, t);
                       for (const k in a) e.setAttribute(k, a[k]); return e; };
let LAB = new Set(__DEFAULT_LABELS__);
const SINGLE = PANELS.length === 1;
const W = SINGLE ? 1120 : 545, H = SINGLE ? 560 : 360;
const M = SINGLE ? {l:66, r:30, t:22, b:52} : {l:52, r:22, t:18, b:46};

// one <div class=panel> per contrast, each with its own svg
const grid = $("grid");
PANELS.forEach((P, i) => {
  const box = document.createElement("div"); box.className = "panel";
  box.innerHTML = `<div class="ptitle">${P.title}</div>`
                + `<div class="pmeta" id="meta${i}"></div>`;
  const sv = document.createElementNS(NS, "svg"); sv.id = "vol" + i;
  box.appendChild(sv); grid.appendChild(box);
  sv.addEventListener("mousemove", ev => hover(ev, i));
  sv.addEventListener("mouseleave", () => { $("gtip").style.display = "none"; });
  sv.addEventListener("click", ev => {
    const p = near(ev, i); if (!p) return;
    LAB.has(p.d.g) ? LAB.delete(p.d.g) : LAB.add(p.d.g); draw();
  });
});

function draw(){ PANELS.forEach((P, i) => drawPanel(P, i)); }

function drawPanel(P, idx){
  const D = P.pts, MED = P.med;
  const pc = +$("bp").value, fc = +$("bfc").value;
  const NULLX = $("cen").checked ? MED : 0;
  const svg = $("vol" + idx); svg.setAttribute("width", W); svg.setAttribute("height", H);
  svg.innerHTML = "";
  const xs = D.map(d => d.x), ys = D.map(d => d.y);
  const xlo = Math.min(...xs), xhi = Math.max(...xs), yhi = Math.max(...ys, 1);
  const pad = (xhi - xlo) * .04 || 1;
  const X = v => M.l + (v - xlo + pad) / (xhi - xlo + 2*pad) * (W - M.l - M.r);
  const Y = v => H - M.b - v / (yhi * 1.06) * (H - M.t - M.b);
  if (pc > 0 && pc < 1){
    const yl = Y(-Math.log10(pc));
    svg.appendChild(el("line", {class:"thr", x1:M.l, x2:W-M.r, y1:yl, y2:yl}));
    const t = el("text", {x:W-M.r, y:yl-5, class:"axl", "text-anchor":"end"});
    t.textContent = "p = " + pc; svg.appendChild(t);
  }
  [NULLX-fc, NULLX+fc].forEach(v => svg.appendChild(
    el("line", {class:"thr", x1:X(v), x2:X(v), y1:M.t, y2:H-M.b})));
  // the median is always drawn for reference, highlighted when it is the active null
  const onMed = $("cen").checked;
  // when the median IS zero the two reference lines coincide -- draw only the null,
  // otherwise the "median" and "null = 0" labels overprint each other.
  const sameLine = Math.abs(MED) < 1e-9;
  if (!sameLine){
    svg.appendChild(el("line", {class:"thr", x1:X(MED), x2:X(MED), y1:M.t, y2:H-M.b,
                                stroke:"#16a34a", "stroke-width": onMed ? 1.6 : 1,
                                "stroke-dasharray": onMed ? "4 3" : "2 4",
                                "stroke-opacity": onMed ? 1 : .55}));
    const tm = el("text", {x:X(MED)+4, y:M.t+10, class:"axl", fill:"#16a34a"});
    tm.textContent = (onMed ? "null = median " : "median ") + MED.toFixed(2);
    svg.appendChild(tm);
  }
  if (!onMed){
    svg.appendChild(el("line", {class:"thr", x1:X(0), x2:X(0), y1:M.t, y2:H-M.b,
                                stroke:"#6b7280", "stroke-dasharray":"4 3"}));
    const t0 = el("text", {x:X(0)+4, y:M.t+10, class:"axl"});
    t0.textContent = "null = 0"; svg.appendChild(t0);
  }
  svg.appendChild(el("line", {class:"ax", x1:M.l, x2:W-M.r, y1:H-M.b, y2:H-M.b}));
  svg.appendChild(el("line", {class:"ax", x1:M.l, x2:M.l, y1:M.t, y2:H-M.b}));

  P.PTS = []; let nu = 0, nd = 0;
  const ringed = [];
  D.forEach(d => {
    const x = X(d.x), y = Y(d.y); P.PTS.push({x:x, y:y, d:d});
    const sig = d.p !== null && d.p < pc && Math.abs(d.x - NULLX) > fc;
    if (sig) (d.x > NULLX ? nu++ : nd++);
    const fill = d.mark ? d.mark
               : (sig ? (d.x > NULLX ? "#b2182b" : "#2166ac") : "#c9ccd1");
    if (d.ring){ ringed.push({x:x, y:y, fill:fill, ring:d.ring, sig:sig}); return; }
    svg.appendChild(el("circle", {cx:x, cy:y, r: sig ? 3.4 : 2, fill: fill,
                                  class: sig ? "sig" : "bg"}));
  });
  // outlined categories last, so their rings are never hidden by neighbouring points
  ringed.forEach(r => svg.appendChild(el("circle", {
    cx:r.x, cy:r.y, r: r.sig ? 4.2 : 3.4, fill: r.fill,
    "fill-opacity": r.sig ? 1 : .75,
    stroke: r.ring, "stroke-width": 2.1})));
  $("meta" + idx).textContent = nu + " up / " + nd + " down of " + D.length
    + "  ·  median " + MED.toFixed(2);
  // place labels: only nudge against labels that actually overlap in BOTH axes,
  // and draw a leader line whenever a label ends up away from its point
  const placed = [];
  D.filter(d => LAB.has(d.g)).sort((a,b) => b.y - a.y).forEach(d => {
    const px = X(d.x), py = Y(d.y);
    const w = d.g.length * 6.2 + 8;                 // approximate rendered text width
    let lx = px + 7, ly = py - 6, dir = -1, step = 0;
    const hits = (x, y) => placed.some(q =>
      Math.abs(y - q.y) < 11 && x < q.x + q.w && q.x < x + w);
    while (hits(lx, ly) && step < 14){               // alternate up/down, growing
      step++; dir = -dir;
      ly = py - 6 + dir * Math.ceil(step / 2) * 11;
    }
    if (hits(lx, ly)){ lx = px - w - 7; ly = py - 6; } // last resort: flip to the left
    placed.push({x: lx, y: ly, w: w});
    if (Math.abs(ly - (py - 6)) > 6 || lx < px){     // displaced -> connect it
      svg.appendChild(el("line", {x1: px, y1: py,
        x2: lx < px ? lx + w + 2 : lx - 2, y2: ly - 3,
        stroke: "#9ca3af", "stroke-width": .7}));
    }
    const t = el("text", {x: lx, y: ly, class: "lab"});
    t.textContent = d.g; svg.appendChild(t);
  });
  let t = el("text", {x:(M.l+W-M.r)/2, y:H-14, class:"axl", "text-anchor":"middle"});
  t.textContent = P.xlabel; svg.appendChild(t);
  t = el("text", {x:16, y:(M.t+H-M.b)/2, class:"axl", "text-anchor":"middle",
                  transform:"rotate(-90 16 " + ((M.t+H-M.b)/2) + ")"});
  t.textContent = "-log10 p"; svg.appendChild(t);
  [xlo, NULLX, xhi].forEach(v => { const q = el("text",
    {x:X(v), y:H-M.b+18, class:"axl", "text-anchor":"middle"});
    q.textContent = v.toFixed(1); svg.appendChild(q); });
  for (let k = 0; k <= yhi; k += Math.max(1, Math.ceil(yhi/(SINGLE?7:5)))){
    const q = el("text", {x:M.l-10, y:Y(k)+4, class:"axl", "text-anchor":"end"});
    q.textContent = k; svg.appendChild(q); }
}
function near(ev, idx){
  const r = $("vol" + idx).getBoundingClientRect();
  const mx = ev.clientX - r.left, my = ev.clientY - r.top;
  let best = null, bd = SINGLE ? 100 : 64;
  for (const p of (PANELS[idx].PTS || [])){
    const dd = (p.x-mx)**2 + (p.y-my)**2; if (dd < bd){ bd = dd; best = p; } }
  return best;
}
function hover(ev, idx){
  const p = near(ev, idx), lab = $("gtip");
  if (!p){ lab.style.display = "none"; return; }
  const a = p.d.x - ($("cen").checked ? PANELS[idx].med : 0);
  lab.textContent = p.d.g + "  " + (a>=0?"+":"") + a.toFixed(2)
                  + ($("cen").checked ? " (adj)" : "")
                  + (p.d.p === null ? "" : "  p=" + p.d.p.toExponential(1));
  lab.style.display = "block";
  lab.style.left = (ev.pageX + 12) + "px"; lab.style.top = (ev.pageY - 28) + "px";
}
$("bp").oninput = $("bfc").oninput = draw;
$("cen").onchange = draw;
draw();
</script>
"""


def _table_html(title, df, caption=None, highlight=None, legend=None):
    """highlight: list of (predicate(row) -> bool, css colour) applied in order."""
    head = "".join(f"<th>{_html.escape(str(c))}</th>" for c in df.columns)
    rows = []
    for _, r in df.iterrows():
        bg = ""
        for pred, colour in (highlight or []):
            if pred(r):
                bg = f' style="background:{colour}"'
                break
        cells = "".join(
            f'<td class="{"n" if isinstance(v,(int,float,np.floating)) and not isinstance(v,bool) else ""}">'
            f'{_html.escape("" if v is None or (isinstance(v,float) and np.isnan(v)) else str(v))}</td>'
            for v in r)
        rows.append(f"<tr{bg}>{cells}</tr>")
    leg = f'<div class="legend">{legend}</div>' if legend else ""
    cap = f'<div class="cap">{caption}</div>' if caption else ""
    return (f"<h2>{_html.escape(title)}</h2>{cap}{leg}"
            f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table>")


def volcano_report(d=None, panels=None, out=None, title="Volcano", subtitle="",
                   xlabel="log2 fold change", center=None, centre_on_median=True,
                   p_threshold=0.05, fc_threshold=1.0, default_labels=(), mark=None,
                   rings=None, tables=(), width="100%", height=1100,
                   show_centre=True, fc_label=None):
    """Build the page.

    d:      a single frame with gene / log2fc / pval, or
    panels: [{"title":..., "d":frame, "xlabel":..., "center":...}, ...] for a grid of
            contrasts sharing one set of controls. Clicking a point labels that gene in
            EVERY panel, which is how you compare a protein across contrasts.
    centre_on_median: initial state of the page's "centre on median" toggle; each panel
            uses its own median.
    rings:  {gene: colour} -- coloured outline, drawn above everything else.
    mark:   {gene: colour} -- replaces the point's fill.
    show_centre: False hides the "centre on median" control (and forces the null to 0).
            Centring assumes most proteins are unchanged between the two arms -- true for
            two lysates, FALSE for a bait-vs-isotype IP, where the design intends most
            proteins to differ. Pass center=0.0 per panel alongside it.
    fc_label: overrides the threshold label, e.g. "log2FC >" when not centring.
    """
    if panels is None:
        if d is None:
            raise ValueError("pass either d= or panels=")
        panels = [{"title": title, "d": d, "xlabel": xlabel, "center": center}]

    def pack(p):
        f = p["d"]
        med = float(f.log2fc.median()) if p.get("center") is None else float(p["center"])
        pts = [{"g": r.gene, "x": round(float(r.log2fc), 3),
                "y": 0.0 if pd.isna(r.pval) else round(float(-np.log10(max(r.pval, 1e-300))), 3),
                "p": None if pd.isna(r.pval) else float(r.pval),
                "mark": (mark or {}).get(r.gene), "ring": (rings or {}).get(r.gene)}
               for r in f.itertuples() if np.isfinite(r.log2fc)]
        return {"title": _html.escape(p["title"]), "xlabel": p.get("xlabel", xlabel),
                "med": round(med, 4), "pts": pts}

    payload = [pack(p) for p in panels]
    html = (_TPL.replace("__DATA__", _json.dumps(payload, separators=(",", ":")))
                .replace("__TITLE__", _html.escape(title))
                .replace("__SUBTITLE__", subtitle)
                .replace("__PTHR__", str(p_threshold))
                .replace("__FCTHR__", str(fc_threshold))
                .replace("__DEFAULT_LABELS__", _json.dumps(list(default_labels)))
                .replace("__CENCHECKED__", "checked" if centre_on_median and show_centre else "")
                .replace("__CENHIDE__", "" if show_centre else 'style="display:none"')
                .replace("__FCLABEL__", _html.escape(
                    fc_label or ("|log2FC − median| >" if show_centre else "log2FC >"))
                    .replace("log2FC", "log<sub>2</sub>FC").replace("&gt;", "&nbsp;&gt;"))
                .replace("__TABLES__", "".join(_table_html(**t) for t in tables)))
    if out:
        open(out, "w").write(html)
        print(f"wrote {out}  ({len(html)//1024} KB)")
    return Report(html, width, height, out)


class Report:
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
        return f"<Report{' -> ' + self.path if self.path else ''}>"
