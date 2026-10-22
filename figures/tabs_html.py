"""Small self-contained tabbed HTML report: a figure plus tables behind named tabs.

    import tabs_html as th
    th.tabbed_report([{"name": "2D150", "html": th.table_html("KEGG", df)},
                      {"name": "DelESEV3", "html": ...}],
                     figure=fig, out="report.html")

Used where two conditions share a layout and you want to flip between them rather than
scroll past both. The figure is embedded as inline SVG, so the page stays one file and
the text in the plot stays selectable and sharp at any zoom.
"""
import html as _html
import io as _io

import numpy as np
import pandas as pd


def fig_svg(fig, strip_size=True) -> str:
    """Matplotlib figure -> inline SVG string."""
    buf = _io.StringIO()
    fig.savefig(buf, format="svg", bbox_inches="tight")
    s = buf.getvalue()
    s = s[s.index("<svg"):]                      # drop the XML/doctype preamble
    if strip_size:                               # let CSS size it instead
        import re
        s = re.sub(r'(<svg[^>]*?)\swidth="[^"]*"', r"\1", s, count=1)
        s = re.sub(r'(<svg[^>]*?)\sheight="[^"]*"', r"\1", s, count=1)
    return s


def table_html(title, df, caption=None, highlight=None, legend=None, index=False):
    """A titled HTML table. `highlight` is [(predicate(row) -> bool, css colour)]."""
    cols = list(df.columns)
    head = "".join(f"<th>{_html.escape(str(c))}</th>" for c in cols)
    rows = []
    for _, r in df.iterrows():
        bg = ""
        for pred, colour in (highlight or []):
            if pred(r):
                bg = f' style="background:{colour}"'
                break
        cells = ""
        for v in r[cols]:
            num = isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, bool)
            txt = "" if v is None or (isinstance(v, float) and np.isnan(v)) else str(v)
            cells += f'<td class="{"n" if num else ""}">{_html.escape(txt)}</td>'
        rows.append(f"<tr{bg}>{cells}</tr>")
    cap = f'<div class="cap">{caption}</div>' if caption else ""
    leg = f'<div class="legend">{legend}</div>' if legend else ""
    ttl = f"<h2>{_html.escape(title)}</h2>" if title else ""
    return (f"{ttl}{cap}{leg}<div class='tw'><table><thead><tr>{head}</tr></thead>"
            f"<tbody>{''.join(rows)}</tbody></table></div>")


def tabbed_report(tabs, out=None, title="Report", subtitle="", figure=None,
                  intro=None, common=None, width="100%", height=1200):
    """tabs: [{"name": str, "html": str}]. `figure` is a matplotlib figure or SVG string.

    `intro` renders above the tabs, `common` below them -- for content that does not
    change with the active tab.
    """
    if figure is not None and not isinstance(figure, str):
        figure = fig_svg(figure)

    btns = "".join(
        f'<button class="tab{" on" if i == 0 else ""}" data-i="{i}">'
        f'{_html.escape(t["name"])}</button>' for i, t in enumerate(tabs))
    panes = "".join(
        f'<section class="pane" data-i="{i}" style="display:{"block" if i == 0 else "none"}">'
        f'{t["html"]}</section>' for i, t in enumerate(tabs))

    page = (_TPL.replace("__TITLE__", _html.escape(title))
                .replace("__SUBTITLE__", subtitle)
                .replace("__FIG__", f'<div class="fig">{figure}</div>' if figure else "")
                .replace("__INTRO__", intro or "")
                .replace("__TABS__", btns)
                .replace("__PANES__", panes)
                .replace("__COMMON__", common or ""))
    if out:
        open(out, "w").write(page)
        print(f"wrote {out}  ({len(page)//1024} KB)")
    return Report(page, width, height, out)


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


_TPL = r"""<!doctype html>
<meta charset="utf-8">
<title>__TITLE__</title>
<style>
 :root{--ink:#1a1d21;--ink2:#4a5056;--mut:#878d94;--line:#e3e6e9;--bg:#fff;--accent:#0072B2}
 *{box-sizing:border-box}
 body{margin:0;padding:20px 24px 48px;background:var(--bg);color:var(--ink);
      font:13.5px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
      max-width:1280px}
 h1{font-size:20px;margin:0 0 3px}
 .sub{color:var(--mut);font-size:12.5px;margin-bottom:16px}
 h2{font-size:14px;margin:22px 0 7px}
 .cap{color:var(--ink2);font-size:12.5px;margin:0 0 9px;max-width:95ch}
 .legend{display:flex;gap:15px;flex-wrap:wrap;font-size:12px;color:var(--ink2);margin:0 0 9px}
 .legend i.key{display:inline-block;width:11px;height:11px;border-radius:3px;
               border:1px solid #cfd4d9;margin-right:5px;vertical-align:-1px}
 .fig{margin:4px 0 20px}
 .fig svg{width:100%;height:auto;display:block}
 .tabs{display:flex;gap:3px;border-bottom:2px solid var(--line);margin:22px 0 0}
 .tab{border:0;background:none;padding:9px 18px;font:600 13.5px inherit;color:var(--mut);
      cursor:pointer;border-bottom:2px solid transparent;margin-bottom:-2px}
 .tab:hover{color:var(--ink2)}
 .tab.on{color:var(--accent);border-bottom-color:var(--accent)}
 .pane{padding-top:4px}
 .tw{overflow:auto;max-height:620px;border:1px solid var(--line);border-radius:7px}
 table{border-collapse:collapse;width:100%;font-size:12px}
 th,td{text-align:left;padding:5px 8px;border-bottom:1px solid var(--line);vertical-align:top}
 th{background:#f7f8f9;color:var(--ink2);font-weight:600;position:sticky;top:0;z-index:1}
 tbody tr:hover{background:#f4f7fa}
 td.n{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
 code{background:#f2f4f6;padding:1px 4px;border-radius:3px;font-size:11.5px}
</style>
<h1>__TITLE__</h1>
<div class="sub">__SUBTITLE__</div>
__FIG__
__INTRO__
<div class="tabs">__TABS__</div>
__PANES__
__COMMON__
<script>
document.querySelectorAll(".tab").forEach(b=>{
  b.onclick=()=>{
    document.querySelectorAll(".tab").forEach(x=>x.classList.toggle("on",x===b));
    document.querySelectorAll(".pane").forEach(p=>{
      p.style.display = p.dataset.i===b.dataset.i ? "block" : "none";
    });
  };
});
</script>
"""


def page_tabs(pages, out=None, title="Figures", subtitle="", height=880,
              width="100%", note=None):
    """Combine several standalone HTML pages into one tabbed document.

    pages: {name: Page | html string | path to .html}

    Each page keeps its own <script>, globals and CSS, so they are embedded as sandboxed
    iframes rather than concatenated -- every one of these pages declares `const D`, and
    merging them into a single document would collide on the first tab.

    Inactive panes are stacked absolutely inside a sized container and faded with
    `opacity`, never `display:none` or `visibility:hidden`. A display:none iframe has
    zero width, so its page boots with a zero-size canvas and lays out nothing; and a
    visibility:hidden frame can still be blank when you switch to it before it has
    painted, with the switch alone not forcing a repaint. At opacity 0 the frame is laid
    out and painted like any other, so every tab is ready the moment it is shown.
    """
    items = []
    for name, p in pages.items():
        html = getattr(p, "html", None)
        if html is None:
            html = (open(p).read() if isinstance(p, str) and p.endswith(".html")
                    and _os.path.exists(p) else str(p))
        items.append((str(name), html))

    btns = "".join(
        f'<button class="ptab{" on" if i == 0 else ""}" data-i="{i}">{_html.escape(n)}</button>'
        for i, (n, _) in enumerate(items))
    panes = "".join(
        f'<div class="ppane{" on" if i == 0 else ""}" data-i="{i}">'
        f'<iframe srcdoc="{_html.escape(h, quote=True)}" sandbox="allow-scripts"></iframe>'
        f'</div>' for i, (_, h) in enumerate(items))

    page = (_PTPL.replace("__TITLE__", _html.escape(title))
                 .replace("__SUBTITLE__", subtitle)
                 .replace("__NOTE__", note or "")
                 .replace("__TABS__", btns)
                 .replace("__PANES__", panes)
                 .replace("__H__", str(int(height))))
    if out:
        open(out, "w").write(page)
        print(f"wrote {out}  ({len(page)//1024} KB)")
    return Report(page, width, height + 90, out)


import os as _os

_PTPL = r"""<!doctype html>
<meta charset="utf-8">
<title>__TITLE__</title>
<style>
 :root{--ink:#1a1d21;--ink2:#4a5056;--mut:#878d94;--line:#e3e6e9;--accent:#0072B2}
 *{box-sizing:border-box}
 body{margin:0;padding:16px 20px 24px;color:var(--ink);background:#fff;
      font:13px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif}
 h1{font-size:19px;margin:0 0 2px}
 .sub{color:var(--mut);font-size:12.5px;margin-bottom:12px}
 .note{color:var(--ink2);font-size:12.5px;margin:0 0 12px;max-width:100ch}
 .tabs{display:flex;gap:3px;border-bottom:2px solid var(--line);margin-bottom:10px}
 .ptab{border:0;background:none;padding:9px 18px;font:600 13.5px inherit;color:var(--mut);
       cursor:pointer;border-bottom:2px solid transparent;margin-bottom:-2px}
 .ptab:hover{color:var(--ink2)}
 .ptab.on{color:var(--accent);border-bottom-color:var(--accent)}
 #stack{position:relative;height:__H__px}
 .ppane{position:absolute;inset:0;opacity:0;pointer-events:none;z-index:0}
 .ppane.on{opacity:1;pointer-events:auto;z-index:1}
 .ppane iframe{width:100%;height:100%;border:1px solid var(--line);border-radius:8px;
                background:#fff}
</style>
<h1>__TITLE__</h1>
<div class="sub">__SUBTITLE__</div>
<div class="note">__NOTE__</div>
<div class="tabs">__TABS__</div>
<div id="stack">__PANES__</div>
<script>
const tabs=[...document.querySelectorAll(".ptab")], panes=[...document.querySelectorAll(".ppane")];
tabs.forEach(b=>b.onclick=()=>{
  tabs.forEach(x=>x.classList.toggle("on",x===b));
  panes.forEach(p=>p.classList.toggle("on", p.dataset.i===b.dataset.i));
});
</script>
"""
