"""Export an HTML figure (alignment, table) to PNG or PDF for slides and documents.

    from export_figure import export_html
    export_html(h, "alignment.png", scale=3)       # h = the IPython HTML object, or a .html path
    export_html("alignment.html", "alignment.pdf") # PDF is vector -> stays crisp at any zoom
"""
import os, shutil, subprocess, tempfile

CHROME_CANDIDATES = [
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "/Applications/Chromium.app/Contents/MacOS/Chromium",
    "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
    shutil.which("google-chrome"), shutil.which("chromium"), shutil.which("chromium-browser"),
]

def find_chrome():
    for c in CHROME_CANDIDATES:
        if c and os.path.exists(c):
            return c
    raise RuntimeError("No Chrome/Chromium found; install one or pass chrome=/path/to/binary")

def export_html(source, out, fmt=None, width=1400, height=1000, scale=3, padding=28,
                background="#ffffff", chrome=None, timeout=90):
    """source: an IPython HTML object, a raw HTML string, or a path to an .html file.
    out:     destination .png or .pdf (fmt inferred from the extension).
    scale:   device pixel ratio for PNG -- 3 gives a crisp image at slide size.
    PDF output is vector: text stays selectable and sharp at any zoom."""
    fmt = fmt or os.path.splitext(out)[1].lstrip(".").lower()
    body = getattr(source, "data", source)
    if isinstance(body, str) and os.path.exists(body):
        html_path, tmp = os.path.abspath(body), None
    else:
        doc = (f"<!doctype html><meta charset='utf-8'>"
               f"<body style='margin:0;background:{background}'>"
               f"<div style='display:inline-block;padding:{padding}px'>{body}</div>")
        tmp = tempfile.NamedTemporaryFile("w", suffix=".html", delete=False)
        tmp.write(doc); tmp.close(); html_path = tmp.name

    exe = chrome or find_chrome()
    common = [exe, "--headless=new", "--disable-gpu", "--no-first-run",
              f"--default-background-color={background.lstrip('#')}"]
    if fmt == "pdf":
        cmd = common + [f"--print-to-pdf={os.path.abspath(out)}", "--no-pdf-header-footer",
                        f"--virtual-time-budget=10000", f"file://{html_path}"]
    else:
        cmd = common + [f"--screenshot={os.path.abspath(out)}",
                        f"--window-size={width},{height}",
                        f"--force-device-scale-factor={scale}",
                        f"--hide-scrollbars", f"--virtual-time-budget=10000",
                        f"file://{html_path}"]
    subprocess.run(cmd, capture_output=True, timeout=timeout)
    if tmp: os.unlink(tmp.name)
    if not os.path.exists(out):
        raise RuntimeError(f"export failed; try a larger width/height or check {exe}")
    return out


# ---------------------------------------------------------------- RTF (editable text) export
def _rtf_escape(s):
    return s.replace("\\", r"\\").replace("{", r"\{").replace("}", r"\}")

def alignment_rtf(aln, colored, out=None, width=60, offsets=None, conservation=True,
                  font="Courier New", fontsize=10, key=None, title=None):
    """Render an alignment as RTF -- editable text with background highlighting, which pastes
    into Word and PowerPoint with formatting intact.

    aln:      Bio.Align.MultipleSeqAlignment
    colored:  {(seq_id, alignment column): "#rrggbb"} (e.g. from glyco_align.feature_columns)
    key:      [(label, "#rrggbb"), ...] rendered as a legend line above the alignment
    out:      path to write; also returns the RTF string
    """
    from msa_print import conservation_char
    palette, ctab = {}, []
    def cidx(hexcol):
        h = hexcol.lstrip("#")
        if h not in palette:
            ctab.append((int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)))
            palette[h] = len(ctab) + 1          # RTF colour indices are 1-based; 1 = auto/black
        return palette[h]
    for c in colored.values(): cidx(c)
    for _, c in (key or []): cidx(c)

    offsets = offsets or {}
    names = [r.id for r in aln]
    pad = max(len(n) for n in names) + 2
    L = aln.get_alignment_length()
    cons = "".join(conservation_char(aln[:, i]) for i in range(L))
    count = {n: offsets.get(n, 0) for n in names}

    body = []
    if title:
        body.append(r"{\b " + _rtf_escape(title) + r"}\par\par")
    if key:
        parts = [r"{\b Key}  "]
        for label, col in key:
            i = cidx(col)
            parts.append(r"{\chcbpat%d\cb%d      }  %s   " % (i, i, _rtf_escape(label)))
        body.append("".join(parts) + r"\par\par")
    for start in range(0, L, width):
        stop = min(start + width, L)
        for rec in aln:
            chunk = str(rec.seq[start:stop])
            count[rec.id] += len(chunk.replace("-", ""))
            line = [_rtf_escape(f"{rec.id + ':':<{pad}}")]
            run_col, run = None, []
            def flush():
                if not run: return
                txt = _rtf_escape("".join(run))
                line.append(f"{{\\chcbpat{palette[run_col.lstrip('#')]}\\cb{palette[run_col.lstrip('#')]} {txt}}}"
                            if run_col else txt)
                run.clear()
            for i, ch in enumerate(chunk):
                col = colored.get((rec.id, start + i))
                if col != run_col:
                    flush(); run_col = col
                run.append(ch)
            flush()
            line.append(_rtf_escape(" " * (width - len(chunk)) + f" {count[rec.id]}"))
            body.append("".join(line) + r"\par")
        if conservation:
            body.append(_rtf_escape(f"{'':<{pad}}{cons[start:stop]}") + r"\line")
        body.append(r"\par")

    colortbl = "".join(f"\\red{r}\\green{g}\\blue{b};" for r, g, b in ctab)
    rtf = (r"{\rtf1\ansi\ansicpg1252\deff0"
           r"{\fonttbl{\f0\fmodern\fcharset0 " + font + r";}}"
           r"{\colortbl;" + colortbl + "}"
           r"\f0\fs" + str(int(fontsize * 2)) + " " + "".join(body) + "}")
    if out:
        open(out, "w", encoding="ascii", errors="xmlcharrefreplace").write(rtf)
    return rtf
