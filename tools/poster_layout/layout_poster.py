"""Clean up the layout of the JHU Building Bridges poster (48 x 36 in, one slide).

    python layout_poster.py            # reads original.pptx, writes poster_v2.pptx (+ copy next to the original)

Grid: 0.75 in margins, 0.5 in gutters. Background 0.75-18.75 | Methods 19.25-47.25 (tolerance 19.25-28.75,
structural features 29.25-47.25) | Results 19.25-37.75 | Citations 38.25-47.25.
Converted figures (groups) are moved as rigid units; where one is resized its text is resized with it.
"""
import copy, os, shutil, subprocess, sys
from lxml import etree
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.oxml.ns import qn
from pptx.oxml import parse_xml
from pptx.opc.package import Part
from pptx.opc.constants import RELATIONSHIP_TYPE as RT

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from references import REFS

SRC = os.path.join(HERE, "original.pptx")
OUT = os.path.join(HERE, "poster_v2.pptx")
FINAL = os.path.expanduser("~/Documents/Work/Presentations/JHU_building_bridges_10.23.26_poster_layout.pptx")
POSTER = "/Users/awlehr/LehrPDI/figures/results/tolerance_model/poster"
HEAD = "2C3947"                        # header colour used throughout the poster
MUTED = "6B6A66"
A_NS = "http://schemas.openxmlformats.org/drawingml/2006/main"

prs = Presentation(SRC)
slide = prs.slides[0]
S = list(slide.shapes)                 # original order: indices below refer to it

# ------------------------------------------------------------------ helpers
I = lambda v: int(round(v * 914400))
inch = lambda v: v / 914400


def place(sh, x=None, y=None, w=None, h=None):
    if x is not None: sh.left = I(x)
    if y is not None: sh.top = I(y)
    if w is not None: sh.width = I(w)
    if h is not None: sh.height = I(h)


def box(sh):
    return inch(sh.left), inch(sh.top), inch(sh.width), inch(sh.height)


def move_by(shapes, dx, dy):
    for sh in shapes:
        if hasattr(sh, "begin_x") and sh.shape_type == 9:   # connector: move both ends
            bx, by, ex, ey = sh.begin_x, sh.begin_y, sh.end_x, sh.end_y
            sh.begin_x, sh.begin_y, sh.end_x, sh.end_y = bx + I(dx), by + I(dy), ex + I(dx), ey + I(dy)
        else:
            sh.left, sh.top = sh.left + I(dx), sh.top + I(dy)


def scale_fonts(el, f):
    for e in el.iter(qn("a:rPr"), qn("a:defRPr"), qn("a:endParaRPr")):
        if e.get("sz"):
            e.set("sz", str(max(100, int(round(int(e.get("sz")) * f)))))


def scale_shape(sh, f, anchor="tl"):
    """Uniform scale (pictures; groups also get their text scaled so labels keep fitting)."""
    x, y, w, h = box(sh)
    sh.width, sh.height = I(w * f), I(h * f)
    if sh.shape_type == 6:
        scale_fonts(sh._element, f)
    if anchor == "c":
        place(sh, x + (w - w * f) / 2, y + (h - h * f) / 2)


def fit_height(sh, h):
    x, y, w, h0 = box(sh)
    sh.width, sh.height = I(w * h / h0), I(h)


def set_text_style(sh, size=None, bold=None, font="Arial", color=None, italic=None):
    for p in sh.text_frame.paragraphs:
        for r in p.runs:
            if font: r.font.name = font
            if size: r.font.size = Pt(size)
            if bold is not None: r.font.bold = bold
            if italic is not None: r.font.italic = italic
            if color: r.font.color.rgb = __import__("pptx").dml.color.RGBColor.from_string(color)
        end = p._p.find(qn("a:endParaRPr"))
        if end is not None and size:
            end.set("sz", str(int(size * 100)))


def no_wrap(sh):
    sh.text_frame.word_wrap = False


def wrap(sh):
    sh.text_frame.word_wrap = True


def rewrite(paragraph, runs, size=None, font="Arial", bold=None, color=None):
    """Replace a paragraph's runs, keeping its paragraph properties. runs: [(text, dict(sup=, i=, b=))]."""
    p = paragraph._p
    old = p.findall(qn("a:r"))
    base = copy.deepcopy(old[0].find(qn("a:rPr"))) if old and old[0].find(qn("a:rPr")) is not None else etree.SubElement(etree.Element("x"), qn("a:rPr"))
    for r in old:
        p.remove(r)
    for e in p.findall(qn("a:br")):
        p.remove(e)
    end = p.find(qn("a:endParaRPr"))
    for text, o in runs:
        r = etree.SubElement(p, qn("a:r"))
        rp = copy.deepcopy(base); r.append(rp)
        if size: rp.set("sz", str(int(size * 100)))
        if bold is not None: rp.set("b", "1" if bold else "0")
        if o.get("b") is not None: rp.set("b", "1" if o["b"] else "0")
        if o.get("i"): rp.set("i", "1")
        if o.get("sup"): rp.set("baseline", "30000")
        if font:
            for tag in ("a:latin", "a:ea", "a:cs"):
                e = rp.find(qn(tag))
                if e is not None: rp.remove(e)
            lat = etree.SubElement(rp, qn("a:latin")); lat.set("typeface", font)
        if color:
            for e in rp.findall(qn("a:solidFill")): rp.remove(e)
            sf = etree.Element(qn("a:solidFill")); c = etree.SubElement(sf, qn("a:srgbClr")); c.set("val", color)
            rp.insert(0, sf)
        t = etree.SubElement(r, qn("a:t")); t.text = text
        if text != text.strip():
            t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
        if end is not None:
            p.remove(end); p.append(end)


def sup(n):
    return (str(n), {"sup": True})


def delete(sh):
    el = sh._element; el.getparent().remove(el)


def add_svg_picture(svg_path, png_path, x, y, w, h, name):
    pic = slide.shapes.add_picture(png_path, I(x), I(y), I(w), I(h))
    pic.name = name
    attach_svg(pic, open(svg_path, "rb").read())
    return pic


def attach_svg(pic, svg_bytes, rid=None):
    if rid is None:
        partname = prs.part.package.next_partname("/ppt/media/image%d.svg")
        svg_part = Part(partname, "image/svg+xml", prs.part.package, svg_bytes)
        rid = slide.part.relate_to(svg_part, RT.IMAGE)
    blip = pic._element.xpath(".//a:blip")[0]
    for e in blip.findall(qn("a:extLst")):
        blip.remove(e)
    blip.append(parse_xml(
        f'<a:extLst xmlns:a="{A_NS}"><a:ext uri="{{96DAC541-7B7A-43D3-8B79-37D633B846F1}}">'
        f'<asvg:svgBlip xmlns:asvg="http://schemas.microsoft.com/office/drawing/2016/SVG/main" '
        f'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" r:embed="{rid}"/></a:ext></a:extLst>'))


def add_png_fallback(pic, png_path):
    """An SVG-only picture gets a PNG main image (kept crop), so every viewer shows it."""
    blip = pic._element.xpath(".//a:blip")[0]
    svg = next(e for e in blip.iter() if e.tag.endswith("}svgBlip"))
    svg_rid = svg.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}embed")
    image_part, rid = slide.part.get_or_add_image_part(png_path)
    blip.set(qn("r:embed"), rid)
    return svg_rid


CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"


def svg_to_png(svg_path, png_path, scale=2.5):
    """Render an SVG at its own aspect ratio (headless Chrome screenshot, cropped to the SVG's exact size)."""
    import re as _re, math, tempfile
    from PIL import Image
    head = open(svg_path).read(600)
    W = float(_re.search(r'width="([\d.]+)', head).group(1)); H = float(_re.search(r'height="([\d.]+)', head).group(1))
    html = os.path.join(tempfile.mkdtemp(), "svg.html")
    open(html, "w").write(f'<html><body style="margin:0;background:white"><img src="file://{svg_path}" '
                          f'style="display:block;width:{W}px;height:{H}px"></body></html>')
    raw = png_path + ".raw.png"
    subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--hide-scrollbars", f"--force-device-scale-factor={scale}",
                    f"--window-size={math.ceil(W)},{math.ceil(H) + 100}", f"--screenshot={raw}", f"file://{html}"],
                   capture_output=True, timeout=60)
    im = Image.open(raw).convert("RGB")
    im.crop((0, 0, int(round(W * scale)), int(round(H * scale)))).save(png_path)
    os.remove(raw)


def textbox(x, y, w, h, runs, size, align="l", color="000000", bold=False, italic=False, name=None):
    tb = slide.shapes.add_textbox(I(x), I(y), I(w), I(h))
    tf = tb.text_frame; tf.word_wrap = True
    for m in ("margin_left", "margin_right", "margin_top", "margin_bottom"):
        setattr(tf, m, 0)
    p = tf.paragraphs[0]
    p.alignment = {"l": 1, "c": 2, "r": 3}[align]
    for text, o in runs:
        r = p.add_run(); r.text = text
        r.font.name = "Arial"; r.font.size = Pt(size); r.font.bold = o.get("b", bold); r.font.italic = o.get("i", italic)
        r.font.color.rgb = __import__("pptx").dml.color.RGBColor.from_string(color)
        if o.get("sup"):
            r._r.get_or_add_rPr().set("baseline", "30000")
    if name: tb.name = name
    return tb


# ------------------------------------------------------------------ named shapes (original order)
(phylo_old, bar_methods, qr, flow, sub_genome, bar_background, logo_brown, logo_nih, gradient, title, questions,
 bar_methods_ext, sub_tolerance, sub_structure, h_surface, h_structured, h_domains, h_interface, hist_esm, hist_tol,
 m0, sub_esm, fig_r2b, fig_r2, m1, m3_tracks, residue_table, r1_table, message, esm_schematic, credit_brandes,
 esm_hist, arrow1, arrow2, text_scores, divider1, divider2, bar_results, bar_citations, citations, hist_title,
 auroc_note, m3b_slice, m3c_png, lab_a, lab_alpha, lab_beta, t_glycine, legend, af_model, aiupred_track,
 annot_track, white_mask, sources_table, domains_text, hist_caption, synapse, net_img, manhattan, pedigree,
 carrow1, carrow2, carrow3, credit_pintacuda, credit_grove, icons) = S

L, R = 0.75, 47.25
BG = (0.75, 18.75)                     # Background column
TOL = (19.25, 28.75)                   # Methods: tolerance
STR = (29.25, 47.25)                   # Methods: structural features
RES = (19.25, 37.75)                   # Results
CIT = (38.25, 47.25)                   # Citations
cx = lambda a: (a[0] + a[1]) / 2

# ------------------------------------------------------------------ header
place(title, x=24.0 - inch(title.width) / 2)
place(gradient, x=L)
place(logo_nih, x=L, y=6.0 - inch(logo_nih.height) / 2)
place(logo_brown, x=L + inch(logo_nih.width) + 0.35, y=6.0 - inch(logo_brown.height) / 2)
place(qr, x=R - inch(qr.width), y=6.0 - inch(qr.height) / 2)
paras = title.text_frame.paragraphs
for p in paras:                        # authors with affiliations; affiliation numbers as superscripts
    for r in list(p.runs):
        if r.text.startswith("Alexander W. Lehr"):
            anchor = r._r
            r.text = "Alexander W. Lehr"
            for text, is_sup in [("1,2", True), (", Katherine W. Roche", False), ("1", True)]:
                new = copy.deepcopy(anchor); new.find(qn("a:t")).text = text
                if is_sup: new.find(qn("a:rPr")).set("baseline", "30000")
                anchor.addnext(new); anchor = new
    for r in p.runs:
        if r.text.strip() in ("1", "2") and len(p.runs) > 1:
            r.text = r.text.strip()
            r._r.get_or_add_rPr().set("baseline", "30000")

# ------------------------------------------------------------------ section bars and column rules
for bar, (a, b) in [(bar_background, BG), (bar_methods, (TOL[0], STR[1])), (bar_results, RES), (bar_citations, CIT)]:
    place(bar, x=a, w=b - a)
    set_text_style(bar, size=48, bold=True)
delete(bar_methods_ext)
for d, xpos, y0, y1 in [(divider1, 19.0, 8.75, 35.25), (divider2, 29.0, 8.75, 28.3)]:
    d.begin_x, d.begin_y, d.end_x, d.end_y = I(xpos), I(y0), I(xpos), I(y1)

# ------------------------------------------------------------------ text styles: level-2 subheads, level-3 heads, body
for sh in (sub_genome, sub_tolerance, sub_structure, sub_esm):
    set_text_style(sh, size=40, bold=True, color=HEAD); no_wrap(sh)
rewrite(sub_esm.text_frame.paragraphs[0], [("Evolutionary Scale Modelling", {}), sup(1)], size=40, bold=True, color=HEAD)
for sh in (h_structured, h_surface):
    ps = sh.text_frame.paragraphs
    for k, p in enumerate(ps):
        for r in p.runs:
            r.font.name = "Arial"; r.font.size = Pt(26 if k == 0 else 20); r.font.bold = k == 0
            r.font.color.rgb = __import__("pptx").dml.color.RGBColor.from_string(HEAD)
    wrap(sh)
rewrite(h_structured.text_frame.paragraphs[1], [("AIUPred3 disorder score prediction", {}), sup(10)], size=20, bold=False, color=HEAD)
for sh in (h_domains, h_interface):
    set_text_style(sh, size=26, bold=True, color=HEAD); wrap(sh)
rewrite(t_glycine.text_frame.paragraphs[0], [("Glycine Binding Interface", {}), sup(17)], size=20, bold=False, color=HEAD)
wrap(t_glycine)

# ------------------------------------------------------------------ Background column
place(sub_genome, x=cx(BG) - inch(sub_genome.width) / 2, y=8.75)
place(flow, x=cx(BG) - inch(flow.width) / 2, y=9.75)
place(questions, x=cx(BG) - inch(questions.width) / 2, y=18.05)
place(sub_esm, x=cx(BG) - 4.75, y=19.4, w=9.5, h=0.77)
fx0, fy0, fw0, fh0 = box(esm_hist)
fit_height(esm_schematic, 2.53 * 10.2 / 12.01); place(esm_schematic, x=L, y=20.4)
place(esm_hist, x=1.3, y=22.75)
hx, hy, hw, hh = box(esm_hist)
for arrow, (u, v), y_start in [(arrow1, (0.756, 0.499), 23.72), (arrow2, (0.498, 0.44), 23.42)]:
    arrow.begin_x, arrow.begin_y = I(5.45), I(y_start)
    arrow.end_x, arrow.end_y = I(hx + u * hw), I(hy + v * hh)
rewrite(text_scores.text_frame.paragraphs[0], [("Scores from ClinVar", {}), sup(3), (" pathogenic mutations are lower than scores from mutations in gnomAD", {}), sup(4)], size=20, color="000000")
place(text_scores, x=5.55, y=23.1, w=5.35, h=1.35); wrap(text_scores)
rewrite(credit_brandes.text_frame.paragraphs[0], [("Adapted from Brandes et al., 2023", {}), sup(2)], size=14, color=MUTED)
for r in credit_brandes.text_frame.paragraphs[0].runs: r.font.italic = True
place(credit_brandes, x=5.55, y=25.62, w=5.35, h=0.32); wrap(credit_brandes)
phylo = add_svg_picture(os.path.join(POSTER, "svg", "I0_conservation_phylogeny_poster.svg"),
                        os.path.join(HERE, "assets", "I0_conservation_phylogeny_poster.png"),
                        BG[1] - 7.5, 20.35, 7.5, 5.85, "Conservation phylogeny (I0, poster size)")
delete(phylo_old)
rewrite(message.text_frame.paragraphs[0], [
    ("Synaptic proteins", {}), sup(5), (" are highly retained in evolution", {}), sup(6),
    (" but highly sensitive to mutations in disease. This indicates that they could be a good avenue for investigating "
     "human-specific mutations that lead to disease. There are many annotated pathogenic genes, and even more unannotated ones.", {})],
    size=26, bold=True, color="000000")
for p in message.text_frame.paragraphs[1:]:
    for r in p.runs: r.text = ""
place(message, x=BG[0] + 0.15, y=26.45, w=BG[1] - BG[0] - 0.3, h=1.45); wrap(message)
# bottom collage: three data sources -> synapse -> disorders
stack_cx = L + 2.75
place(net_img, x=stack_cx - inch(net_img.width) / 2, y=28.3)
place(manhattan, x=stack_cx - inch(manhattan.width) / 2, y=31.3)
place(pedigree, x=max(L, stack_cx - inch(pedigree.width) / 2), y=33.78)
for cr, img, n_, label in [(credit_pintacuda, net_img, 7, "Pintacuda et al., 2023"), (credit_grove, manhattan, 8, "Grove et al., 2019")]:
    rewrite(cr.text_frame.paragraphs[0], [(label, {}), sup(n_)], size=12, color=MUTED)
    for r in cr.text_frame.paragraphs[0].runs: r.font.italic = True
    x_, y_, w_, h_ = box(img)
    place(cr, x=x_, y=y_ + h_ + 0.02, w=w_, h=0.25); wrap(cr)
    cr.text_frame.paragraphs[0].alignment = 3
place(synapse, x=6.85, y=28.75)
sx, sy, sw, sh_ = box(synapse)
textbox(sx, sy + sh_ + 0.04, sw, 0.25, [("Bemben et al., 2015", {}), sup(9)], 12, align="c", color=MUTED, italic=True,
        name="Credit synapse (check)")
place(icons, x=BG[1] - inch(icons.width) - 0.1, y=31.85 - inch(icons.height) / 2)
nx, ny, nw, nh = box(net_img); mx, my, mw, mh = box(manhattan); px_, py_, pw_, ph_ = box(pedigree)
for c, (b0, b1, e0, e1) in [(carrow1, (nx + nw + 0.05, ny + nh * 0.55, sx - 0.02, sy + 2.75)),
                            (carrow2, (mx + mw + 0.05, my + mh * 0.5, sx - 0.02, sy + 3.15)),
                            (carrow3, (px_ + pw_ + 0.05, py_ + ph_ * 0.45, sx + 0.04, sy + 3.6))]:
    c.begin_x, c.begin_y, c.end_x, c.end_y = I(b0), I(b1), I(e0), I(e1)
# disease labels: Arial, nudged below each icon's own small printed caption
g_off = icons._element.find(qn("p:grpSpPr")).find(qn("a:xfrm"))
ch_ext = g_off.find(qn("a:chExt")); ext = g_off.find(qn("a:ext"))
k_y = int(ch_ext.get("cy")) / int(ext.get("cy"))            # child units per slide unit
k_x = int(ch_ext.get("cx")) / int(ext.get("cx"))
for lab in icons.shapes:                # each label centred under its icon, wrapping to two lines if long
    if lab.has_text_frame and lab.text_frame.text.strip():
        set_text_style(lab, size=14, color="1F1F1F")
        wrap(lab)
        c_ = lab.left + lab.width // 2
        lab.width = int(I(1.75) * k_x); lab.left = c_ - lab.width // 2
        lab.top = lab.top + int(I(0.13) * k_y)
        for para in lab.text_frame.paragraphs: para.alignment = 2

# ------------------------------------------------------------------ Methods: tolerance column
place(sub_tolerance, x=cx(TOL) - inch(sub_tolerance.width) / 2, y=8.75)
place(m0, x=cx(TOL) - inch(m0.width) / 2, y=9.7)
place(m1, x=cx(TOL) - inch(m1.width) / 2, y=13.62)
set_text_style(hist_title, size=18, bold=True, color="000000"); wrap(hist_title)
hist_title.text_frame.paragraphs[0].alignment = 2
place(hist_title, x=TOL[0], y=16.4, w=TOL[1] - TOL[0], h=0.7)
gap = 0.3
place(hist_tol, x=cx(TOL) - (inch(hist_tol.width) + inch(hist_esm.width) + gap) / 2, y=17.15)
place(hist_esm, x=inch(hist_tol.left) + inch(hist_tol.width) + gap, y=17.15)
set_text_style(hist_caption, size=16, color="000000"); wrap(hist_caption)
place(hist_caption, x=TOL[0], y=20.15, w=TOL[1] - TOL[0], h=0.92)
place(m3_tracks, x=cx(TOL) - inch(m3_tracks.width) / 2, y=21.17)

# ------------------------------------------------------------------ Methods: structural features
place(sub_structure, x=cx(STR) - inch(sub_structure.width) / 2, y=8.75)
COL_A = (STR[0], STR[0] + 7.3)
scale_shape(af_model, 0.93)
COL_B = (COL_A[1] + 0.15, COL_A[1] + 0.15 + inch(af_model.width))
COL_C = (COL_B[1] + 0.12, STR[1])
# column A: structured vs unstructured, then domains & motifs (track + mask + table move together)
place(h_structured, x=COL_A[0], y=9.7, w=COL_A[1] - COL_A[0], h=0.95)
place(aiupred_track, x=cx(COL_A) - inch(aiupred_track.width) / 2, y=10.75)
place(h_domains, x=COL_A[0], y=12.4, w=COL_A[1] - COL_A[0], h=0.5)
rewrite(domains_text.text_frame.paragraphs[0], [
    ("To identify critical residues in unstructured regions, we used domain annotation and motif detecting datasets", {}),
    sup("11–15"), (" to corroborate potentially critical unstructured regions, where structure doesn't readily reveal function", {})],
    size=16, bold=False, color=HEAD)
place(domains_text, x=COL_A[0], y=12.93, w=COL_A[1] - COL_A[0], h=1.2); wrap(domains_text)
ax0, ay0, _, _ = box(annot_track)
dx, dy = (cx(COL_A) - inch(annot_track.width) / 2) - ax0, 14.3 - ay0
move_by([annot_track, white_mask, sources_table], dx, dy)
# column B: AlphaFold model (+ citation on its title)
place(af_model, x=COL_B[0], y=9.7)
for sh in af_model.shapes:
    pass
def add_sup_after(group, needle, n_):
    for el in group._element.iter(qn("a:r")):
        t = el.find(qn("a:t"))
        if t is not None and t.text and needle in t.text:
            new = copy.deepcopy(el); new.find(qn("a:t")).text = str(n_)
            rp = new.find(qn("a:rPr"))
            if rp is None: rp = etree.SubElement(new, qn("a:rPr")); new.insert(0, rp)
            rp.set("baseline", "30000")
            el.addnext(new)
            return True
    return False
add_sup_after(af_model, "AlphaFold model", 16)
# column C: surface vs core, then interface (glycine) with its labels and legend moving together
place(h_surface, x=COL_C[0], y=9.7, w=COL_C[1] - COL_C[0], h=1.25)
place(m3b_slice, x=cx(COL_C) - inch(m3b_slice.width) / 2, y=11.05)
place(h_interface, x=COL_C[0], y=14.82, w=COL_C[1] - COL_C[0], h=0.5)
place(t_glycine, x=COL_C[0], y=15.3, w=COL_C[1] - COL_C[0], h=0.38)
t_glycine.text_frame.paragraphs[0].alignment = 2
px0, py0, _, _ = box(m3c_png)
dx, dy = (cx(COL_C) - inch(m3c_png.width) / 2) - px0, 15.85 - py0
move_by([m3c_png, lab_a, lab_alpha, lab_beta, legend], dx, dy)
for lab, w_ in ((lab_a, 1.4), (lab_alpha, 1.35), (lab_beta, 1.3)):
    set_text_style(lab, size=12); no_wrap(lab)
    lab.width = I(w_)
# encoded-features table: enlarged to fill the band under the three columns
scale_shape(residue_table, 1.22)
top_rt = max(inch(sources_table.top) + inch(sources_table.height), inch(legend.top) + inch(legend.height),
             inch(af_model.top) + inch(af_model.height)) + 0.3
place(residue_table, x=cx(STR) - inch(residue_table.width) / 2, y=(top_rt + 28.3) / 2 - inch(residue_table.height) / 2)

# ------------------------------------------------------------------ Results and Citations
scale_shape(r1_table, 0.94)
place(r1_table, x=RES[0], y=30.15)
set_text_style(auroc_note, size=16, color="000000"); wrap(auroc_note)
place(auroc_note, x=RES[0], y=30.15 + inch(r1_table.height) + 0.08, w=inch(r1_table.width), h=0.6)
work = os.path.join(HERE, "assets")
for pic, svg_name in [(fig_r2, "image12.svg"), (fig_r2b, "image11.svg")]:
    src_svg = os.path.join(HERE, "unpacked", "ppt", "media", svg_name)
    png = os.path.join(work, svg_name.replace(".svg", "_exact.png"))
    if not os.path.exists(png):
        svg_to_png(src_svg, png)
    add_png_fallback(pic, png)
    fit_height(pic, 3.75)
x_fig = inch(r1_table.left) + inch(r1_table.width) + 0.3
place(fig_r2, x=x_fig, y=30.2)
place(fig_r2b, x=RES[1] - inch(fig_r2b.width), y=30.2)
# citations: numbered, APA style
tf = citations.text_frame
place(citations, x=CIT[0] + 0.1, y=30.1, w=CIT[1] - CIT[0] - 0.15, h=5.1)
wrap(citations)
body = tf._txBody
for p in body.findall(qn("a:p")):
    body.remove(p)
for k, runs in enumerate(REFS):
    p = etree.SubElement(body, qn("a:p"))
    pPr = etree.SubElement(p, qn("a:pPr")); pPr.set("marL", str(I(0.32))); pPr.set("indent", str(-I(0.32)))
    spc = etree.SubElement(pPr, qn("a:spcAft")); etree.SubElement(spc, qn("a:spcPts")).set("val", "200")
    bf = etree.SubElement(pPr, qn("a:buFont")); bf.set("typeface", "Arial")
    an = etree.SubElement(pPr, qn("a:buAutoNum")); an.set("type", "arabicPeriod")
    for text, it in runs:
        r = etree.SubElement(p, qn("a:r")); rp = etree.SubElement(r, qn("a:rPr"))
        rp.set("lang", "en-US"); rp.set("sz", "1000"); rp.set("dirty", "0")
        if it: rp.set("i", "1")
        sf = etree.SubElement(rp, qn("a:solidFill")); etree.SubElement(sf, qn("a:srgbClr")).set("val", "1F1F1F")
        etree.SubElement(rp, qn("a:latin")).set("typeface", "Arial")
        t = etree.SubElement(r, qn("a:t")); t.text = text
        if text != text.strip():
            t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
bodyPr = body.find(qn("a:bodyPr"))
for a_ in ("lIns", "rIns", "tIns", "bIns"):
    bodyPr.set(a_, "0")

prs.save(OUT)
print("saved", OUT)
if "--final" in sys.argv:
    shutil.copy(OUT, FINAL); print("copied to", FINAL)
