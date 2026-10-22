"""Residue track: a reference sequence as a row of cells, with variants stacked above and below.

    from tracks import variant_track
    variant_track(seq[100:140], start=101,
                  above=["p.R112C", "p.R112H", 118],
                  below=["p.G130S"])
"""
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

REF_COLOR   = "#c9f0dd"      # light mint green -- the reference row
ABOVE_COLOR = "#6495ed"      # cornflower blue  -- the first variant list
BELOW_COLOR = "#e06666"      # red              -- the second variant list
EDGE, LW    = "#333333", 0.8

def _positions(items, start=1):
    """Variant positions -> (column, substituted residue, full label).

    Accepts ints, HGVS strings ('p.R112C', 'c.334C>T', 'R112C'), or (position, text) pairs.
    The substituted residue is taken from the HGVS alt field when present; for a
    (position, text) pair the text is used, so you can supply your own cell character."""
    from exons import parse_variant
    out = []
    for it in items or []:
        given = None
        if isinstance(it, (tuple, list)) and len(it) == 2:
            it, given = it
        if isinstance(it, (int, float)):
            pos, alt, label = int(it), (given or ""), (given or str(int(it)))
        else:
            q = parse_variant(str(it))
            pos = q["aa_pos"]
            if pos is None:
                continue
            alt = given if given is not None else (q["alt"] or "")
            label = str(it)
        out.append((pos - start, str(alt), label))
    return out

def variant_track(sequence, above=(), below=(), start=1, cell=0.34,
                  ref_color=REF_COLOR, above_color=ABOVE_COLOR, below_color=BELOW_COLOR,
                  show_letters=True, show_variant_letters=True, letter_size=11,
                  variant_letter_color="white", edgecolor=EDGE, linewidth=LW,
                  tick_every=0, tick_size=7, tick_marks=True,
                  label_above=False, label_below=False, label_size=7,
                  ax=None, savepath=None, title=None):
    """Reference row of `sequence`, with `above`/`below` variants stacked outward.

    sequence:   reference amino-acid string (slice it yourself for a window)
    above/below: variant positions -- ints, HGVS strings, or (position, label) pairs
    start:      residue number of sequence[0], so positions are real protein coordinates
    cell:       cell size in inches
    show_letters:         draw the reference residue in each reference cell
    show_variant_letters: draw the substituted residue inside each variant cell
    tick_every:           label every Nth residue position under the track (0 = off, default)
    tick_size/tick_marks: font size of those numbers, and whether to draw a short rule from
                          the track down to each one
    label_above/below:    also write the full variant name as a rotated label
    title:                optional; omitted by default
    """
    n = len(sequence)
    up, dn = _positions(above, start), _positions(below, start)
    up = [t for t in up if 0 <= t[0] < n]
    dn = [t for t in dn if 0 <= t[0] < n]

    def stack(items):
        levels, seen = [], {}
        for c, alt, lab in sorted(items):
            k = seen.get(c, 0); seen[c] = k + 1
            levels.append((c, k, alt, lab))
        return levels, (max(seen.values()) if seen else 0)

    up_cells, n_up = stack(up)
    dn_cells, n_dn = stack(dn)
    height = n_up + 1 + n_dn

    if ax is None:
        fig, ax = plt.subplots(figsize=(n * cell + 0.4, height * cell + 0.9))
        fig.patch.set_alpha(0)                       # clear background
    else:
        fig = ax.figure
    ax.patch.set_alpha(0)
    y_ref = n_dn                                     # reference row sits above the lower stack

    for c, ch in enumerate(sequence):                # reference row
        ax.add_patch(Rectangle((c, y_ref), 1, 1, facecolor=ref_color,
                               edgecolor=edgecolor, linewidth=linewidth, joinstyle="miter"))
        if show_letters:
            ax.text(c + 0.5, y_ref + 0.5, ch, ha="center", va="center",
                    fontsize=letter_size, family="DejaVu Sans Mono")

    fs = letter_size
    for c, k, alt, lab in up_cells:                  # variants above
        y = y_ref + 1 + k
        ax.add_patch(Rectangle((c, y), 1, 1, facecolor=above_color,
                               edgecolor=edgecolor, linewidth=linewidth, joinstyle="miter"))
        if show_variant_letters and alt:
            ax.text(c + 0.5, y + 0.5, alt, ha="center", va="center", fontsize=fs,
                    color=variant_letter_color, family="DejaVu Sans Mono", weight="bold")
        if label_above and lab:
            ax.text(c + 0.5, y_ref + 1 + n_up + 0.15, lab, rotation=45, rotation_mode="anchor",
                    ha="left", va="bottom", fontsize=label_size, color=above_color)
    for c, k, alt, lab in dn_cells:                  # variants below
        y = y_ref - 1 - k
        ax.add_patch(Rectangle((c, y), 1, 1, facecolor=below_color,
                               edgecolor=edgecolor, linewidth=linewidth, joinstyle="miter"))
        if show_variant_letters and alt:
            ax.text(c + 0.5, y + 0.5, alt, ha="center", va="center", fontsize=fs,
                    color=variant_letter_color, family="DejaVu Sans Mono", weight="bold")
        if label_below and lab:
            ax.text(c + 0.5, y_ref - n_dn - 0.15, lab, rotation=-45, rotation_mode="anchor",
                    ha="right", va="top", fontsize=label_size, color=below_color)

    if tick_every:
        y_base = y_ref - n_dn
        for c in range(n):
            pos = start + c
            if pos % tick_every == 0:
                if tick_marks:
                    ax.plot([c + 0.5, c + 0.5], [y_base, y_base - 0.35], color="#555",
                            lw=0.8, clip_on=False, zorder=1)
                ax.text(c + 0.5, y_base - 0.45, str(pos), ha="center", va="top",
                        fontsize=tick_size, color="#555")

    ax.set_xlim(-0.5, n + 0.5)
    ax.set_ylim(y_ref - n_dn - 1.2, y_ref + n_up + 1.2)
    ax.set_aspect("equal"); ax.axis("off")
    if title:
        ax.set_title(title, fontsize=10, loc="left", pad=6)
    if savepath:
        fig.savefig(savepath, bbox_inches="tight", pad_inches=0.05, transparent=True)
    return ax


# ---------------------------------------------------------------- lollipop variant
LOLLI_ABOVE_STEM, LOLLI_ABOVE_HEAD = "black", "#4a7fe0"      # black shaft, blue head
LOLLI_BELOW_STEM, LOLLI_BELOW_HEAD = "black", "#d62728"      # black shaft, red head

def lollipop_track(sequence, above=(), below=(), start=1, cell=0.34,
                   ref_color=REF_COLOR, show_letters=True, letter_size=11,
                   above_stem=LOLLI_ABOVE_STEM, above_head=LOLLI_ABOVE_HEAD,
                   below_stem=LOLLI_BELOW_STEM, below_head=LOLLI_BELOW_HEAD,
                   stem_width=1.4, head_size=90, scale_by_count=False, unit=0.8,
                   min_stem=0.9, edgecolor=EDGE, linewidth=LW, head_edge="black",
                   head_edge_width=0.6, tick_every=0, tick_size=7, tick_marks=True,
                   headroom=0.0, footroom=0.0,
                   ax=None, savepath=None, title=None):
    """Reference sequence as a row of cells, with a lollipop at every variant position.

    above/below: variant positions -- ints, HGVS strings, or (position, label) pairs
    scale_by_count: when True the stem lengthens with the number of variants at that residue;
                    default False -- every lollipop is the same height
    unit:        extra stem length per additional variant at a position
    head_size:   marker area for the ball
    headroom/footroom: extra vertical space (in cell units) reserved above/below the
                 lollipops, so a density curve drawn into the same axes afterwards is
                 not clipped by the y-limits this function sets
    ax:          draw into an existing axes instead of creating one; the axes is
                 returned either way (`fig = ax.figure` for the figure)
    """
    n = len(sequence)
    up, dn = _positions(above, start), _positions(below, start)
    up = [t for t in up if 0 <= t[0] < n]
    dn = [t for t in dn if 0 <= t[0] < n]

    def counts(items):
        d = {}
        for c, _, _ in items:
            d[c] = d.get(c, 0) + 1
        return d

    c_up, c_dn = counts(up), counts(dn)
    h_up = {c: min_stem + (k - 1) * unit if scale_by_count else min_stem
            for c, k in c_up.items()}
    h_dn = {c: min_stem + (k - 1) * unit if scale_by_count else min_stem
            for c, k in c_dn.items()}
    top = max(h_up.values(), default=0)
    bot = max(h_dn.values(), default=0)

    if ax is None:
        fig, ax = plt.subplots(figsize=(n * cell + 0.4,
                                        (top + bot + headroom + footroom + 1) * cell + 0.8))
        fig.patch.set_alpha(0)
    else:
        fig = ax.figure
    ax.patch.set_alpha(0)

    for c, ch in enumerate(sequence):                      # reference row, y = 0..1
        ax.add_patch(Rectangle((c, 0), 1, 1, facecolor=ref_color, edgecolor=edgecolor,
                               linewidth=linewidth, joinstyle="miter"))
        if show_letters:
            ax.text(c + 0.5, 0.5, ch, ha="center", va="center",
                    fontsize=letter_size, family="DejaVu Sans Mono")

    for c, h in sorted(h_up.items()):                      # above: stem then head
        ax.plot([c + 0.5, c + 0.5], [1, 1 + h], color=above_stem,
                lw=stem_width, solid_capstyle="round", zorder=2)
        ax.scatter([c + 0.5], [1 + h], s=head_size, color=above_head,
                   edgecolors=head_edge, linewidths=head_edge_width, zorder=3)
    for c, h in sorted(h_dn.items()):                      # below
        ax.plot([c + 0.5, c + 0.5], [0, -h], color=below_stem,
                lw=stem_width, solid_capstyle="round", zorder=2)
        ax.scatter([c + 0.5], [-h], s=head_size, color=below_head,
                   edgecolors=head_edge, linewidths=head_edge_width, zorder=3)

    if tick_every:
        for c in range(n):
            if (start + c) % tick_every == 0:
                if tick_marks:
                    ax.plot([c + 0.5, c + 0.5], [-bot - 0.55, -bot - 0.9], color="#555",
                            lw=0.8, clip_on=False, zorder=1)
                ax.text(c + 0.5, -bot - 1.0, str(start + c), ha="center", va="top",
                        fontsize=tick_size, color="#555")

    # leave room for the marker itself: heads sit at 1+top and -bot, and the circle has a
    # radius of roughly sqrt(head_size/pi) points, converted to data units via the cell size
    import math
    head_r = math.sqrt(head_size / math.pi) / (72 * cell) + 0.1
    ax.set_xlim(-0.5, n + 0.5)
    ax.set_ylim(-bot - head_r - footroom, 1 + top + head_r + headroom)
    ax.set_aspect("equal"); ax.axis("off")
    if title:
        ax.set_title(title, fontsize=10, loc="left", pad=6)
    if savepath:
        fig.savefig(savepath, bbox_inches="tight", pad_inches=0.1, transparent=True)
    return ax


def lollipops_on_axis(ax, variants, up=True, start=1, height=0.35,
                      scale_by_weight=False, matrix="grantham", max_height=0.8,
                      head_color="#4a7fe0", stem_color="black", head_size=28,
                      stem_lw=1.0, baseline=0.0, dedupe="max", zorder=4, label=None):
    """Draw missense variants as lollipops on an existing data-coordinate axes.

    Intended as an overlay on a convolve_v4 density plot: the stems rise from `baseline`
    (the zero line) to +height, or fall to -height when up=False.

    scale_by_weight: False -> every lollipop reaches `height`.
                     True  -> height is the substitution's normalised severity (0-1 from
                              `matrix`, the same weight convolve_v4 applies) times
                              `max_height`, so a conservative change sits low and a
                              radical one reaches the full extent.
    dedupe: how to collapse several variants at one residue -- 'max' (tallest wins),
            'sum' (heights add, clipped to max_height), or None (all drawn, overlapping).
    Non-missense entries (nonsense, synonymous, bare positions with no ref/alt) are
    skipped when scale_by_weight is on, since they have no matrix value.
    """
    from convolve import parse_substitution, substitution_weight

    sign = 1 if up else -1
    per_pos = {}
    for v in variants or ():
        ref, pos, alt = parse_substitution(v)
        if pos is None:
            continue
        if scale_by_weight:
            if ref is None or alt is None or alt in ("*", "="):
                continue
            h = substitution_weight(ref, alt, matrix=matrix, normalize=True) * max_height
        else:
            h = height
        per_pos.setdefault(pos, []).append(h)

    pts = []
    for pos, hs in sorted(per_pos.items()):
        if dedupe == "max":
            pts.append((pos, max(hs)))
        elif dedupe == "sum":
            pts.append((pos, min(sum(hs), max_height)))
        else:
            pts.extend((pos, h) for h in hs)
    if not pts:
        return []

    xs = [p for p, _ in pts]
    tops = [baseline + sign * h for _, h in pts]
    ax.vlines(xs, baseline, tops, color=stem_color, lw=stem_lw, zorder=zorder)
    ax.scatter(xs, tops, s=head_size, color=head_color, edgecolors="none",
               zorder=zorder + 1, label=label, clip_on=False)
    return pts


def density_on_track(ax, values, start=1, track_start=1, track_len=None, up=True,
                     baseline=None, amplitude=3.0, vmax=None, color=None, lw=1.6,
                     fill=True, fill_alpha=0.15, fill_color=None, zorder=1.5,
                     expand_ylim=True):
    """Draw a convolved density curve into a track axes, in the track's own coordinates.

    The track puts residue `track_start` in column 0 and the reference row between y=0 and
    y=1, so a curve is aligned simply by plotting at x = pos - track_start + 0.5.  Nothing
    about the axes is changed, which is what keeps the curve registered to the cells.

    values:      density, one entry per residue, values[0] being residue `start`
                 (pass the whole protein; it is sliced to the track window)
    track_start: residue number of the track's first column
    track_len:   number of columns; defaults to the axes' current x-limit
    up:          True  -> rises from the top of the reference row (y=1)
                 False -> inverted, falling from the bottom of the reference row (y=0)
    amplitude:   full height of the curve in cell units when a value equals `vmax`
    vmax:        value mapped to `amplitude`; default is this curve's own max.  Pass the
                 same vmax to both curves to put them on one comparable scale -- otherwise
                 each is independently normalised and their heights are not comparable.
    expand_ylim: grow the axes' y-limits if the curve would fall outside them.  The track
                 fixes its limits around the lollipops, so without this a curve drawn
                 afterwards is silently clipped unless headroom/footroom was reserved.
    """
    import warnings
    import numpy as _np

    v = _np.asarray(values, dtype=float)
    if track_len is None:
        track_len = int(round(ax.get_xlim()[1] - 0.5))
    lo = track_start - start                       # index into `values` of the first column
    idx = _np.arange(lo, lo + track_len)
    keep = (idx >= 0) & (idx < len(v))
    if not keep.any():
        warnings.warn(f"density_on_track: no overlap -- `values` covers residues {start}"
                      f"-{start + len(v) - 1} but the track covers {track_start}"
                      f"-{track_start + track_len - 1}. Check `start` (the residue of "
                      f"values[0]) against `track_start`.", stacklevel=2)
        return None
    xs = _np.arange(track_len)[keep] + 0.5
    vv = v[idx[keep]]

    top = float(vmax) if vmax is not None else float(_np.nanmax(vv))
    if not top:
        warnings.warn("density_on_track: vmax is 0, nothing to scale -- the curve is "
                      "flat zero over this window.", stacklevel=2)
        return None
    if _np.nanmax(_np.abs(vv)) / top < 0.02:
        warnings.warn(f"density_on_track: this window peaks at {_np.nanmax(vv):.4g}, "
                      f"only {100 * _np.nanmax(vv) / top:.1f}% of vmax={top:.4g} -- the "
                      f"curve will sit on the baseline. Drop the shared vmax or raise "
                      f"amplitude.", stacklevel=2)
    sign = 1 if up else -1
    base = baseline if baseline is not None else (1.0 if up else 0.0)
    ys = base + sign * (vv / top) * amplitude

    color = color or (LOLLI_ABOVE_HEAD if up else LOLLI_BELOW_HEAD)
    line, = ax.plot(xs, ys, color=color, lw=lw, zorder=zorder, solid_capstyle="round")
    if fill:
        ax.fill_between(xs, ys, base, color=fill_color or color, alpha=fill_alpha,
                        lw=0, zorder=zorder - 0.1)
    if expand_ylim:                                   # the track pins its limits; grow them
        y0, y1 = ax.get_ylim()
        pad = 0.15
        ax.set_ylim(min(y0, float(ys.min()) - pad), max(y1, float(ys.max()) + pad))
    return line
