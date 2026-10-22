"""AlphaFold (or PDB) structure coloured by a per-residue score, e.g. a variant-density delta.

    from af_viewer import view_scores
    view_scores(delta, gene="NLGN3", start=1).show()

The score is mapped through a diverging colormap pinned at `center` (0 by default, so the
sign of a difference is what the colour reports) and handed to 3Dmol as a
{residue number: '#rrggbb'} map via `colorscheme={'prop':'resi','map':...}`.
"""
import json as _json
import matplotlib.colors as _mcolors
import matplotlib.pyplot as _plt
import py3Dmol

from netutil import urlopen_retry as _urlopen, pkg_path as _pkg_path


def uniprot_id(gene, organism=9606):
    """Gene symbol -> reviewed Swiss-Prot accession for that organism."""
    url = ("https://rest.uniprot.org/uniprotkb/search?query="
           f"gene_exact:{gene}+AND+organism_id:{organism}+AND+reviewed:true"
           "&fields=accession&format=json&size=1")
    results = _json.loads(_urlopen(url))["results"]
    if not results:
        raise ValueError(f"no reviewed UniProt entry for {gene!r} (organism {organism})")
    return results[0]["primaryAccession"]


def alphafold_cif(accession, cache_dir="af_cache"):
    """AlphaFold DB model for a UniProt accession, as CIF text. Cached on disk."""
    import os
    d = _pkg_path(cache_dir)
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, f"{accession}.cif")
    if os.path.exists(path):
        return open(path).read()
    meta = _json.loads(_urlopen(f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"))
    cif = _urlopen(meta[0]["cifUrl"])
    cif = cif.decode() if isinstance(cif, bytes) else cif
    open(path, "w").write(cif)
    return cif


def score_colors(values, start=1, cmap="seismic_r", center=0.0, vmin=None, vmax=None,
                 n_steps=256, nan_color="#d9d9d9"):
    """Per-residue scores -> {residue number: '#rrggbb'} for 3Dmol's resi colorscheme.

    start:  residue number of values[0] (1 for a whole-protein convolve_v4 output)
    cmap:   'seismic_r' runs red -> white -> blue, so with center=0 a negative score
            (ClinVar-enriched, constrained) is red and a positive one (gnomAD-enriched,
            tolerant) is blue.  Use 'seismic' to swap, or any matplotlib colormap name.
    center: value pinned to the middle of the colormap.  None uses a plain min-max
            stretch instead, which is what you want for an unsigned score.
    """
    import numpy as np

    v = np.asarray(list(values), dtype=float)
    lo = float(np.nanmin(v)) if vmin is None else float(vmin)
    hi = float(np.nanmax(v)) if vmax is None else float(vmax)
    if center is None:
        norm = _mcolors.Normalize(vmin=lo, vmax=hi)
    elif lo < center < hi:
        norm = _mcolors.TwoSlopeNorm(vcenter=float(center), vmin=lo, vmax=hi)
    else:                       # all one side of centre -- keep centre at the midpoint
        half = max(abs(lo - center), abs(hi - center)) or 1.0
        norm = _mcolors.Normalize(vmin=center - half, vmax=center + half)
    cm = _plt.get_cmap(cmap, n_steps)
    out = {}
    for i, val in enumerate(v):
        resi = start + i
        out[resi] = nan_color if np.isnan(val) else _mcolors.to_hex(cm(norm(val)))
    return out


def view_scores(values, gene=None, accession=None, cif=None, pdb=None, start=1,
                cmap="seismic_r", center=0.0, vmin=None, vmax=None, chain=None,
                style="cartoon", background="white", width=800, height=600,
                return_colors=False):
    """Load a structure and colour every residue by `values`.

    Give one of: gene (looked up on UniProt, then AlphaFold), accession (AlphaFold),
    cif (CIF text you already have), or pdb (a 4-character PDB id).
    """
    if cif is None and pdb is None:
        accession = accession or uniprot_id(gene)
        cif = alphafold_cif(accession)

    colors = score_colors(values, start=start, cmap=cmap, center=center,
                          vmin=vmin, vmax=vmax)

    v = py3Dmol.view(width=width, height=height)
    if pdb is not None:
        from structure import fetch_pdb
        v.addModel(fetch_pdb(pdb), "pdb")
    else:
        v.addModel(cif, "cif")

    sel = {} if chain is None else {"chain": chain}
    v.setStyle(sel, {style: {"colorscheme": {"prop": "resi", "map": colors}}})
    v.zoomTo(sel)
    return (v, colors) if return_colors else v


def score_colorbar(values, cmap="seismic_r", center=0.0, vmin=None, vmax=None,
                   label="score", figsize=(4, 0.45)):
    """Standalone colourbar matching `score_colors`, for a figure legend."""
    import numpy as np
    v = np.asarray(list(values), dtype=float)
    lo = float(np.nanmin(v)) if vmin is None else float(vmin)
    hi = float(np.nanmax(v)) if vmax is None else float(vmax)
    if center is None or not (lo < center < hi):
        norm = _mcolors.Normalize(vmin=lo, vmax=hi)
    else:
        norm = _mcolors.TwoSlopeNorm(vcenter=float(center), vmin=lo, vmax=hi)
    fig, ax = _plt.subplots(figsize=figsize)
    fig.colorbar(_plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax,
                 orientation="horizontal", label=label)
    return ax
