"""py3Dmol view of a structure coloured by the consolidated core/surface/interface annotation.

    from viewer import view_annotated
    view_annotated(results, "NLGN1", pdb="5OJ6", chain="A").show()
    view_annotated(results, "NLGN1", pdb="3B3Q", show_dimer=True,
                   partners={"dimer": ["B"], "NRXN1B": ["E"]}).show()
"""
import py3Dmol
from structure import fetch_pdb, sifts_chains

IFACE_COLORS = {"iface_NRXN1B": "green", "iface_dimer": "blue", "iface_MDGA1": "purple"}


def _residue_sets(results, gene, pdb, chain, iface_colors):
    """Core / interface / other residue lists (PDB numbering) for one annotated chain."""
    r = results[gene]
    per = r["per_structure"]
    sub = per[per.structure == f"{pdb}:{chain}"]
    if sub.empty:
        raise ValueError(f"{pdb}:{chain} not in results[{gene!r}]; available: {sorted(per.structure.unique())}")
    m = sub[["pdb_resnum", "human_resnum"]].merge(r["table"], left_on="human_resnum", right_index=True)
    core = [int(x) for x in m.pdb_resnum[m["class"] == "core"]]
    exposed = m[m["class"] != "core"]
    groups, used = {}, set()
    for c in [c for c in iface_colors if c in m.columns]:     # first listed colour wins on overlap
        resi = [int(x) for x in exposed.pdb_resnum[exposed[c]] if x not in used]
        groups[c] = resi; used.update(resi)
    other = [int(x) for x in exposed.pdb_resnum if x not in used]
    return core, groups, other


def _style_chain(v, ch, core, groups, other, iface_colors, core_color, backbone_color,
                 surface_color, opacity, outline, outline_width, outline_backbone):
    """Apply the full annotation styling (core cartoon + translucent surfaces) to one chain."""
    ol = {"color": "black", "width": outline_width}
    bb = {"color": backbone_color, "thickness": 0.2}
    if outline and outline_backbone: bb["outline"] = ol
    v.setStyle({"chain": ch}, {"cartoon": bb})
    cc = {"color": core_color}
    if outline: cc["outline"] = ol
    v.setStyle({"chain": ch, "resi": core}, {"cartoon": cc})
    v.addSurface(py3Dmol.VDW, {"opacity": opacity, "color": surface_color}, {"chain": ch, "resi": other})
    for c, resi in groups.items():
        if resi:
            v.addSurface(py3Dmol.VDW, {"opacity": opacity, "color": iface_colors[c]}, {"chain": ch, "resi": resi})


def _highlight_residues(v, ch, mapping, highlight, sidechain, label, stick_radius, seq=None):
    """mapping: {human_resnum: pdb_resnum}. Colours each highlighted residue, optionally drawing
    its side chain as sticks and labelling it."""
    drawn = {}
    for pos, color in highlight.items():
        pdb_resi = mapping.get(pos)
        if pdb_resi is None:
            continue
        sel = {"chain": ch, "resi": [int(pdb_resi)]}
        v.addStyle(sel, {"cartoon": {"color": color}})
        if sidechain:                                  # backbone N/C/O dropped so the R group reads
            v.addStyle({**sel, "not": {"atom": ["N", "C", "O"]}},
                       {"stick": {"color": color, "radius": stick_radius}})
        if label:
            aa = seq[pos - 1] if seq and 0 < pos <= len(seq) else ""
            v.addLabel(f"{aa}{pos}", {"position": {}, "backgroundColor": color, "fontColor": "white",
                                      "fontSize": 12, "backgroundOpacity": 0.85,
                                      "alignment": "center"}, sel)
        drawn[pos] = int(pdb_resi)
    return drawn


def view_annotated(results, gene, pdb, chain="A", iface_colors=IFACE_COLORS,
                   core_color="brown", backbone_color="#d9d9d9", surface_color="white",
                   opacity=0.6, outline=True, outline_width=0.06, outline_backbone=False,
                   show_dimer=False, dimer_chains=None, dimer_core_color=None,
                   show_partners=True, partners=None, partner_colors=None,
                   partner_color="#bbbbbb", partner_opacity=0.8,
                   highlight=None, highlight_sidechain=True, highlight_label=True,
                   highlight_stick_radius=0.22, width=800, height=600):
    """Core residues: cartoon (core_color). Surface + boundary residues: translucent surface
    mesh, coloured by interface membership (iface_colors), else surface_color.

    outline: black edge on the cartoon only (surfaces stay clean); outline_backbone extends it
      to the thin backbone cartoon as well.
    show_dimer: render the other NLGN copy/copies with the same full annotation styling (the
      dimer shares residue numbering with `chain`). Chains come from dimer_chains, else
      partners["dimer"], else every chain SIFTS maps to the same UniProt entry as `chain`.
      dimer_core_color: core colour for the second copy (default: same as core_color).
    show_partners: draw the remaining chains as thin cartoons.
    partners: {label: [chains]} as in structure_pipeline.STRUCTURES -> each partner chain is
      coloured to match its interface colour (iface_<label> in iface_colors).
    partner_colors: {chain: colour} explicit override; unlisted chains get partner_color.
    highlight: {residue number: colour} in HUMAN UniProt numbering (the same numbering as
      results[gene]["table"]); the residue's cartoon is recoloured, its side chain drawn as
      sticks (highlight_sidechain), and it is labelled (highlight_label).
    Residue numbers are taken from `results[gene]["per_structure"]` for this pdb:chain so the
    PDB's own numbering is used, whatever species it is."""
    core, groups, other = _residue_sets(results, gene, pdb, chain, iface_colors)
    style = dict(iface_colors=iface_colors, backbone_color=backbone_color, surface_color=surface_color,
                 opacity=opacity, outline=outline, outline_width=outline_width, outline_backbone=outline_backbone)

    v = py3Dmol.view(width=width, height=height, options={"shadow": True})
    v.addModel(fetch_pdb(pdb), "pdb")
    v.setStyle({}, {})                                                        # hide everything
    _style_chain(v, chain, core, groups, other, core_color=core_color, **style)

    # human resnum -> this structure's PDB numbering
    sub = results[gene]["per_structure"]
    sub = sub[sub.structure == f"{pdb}:{chain}"]
    h2p = dict(zip(sub.human_resnum, sub.pdb_resnum))
    hseq = results[gene].get("human_sequence")
    if highlight:
        _highlight_residues(v, chain, h2p, highlight, highlight_sidechain, highlight_label,
                            highlight_stick_radius, hseq)

    dimer = []
    if show_dimer:
        if dimer_chains: dimer = list(dimer_chains)
        elif partners and "dimer" in partners: dimer = list(partners["dimer"])
        else:
            info = sifts_chains(pdb)
            dimer = [c for c in info.index if c != chain and info.loc[c, "accession"] == info.loc[chain, "accession"]]
        for ch in dimer:
            _style_chain(v, ch, core, groups, other, core_color=dimer_core_color or core_color, **style)
            if highlight:
                _highlight_residues(v, ch, h2p, highlight, highlight_sidechain, False,
                                    highlight_stick_radius, hseq)

    if show_partners:                                                         # remaining chains, thin cartoons
        colors = {}
        for lab, chs in (partners or {}).items():
            for ch in chs: colors[ch] = iface_colors.get(f"iface_{lab}", partner_color)
        colors.update(partner_colors or {})
        shown = [chain, *dimer]
        rest = {"not": {"or": [{"chain": c} for c in shown]}}
        v.setStyle(rest, {"cartoon": {"color": partner_color, "thickness": 0.15, "opacity": partner_opacity}})
        for ch, col in colors.items():
            if ch not in shown:
                v.setStyle({"chain": ch}, {"cartoon": {"color": col, "thickness": 0.15, "opacity": partner_opacity}})

    v.zoomTo({"or": [{"chain": c} for c in [chain, *dimer]]})
    return v
