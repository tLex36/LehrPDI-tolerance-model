"""Reference species tree (TimeTree median divergence ages, Mya) for the NLGN panel.

    tree = build_tree()                       # Bio.Phylo tree, branch lengths in Mya
    draw_species_tree(tree, identity=res["identity_matrix"]["Human"])
"""
from Bio.Phylo.BaseTree import Tree, Clade
from Bio import Phylo
import matplotlib.pyplot as plt

# (age_Mya, [children]) for internal nodes; bare strings for tips.
# Ages are TimeTree.org medians, rounded. Edit freely; branch lengths are derived.
TREE = (684, [
    (620, [
        (600, ["Drosophila", "C. elegans"], "Ecdysozoa"),
        "Octopus",
    ], "Protostomia"),
    (435, [
        "Zebrafish",
        (352, [
            "Xenopus tropicalis",
            (312, [
                (95, ["Chicken", "Crow"], "Aves"),
                (102, [
                    "Elephant",
                    (94, [
                        (78, [
                            (77, ["Dog", "Horse"]),
                            (62, ["Pig", (56, ["Cow", "Dolphin"])], "Artiodactyla"),
                        ], "Laurasiatheria"),
                        (88, [
                            (72, ["Guinea pig", (12, ["Mouse", "Rat"])], "Rodentia"),
                            (43, ["Marmoset", (29, ["Rhesus macaque", (6, ["Human", "Chimpanzee"])])], "Primates"),
                        ], "Euarchontoglires"),
                    ], "Boreoeutheria"),
                ], "Mammalia"),
            ], "Amniota"),
        ], "Tetrapoda"),
    ], "Vertebrata"),
], "Bilateria")

def _clade(node, parent_age):
    if isinstance(node, str):
        return Clade(branch_length=parent_age, name=node)
    age, children, *name = node                      # optional 3rd element = clade name
    c = Clade(branch_length=parent_age - age, clades=[_clade(ch, age) for ch in children])
    c.clade_name = name[0] if name else None         # kept separate from .name so Phylo.draw ignores it
    return c

def build_tree(spec=TREE, keep=None):
    """keep: optional set of tip names to retain (others pruned)."""
    root_age = spec[0]
    tree = Tree(root=_clade(spec, root_age), rooted=True)
    if keep:
        for tip in tree.get_terminals():
            if tip.name not in keep:
                tree.prune(tip)
    tree.root.branch_length = 0          # pruning can leave a stub on the root; drop it
    return tree

def draw_species_tree(tree, identity=None, title="Species tree", ax=None, labels=True):
    """identity: Series/dict {species: % identity to reference}, appended to tip labels."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 0.35 * len(tree.get_terminals()) + 1))
    def label(c):
        if not c.is_terminal() or not labels: return None
        if identity is not None and c.name in identity:
            return f"{c.name}  ({identity[c.name]:.0f}%)"
        return c.name
    Phylo.draw(tree, axes=ax, do_show=False, label_func=label)
    root_age = tree.distance(tree.root, tree.get_terminals()[0])
    ax.set_xlabel("Million years ago")
    step = 100 if root_age < 500 else 200
    ticks = [t for t in range(0, int(root_age) + 1, step)]
    ax.set_xticks([root_age - t for t in ticks]); ax.set_xticklabels(ticks)
    ax.set_ylabel(""); ax.set_yticks([]); ax.set_title(title)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    return ax


def node_positions(tree):
    """{clade: (x, y)} in Phylo.draw's coordinates: x = distance from root, tips at y = 1..n."""
    tips = tree.get_terminals()
    ypos = {t: i + 1 for i, t in enumerate(tips)}
    pos = {}
    for c in tree.find_clades(order="postorder"):
        y = ypos[c] if c.is_terminal() else sum(pos[ch][1] for ch in c.clades) / len(c.clades)
        pos[c] = (tree.distance(tree.root, c), y)
    return pos

def label_nodes(tree, ax, names=None, fontsize=8, dx=-4, dy=-0.15, **text_kw):
    """Write clade names at their nodes (right-aligned, just above-left of the node so they sit
    over the parent branch, not the descendants). names: subset to show (default: all named)."""
    pos = node_positions(tree)
    for c, (x, y) in pos.items():
        nm = getattr(c, "clade_name", None)
        if nm and (names is None or nm in names):
            ax.text(x + dx, y + dy, nm, fontsize=fontsize, va="bottom", ha="right",
                    style="italic", color="#444",
                    bbox=dict(facecolor="white", edgecolor="none", pad=1), **text_kw)

def clade_brackets(tree, ax, names=None, x=0.9, gap=0.15, fontsize=8, lw=1.2, pad=0.015):
    """Vertical bracket beside the tips for each named clade, label to its right.
    x and pad are in axes fraction (0 = left edge, 1 = right edge); widen xlim if it clips."""
    import matplotlib.transforms as mt
    tr = mt.blended_transform_factory(ax.transAxes, ax.transData)   # x: axes frac, y: data
    pos = node_positions(tree)
    for c in tree.find_clades():
        nm = getattr(c, "clade_name", None)
        if not nm or (names is not None and nm not in names):
            continue
        ys = [pos[t][1] for t in c.get_terminals()]
        y0, y1 = min(ys) - 0.5 + gap, max(ys) + 0.5 - gap
        ax.plot([x, x], [y0, y1], color="#444", lw=lw, transform=tr, clip_on=False)
        ax.text(x + pad, (y0 + y1) / 2, nm, fontsize=fontsize, va="center", ha="left",
                style="italic", color="#444", transform=tr, clip_on=False)
