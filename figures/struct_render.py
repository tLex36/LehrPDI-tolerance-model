"""Protein structures drawn with matplotlib: shaded cartoon ribbons and spheres, depth-sorted (painter's
algorithm), orthographic view. Everything else in the figure (labels, arrows, legends) stays ordinary editable
matplotlib text; the molecule itself can be rasterised so the SVG stays light.

    from struct_render import chain_residues, secondary_structure, view_matrix, Scene
    res = chain_residues("AF-P23415-2-F1.pdb", "A")
    ss = secondary_structure(res)                       # 'H' helix, 'E' strand, 'C' coil (simplified DSSP)
    sc = Scene(view_matrix(up=[0, 0, 1], front=[1, 0, 0]), center=res.ca.mean(0))
    sc.cartoon(res, ss, colors)                         # one colour per residue
    sc.spheres(xyz, radii, colors)
    sc.draw(ax)
"""
import numpy as np
from dataclasses import dataclass
from matplotlib.collections import PolyCollection
from matplotlib.colors import to_rgb

THREE_TO_ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
                "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
                "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M"}
VDW = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8}


@dataclass
class Residues:
    """Backbone and side-chain coordinates of one chain (only residues with N, CA, C and O)."""
    num: np.ndarray          # residue numbers as in the file (author numbering)
    aa: str                  # one-letter sequence
    n: np.ndarray            # (L, 3) backbone atoms
    ca: np.ndarray
    c: np.ndarray
    o: np.ndarray
    side: list               # per residue: (k, 3) side-chain heavy atoms (CB onwards; empty for Gly)
    side_el: list            # per residue: element letters of those atoms

    def __len__(self):
        return len(self.num)

    def centroid(self):
        """Side-chain centroid per residue (CA for glycine)."""
        return np.array([s.mean(0) if len(s) else a for s, a in zip(self.side, self.ca)])


def read_structure(path):
    from Bio.PDB import PDBParser, MMCIFParser
    parser = MMCIFParser(QUIET=True) if str(path).endswith((".cif", ".mmcif")) else PDBParser(QUIET=True)
    return parser.get_structure("s", str(path))


def chain_residues(path_or_structure, chain):
    """Residues of one chain of the first model, standard amino acids only."""
    s = read_structure(path_or_structure) if isinstance(path_or_structure, str) else path_or_structure
    num, aa, n, ca, c, o, side, side_el = [], [], [], [], [], [], [], []
    for r in s[0][chain]:
        if r.id[0] != " " or r.get_resname() not in THREE_TO_ONE or not all(a in r for a in ("N", "CA", "C", "O")):
            continue
        num.append(r.id[1]); aa.append(THREE_TO_ONE[r.get_resname()])
        n.append(r["N"].coord); ca.append(r["CA"].coord); c.append(r["C"].coord); o.append(r["O"].coord)
        sc = [a for a in r if a.get_id() not in ("N", "CA", "C", "O", "OXT") and a.element != "H"]
        side.append(np.array([a.coord for a in sc], float).reshape(-1, 3)); side_el.append([a.element for a in sc])
    f = lambda v: np.asarray(v, float)
    return Residues(np.asarray(num), "".join(aa), f(n), f(ca), f(c), f(o), side, side_el)


def secondary_structure(res, min_helix=4, min_strand=2):
    """'H' / 'E' / 'C' per residue from backbone hydrogen bonds (Kabsch & Sander DSSP, simplified:
    alpha helices from two consecutive i->i+4 turns, strands from ladders of beta bridges)."""
    L = len(res)
    N, C, O = res.n, res.c, res.o
    # amide H placed 1 A from N, opposite the previous carbonyl (no H on the first residue, on prolines, or after gaps)
    H = np.full_like(N, np.nan)
    d = C[:-1] - O[:-1]
    H[1:] = N[1:] + d / np.linalg.norm(d, axis=1, keepdims=True)
    gap = np.r_[True, np.linalg.norm(N[1:] - C[:-1], axis=1) > 2.0]
    H[gap | np.array([a == "P" for a in res.aa])] = np.nan
    dist = lambda A, B: np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)
    with np.errstate(invalid="ignore", divide="ignore"):
        E = 0.084 * 332 * (1 / dist(O, N) + 1 / dist(C, H) - 1 / dist(O, H) - 1 / dist(C, N))   # E[i, j]: CO(i) -> NH(j)
    hb = np.nan_to_num(E, nan=0.0) < -0.5
    idx = np.arange(L)
    hb[np.abs(idx[:, None] - idx[None, :]) < 2] = False
    ss = np.array(["C"] * L)
    turn4 = np.array([i + 4 < L and hb[i, i + 4] for i in range(L)])
    for i in range(1, L - 4):
        if turn4[i - 1] and turn4[i]:
            ss[i:i + 4] = "H"
    hbf = lambda i, j: 0 <= i < L and 0 <= j < L and hb[i, j]
    bridge = np.zeros(L, bool)
    for i in range(1, L - 1):
        for j in range(L):
            if abs(i - j) < 3 or not (0 < j < L - 1):
                continue
            par = (hbf(i - 1, j) and hbf(j, i + 1)) or (hbf(j - 1, i) and hbf(i, j + 1))
            anti = (hbf(i, j) and hbf(j, i)) or (hbf(i - 1, j + 1) and hbf(j - 1, i + 1))
            if par or anti:
                bridge[i] = True
                break
    for i in range(L):                                            # ladders: a bridge next to another bridge
        if bridge[i] and ss[i] != "H" and ((i > 0 and bridge[i - 1]) or (i + 1 < L and bridge[i + 1])):
            ss[i] = "E"
    for kind, shortest in (("H", min_helix), ("E", min_strand)):  # drop fragments that are too short
        i = 0
        while i < L:
            if ss[i] == kind:
                j = i
                while j < L and ss[j] == kind:
                    j += 1
                if j - i < shortest:
                    ss[i:j] = "C"
                i = j
            else:
                i += 1
    return "".join(ss)


def view_matrix(up, front):
    """Rotation whose rows are screen-right, screen-up and toward-the-viewer, from two directions in the
    structure's frame: `up` (made exactly vertical) and `front` (made to point at the viewer)."""
    u = np.asarray(up, float); u /= np.linalg.norm(u)
    f = np.asarray(front, float); f = f - (f @ u) * u; f /= np.linalg.norm(f)
    r = np.cross(u, f)
    return np.vstack([r, u, f])


def rotate_about(R, axis_screen, degrees):
    """Turn the view: rotate by `degrees` about a screen axis ('x', 'y' or 'z')."""
    t = np.radians(degrees); c, s = np.cos(t), np.sin(t)
    M = {"x": [[1, 0, 0], [0, c, -s], [0, s, c]], "y": [[c, 0, s], [0, 1, 0], [-s, 0, c]],
         "z": [[c, -s, 0], [s, c, 0], [0, 0, 1]]}[axis_screen]
    return np.asarray(M) @ R


def _catmull_rom(P, k):
    """k samples per interval through points P (m, d); returns points and tangents ((m-1)*k + 1, d)."""
    P = np.asarray(P, float)
    ext = np.vstack([2 * P[0] - P[1], P, 2 * P[-1] - P[-2]])
    t = np.linspace(0, 1, k, endpoint=False)[:, None]
    pts, tan = [], []
    for i in range(len(P) - 1):
        p0, p1, p2, p3 = ext[i], ext[i + 1], ext[i + 2], ext[i + 3]
        a, b, c3 = -p0 + p2, 2 * p0 - 5 * p1 + 4 * p2 - p3, -p0 + 3 * p1 - 3 * p2 + p3
        pts.append(0.5 * (2 * p1 + a * t + b * t ** 2 + c3 * t ** 3))
        tan.append(0.5 * (a + 2 * b * t + 3 * c3 * t ** 2))
    pts.append(P[-1][None]); tan.append((0.5 * (ext[-1] - ext[-3]))[None])
    return np.vstack(pts), np.vstack(tan)


class Scene:
    """Collects polygons in view space, then draws them back to front on a matplotlib axis."""

    def __init__(self, rotation, center):
        self.R = np.asarray(rotation, float)
        self.c0 = np.asarray(center, float)
        self.polys, self.depth, self.rgb, self.edge, self.kind = [], [], [], [], []

    def to_view(self, xyz):
        return (np.asarray(xyz, float) - self.c0) @ self.R.T

    def project(self, xyz):
        """Screen (x, y) of structure coordinates: use it to place labels and arrows."""
        return self.to_view(np.atleast_2d(xyz))[:, :2]

    # ------------------------------------------------------------------ cartoon
    def cartoon(self, res, ss, colors, helix=1.35, strand=1.15, arrow=1.9, coil=0.32, samples=8, keep=None, scale=None):
        """Ribbons for helices and strands, a thin tube for coil. colors: one colour per residue.
        keep: optional boolean mask of residues to draw (the chain is still split at gaps).
        scale: optional per-residue width factor, e.g. 2.5 to make annotated stretches stand out."""
        keep = np.ones(len(res), bool) if keep is None else np.asarray(keep, bool)
        scale = np.ones(len(res)) if scale is None else np.asarray(scale, float)
        ca = self.to_view(res.ca)
        ori = self.to_view(res.o) - self.to_view(res.c)                       # carbonyl directions
        ori /= np.linalg.norm(ori, axis=1, keepdims=True)
        brk = np.r_[np.linalg.norm(np.diff(res.ca, axis=0), axis=1) > 4.3, True]
        cols = [np.array(to_rgb(c)) for c in colors]
        start = None
        for i in range(len(res)):
            if keep[i] and start is None:
                start = i
            last = start is not None and (brk[i] or i + 1 == len(res) or not keep[i + 1])
            if last:
                if i - start >= 1:
                    self._segment(ca[start:i + 1], ori[start:i + 1], ss[start:i + 1], cols[start:i + 1],
                                  helix, strand, arrow, coil, samples, scale[start:i + 1])
                start = None

    def _segment(self, P, O, ss, cols, helix, strand, arrow, coil, k, scale):
        m = len(P)
        O = O.copy()
        for i in range(1, m):                                                  # keep the ribbon from flipping
            if O[i] @ O[i - 1] < 0:
                O[i] = -O[i]
        Os = O.copy()
        for i in range(m):                                                     # smooth strand pleating
            if ss[i] == "E":
                lo, hi = max(0, i - 1), min(m, i + 2)
                Os[i] = O[lo:hi].mean(0)
        pts, tan = _catmull_rom(P, k)
        tan /= np.linalg.norm(tan, axis=1, keepdims=True)
        n_s = len(pts)
        res_of = np.minimum((np.arange(n_s) + k // 2) // k, m - 1)             # nearest residue of each sample
        base = {"H": helix, "E": strand, "C": coil}
        half = np.empty(n_s); flat = np.empty(n_s, bool)
        for s in range(n_s):
            i, t = s // k, (s % k) / k
            j = min(i + 1, m - 1)
            a, b = base[ss[i]], base[ss[j]]
            if ss[i] == "E" and (j == m - 1 or ss[j] != "E") and j != i:       # arrowhead after a strand's last residue
                a, b = arrow, coil
            half[s] = (a + (b - a) * t) * ((1 - t) * scale[i] + t * scale[j])
            flat[s] = ss[res_of[s]] != "C" or (ss[i] == "E" and ss[j] != "E")
        W = np.empty_like(pts)
        for s in range(n_s):
            i, t = s // k, (s % k) / k
            j = min(i + 1, m - 1)
            o = (1 - t) * Os[i] + t * Os[j]
            w = o - (o @ tan[s]) * tan[s]
            cam = np.cross(tan[s], [0.0, 0.0, 1.0])
            if not flat[s] or np.linalg.norm(w) < 1e-6:
                w = cam
            W[s] = w / (np.linalg.norm(w) or 1.0)
        for s in range(n_s - 1):
            p, q = pts[s], pts[s + 1]
            quad = np.array([p + half[s] * W[s], q + half[s + 1] * W[s + 1], q - half[s + 1] * W[s + 1], p - half[s] * W[s]])
            nrm = np.cross(q - p, W[s]); nn = np.linalg.norm(nrm)
            tube = not (flat[s] and flat[s + 1])
            self._add(quad, quad[:, 2].mean(), cols[res_of[s]], None if tube else (nrm / nn if nn else None), "ribbon")

    # ------------------------------------------------------------------ spheres
    def spheres(self, xyz, radii, colors, n=28):
        """Shaded spheres (a disc plus two highlights each)."""
        V = self.to_view(xyz)
        radii = np.broadcast_to(np.asarray(radii, float), (len(V),))
        ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
        circ = np.c_[np.cos(ang), np.sin(ang)]
        for (x, y, z), r, col in zip(V, radii, colors):
            col = np.array(to_rgb(col))
            self._add(np.c_[x + r * circ[:, 0], y + r * circ[:, 1], np.full(n, z)], z + r * 0.0, col, "sphere", "sphere")
            for off, rr, mix in ((0.28, 0.62, 0.22), (0.38, 0.28, 0.5)):
                cx, cy = x - off * r, y + off * r
                self._add(np.c_[cx + rr * r * circ[:, 0], cy + rr * r * circ[:, 1], np.full(n, z)],
                          z + 1e-3 * (1 + mix), col + (1 - col) * mix, "flat", "highlight")

    # ------------------------------------------------------------------ drawing
    def _add(self, poly3, depth, rgb, normal, kind):
        self.polys.append(np.asarray(poly3)[:, :2]); self.depth.append(depth)
        self.rgb.append((np.asarray(rgb, float), normal)); self.kind.append(kind)

    def draw(self, ax, light=(-0.4, 0.55, 1.0), ambient=0.42, fog=0.38, edge_lw=0.35, rasterized=True, zorder=1,
             adjustable="datalim"):
        """Add every polygon to ax, back to front. fog: how much the farthest parts fade toward white.
        adjustable="datalim" keeps the axis box at its grid size and widens the data range instead (equal aspect)."""
        if not self.polys:
            return None
        L = np.asarray(light, float); L /= np.linalg.norm(L)
        d = np.asarray(self.depth)
        far, near = d.min(), d.max()
        order = np.argsort(d, kind="stable")
        faces, edges, lws = [], [], []
        for i in order:
            rgb, normal = self.rgb[i]
            kind = self.kind[i]
            if normal is None:
                shade = 0.88
            elif isinstance(normal, str):                                       # sphere disc
                shade = 0.78
            else:
                shade = ambient + (1 - ambient) * abs(normal @ L)
            col = np.clip(rgb * shade, 0, 1) if kind != "highlight" else np.clip(rgb, 0, 1)
            fade = fog * (near - d[i]) / ((near - far) or 1.0)
            col = col + (1 - col) * fade
            faces.append(col)
            if kind == "sphere":
                edges.append(np.clip(col * 0.62, 0, 1)); lws.append(0.45)
            elif kind == "highlight":
                edges.append(col); lws.append(0.0)
            else:
                edges.append(col); lws.append(edge_lw)
        pc = PolyCollection([self.polys[i] for i in order], facecolors=faces, edgecolors=edges, linewidths=lws,
                            joinstyle="round", zorder=zorder)
        pc.set_rasterized(rasterized)
        ax.add_collection(pc)
        allp = np.vstack(self.polys)
        ax.update_datalim(allp)
        ax.autoscale_view()
        ax.set_aspect("equal", adjustable=adjustable)
        return pc
