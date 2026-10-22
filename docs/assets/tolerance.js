/* The tolerance recipe and the statistics used on the page, in plain JavaScript.
 *
 * tolerance() is a line-by-line port of figures/convolve.py (convolve_v4) and figures/tolerance.py
 * (tolerance_track), so the page can rebuild the score for any gene with different settings. With the default
 * settings it reproduces the official tolerance_v4 values (checked for every gene by tools/explorer/selftest.html).
 */
(function () {
  "use strict";
  const TOL = (window.TOL = window.TOL || {});

  // ---------------------------------------------------------------- substitution severity
  // BLOSUM62 as loaded by Biopython (substitution_matrices.load("BLOSUM62")); min -4, max 11 over the whole matrix
  const AA = "ARNDCQEGHILKMFPSTWYV";
  const BLOSUM62 = [
    [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0],
    [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],
    [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3],
    [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3],
    [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],
    [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2],
    [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2],
    [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3],
    [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3],
    [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3],
    [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1],
    [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2],
    [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1],
    [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1],
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2],
    [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2],
    [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0],
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3],
    [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1],
    [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4],
  ];
  const BLOSUM_LO = -4, BLOSUM_HI = 11;
  // Grantham's distance (0-215), upper triangle in this residue order, as in convolve.py
  const G_ORDER = "SRLPTAVGIFYCHQNKDEMW";
  const G_UPPER = [
    [0, 110, 145, 74, 58, 99, 124, 56, 142, 155, 144, 112, 89, 68, 46, 121, 65, 80, 135, 177],
    [0, 0, 102, 103, 71, 112, 96, 125, 97, 97, 77, 180, 29, 43, 86, 26, 96, 54, 91, 101],
    [0, 0, 0, 98, 92, 96, 32, 138, 5, 22, 36, 198, 99, 113, 153, 107, 172, 138, 15, 61],
    [0, 0, 0, 0, 38, 27, 68, 42, 95, 114, 110, 169, 77, 76, 91, 103, 108, 93, 87, 147],
    [0, 0, 0, 0, 0, 58, 69, 59, 89, 103, 92, 149, 47, 42, 65, 78, 85, 65, 81, 128],
    [0, 0, 0, 0, 0, 0, 64, 60, 94, 113, 112, 195, 86, 91, 111, 106, 126, 107, 84, 148],
    [0, 0, 0, 0, 0, 0, 0, 109, 29, 50, 55, 192, 84, 96, 133, 97, 152, 121, 21, 88],
    [0, 0, 0, 0, 0, 0, 0, 0, 135, 153, 147, 159, 98, 87, 80, 127, 94, 98, 127, 184],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 21, 33, 198, 94, 109, 149, 102, 168, 134, 10, 61],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 22, 205, 100, 116, 158, 102, 177, 140, 28, 40],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 194, 83, 99, 143, 85, 160, 122, 36, 37],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 174, 154, 139, 202, 154, 170, 196, 215],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 24, 68, 32, 81, 40, 87, 115],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 46, 53, 61, 29, 101, 130],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 94, 23, 42, 142, 174],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 101, 56, 95, 110],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 45, 160, 181],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 126, 152],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 67],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
  ];

  /** BLOSUM62 severity, 0 (similar) to 1 (dissimilar): (max - score) / (max - min). Unknown residues -> 1. */
  function blosumSeverity(ref, alt) {
    const i = AA.indexOf(ref), j = AA.indexOf(alt);
    if (i < 0 || j < 0) return 1.0;            // '*' scores -4 in Biopython -> (11 + 4) / 15 = 1
    return (BLOSUM_HI - BLOSUM62[i][j]) / (BLOSUM_HI - BLOSUM_LO);
  }

  /** Grantham severity: distance / 215 (0 for no change). Unknown residues -> 1. */
  function granthamSeverity(ref, alt) {
    const i = G_ORDER.indexOf(ref), j = G_ORDER.indexOf(alt);
    if (i < 0 || j < 0) return 1.0;
    if (i === j) return 0;
    return (G_UPPER[i][j] || G_UPPER[j][i]) / 215;
  }

  // ---------------------------------------------------------------- the smoothing
  const kernels = new Map();
  /** Binomial(2w, 0.5) probabilities: the smoothing window spans 2w + 1 residues. */
  function kernel(w) {
    if (!kernels.has(w)) {
      const n = 2 * w, k = new Float64Array(n + 1);
      k[0] = Math.pow(0.5, n);
      for (let r = 1; r <= n; r++) k[r] = (k[r - 1] * (n - r + 1)) / r;
      kernels.set(w, k);
    }
    return kernels.get(w);
  }

  /** Severity-weighted variant density along a protein (convolve_v4 with matrix_from="variant",
   *  edge="renormalize").
   *  variants: [{ref, pos, alt}], 1-based positions; disordered: boolean per residue (disorder >= threshold);
   *  wOrdered / wDisordered: window w where the residue is ordered / disordered;
   *  aggregate: how several variants at one residue combine: "max" (official), "sum", "mean" or "binary". */
  function density(variants, L, disordered, wOrdered, wDisordered, aggregate) {
    const vals = new Map();
    for (const v of variants) {
      const i = v.pos - 1;
      if (i < 0 || i >= L) continue;
      // the severity table is chosen by the variant's OWN residue: Grantham where disordered, BLOSUM62 elsewhere
      const s = disordered[i] ? granthamSeverity(v.ref, v.alt) : blosumSeverity(v.ref, v.alt);
      if (!vals.has(i)) vals.set(i, []);
      vals.get(i).push(s);
    }
    const weight = new Float64Array(L);
    for (const [i, list] of vals) {
      if (aggregate === "sum") weight[i] = list.reduce((a, b) => a + b, 0);
      else if (aggregate === "mean") weight[i] = list.reduce((a, b) => a + b, 0) / list.length;
      else if (aggregate === "binary") weight[i] = 1;
      else weight[i] = Math.max(...list);
    }
    const out = new Float64Array(L);
    for (let i = 0; i < L; i++) {
      const w = disordered[i] ? wDisordered : wOrdered, k = kernel(w);
      let raw = 0, kept = 0;
      for (let j = 0; j < k.length; j++) {
        const idx = i - w + j;
        if (idx >= 0 && idx < L) { raw += weight[idx] * k[j]; kept += k[j]; }
      }
      out[i] = kept ? raw / kept : raw;      // renormalise by the part of the window inside the protein
    }
    return { density: out, weight };
  }

  /** Min-max scale to [0, 1]; a flat track (e.g. no variants) becomes all zeros, as in tolerance.normalize. */
  function normalize(y) {
    let mn = Infinity, mx = -Infinity;
    for (const v of y) { if (v < mn) mn = v; if (v > mx) mx = v; }
    const out = new Float64Array(y.length);
    if (mx === mn) return out;
    for (let i = 0; i < y.length; i++) out[i] = (y[i] - mn) / (mx - mn);
    return out;
  }

  /** The tolerance score of one gene: normalised gnomAD density - normalised ClinVar density.
   *  opts: {wOrdered: 31, wDisordered: 5, threshold: 0.5, aggregate: "max"} (the official recipe). */
  function tolerance(gnomad, clinvar, disorder, opts) {
    const o = Object.assign({ wOrdered: 31, wDisordered: 5, threshold: 0.5, aggregate: "max" }, opts || {});
    const L = disorder.length, disordered = new Uint8Array(L);
    for (let i = 0; i < L; i++) disordered[i] = disorder[i] >= o.threshold ? 1 : 0;
    const g = density(gnomad, L, disordered, o.wOrdered, o.wDisordered, o.aggregate);
    const c = density(clinvar, L, disordered, o.wOrdered, o.wDisordered, o.aggregate);
    const gn = normalize(g.density), cn = normalize(c.density), tol = new Float64Array(L);
    for (let i = 0; i < L; i++) tol[i] = gn[i] - cn[i];
    return { tol, gnomad: gn, clinvar: cn, gnomadWeight: g.weight, clinvarWeight: c.weight, disordered };
  }

  // ---------------------------------------------------------------- statistics
  function median(a) {
    const s = Array.from(a).filter((x) => Number.isFinite(x)).sort((x, y) => x - y);
    if (!s.length) return NaN;
    const m = s.length >> 1;
    return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
  }

  /** Ranks 1..n with ties averaged. */
  function ranks(a) {
    const idx = Array.from(a, (_, i) => i).sort((i, j) => a[i] - a[j]);
    const r = new Float64Array(a.length);
    for (let s = 0; s < idx.length; ) {
      let e = s;
      while (e + 1 < idx.length && a[idx[e + 1]] === a[idx[s]]) e++;
      const avg = (s + e) / 2 + 1;
      for (let k = s; k <= e; k++) r[idx[k]] = avg;
      s = e + 1;
    }
    return r;
  }

  function pearson(x, y) {
    const n = x.length;
    let mx = 0, my = 0;
    for (let i = 0; i < n; i++) { mx += x[i]; my += y[i]; }
    mx /= n; my /= n;
    let sxy = 0, sxx = 0, syy = 0;
    for (let i = 0; i < n; i++) {
      const dx = x[i] - mx, dy = y[i] - my;
      sxy += dx * dy; sxx += dx * dx; syy += dy * dy;
    }
    return sxx && syy ? sxy / Math.sqrt(sxx * syy) : NaN;
  }

  const spearman = (x, y) => pearson(ranks(x), ranks(y));

  /** AUROC: the chance that a random "positive" scores higher than a random "negative" (ties count half).
   *  Here positives are benign variants and every score is "higher = more benign". */
  function auroc(pos, neg) {
    const n1 = pos.length, n0 = neg.length;
    if (!n1 || !n0) return NaN;
    const all = Float64Array.from([...pos, ...neg]), r = ranks(all);
    let sum = 0;
    for (let i = 0; i < n1; i++) sum += r[i];
    return (sum - (n1 * (n1 + 1)) / 2) / (n1 * n0);
  }

  /** AUROC from two histograms on the same bins: P(a > b) + half the ties within a bin. */
  function binnedAuroc(a, b) {
    const na = a.reduce((s, x) => s + x, 0), nb = b.reduce((s, x) => s + x, 0);
    if (!na || !nb) return NaN;
    let below = 0, acc = 0;
    for (let i = 0; i < a.length; i++) {
      acc += a[i] * (below + 0.5 * b[i]);
      below += b[i];
    }
    return acc / (na * nb);
  }

  /** ROC curve for "higher score = benign": points (share of pathogenic called benign, share of benign called
   *  benign) as the cut-off sweeps from the top score down. */
  function rocCurve(pos, neg) {
    const items = pos.map((s) => [s, 1]).concat(neg.map((s) => [s, 0])).sort((a, b) => b[0] - a[0]);
    const pts = [[0, 0]];
    let tp = 0, fp = 0;
    for (let i = 0; i < items.length; ) {
      const s = items[i][0];
      while (i < items.length && items[i][0] === s) { if (items[i][1]) tp++; else fp++; i++; }
      pts.push([fp / neg.length, tp / pos.length]);
    }
    return pts;
  }

  TOL.recipe = { AA, blosumSeverity, granthamSeverity, kernel, density, normalize, tolerance,
                 DEFAULTS: { wOrdered: 31, wDisordered: 5, threshold: 0.5, aggregate: "max",
                             gnomad: "ac2", clinvar: "all" } };
  TOL.stats = { median, ranks, pearson, spearman, auroc, binnedAuroc, rocCurve };
})();
