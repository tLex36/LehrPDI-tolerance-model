/* "Explore a gene": tracks along one protein, with the tolerance recipe rebuilt in the browser. */
(function () {
  "use strict";
  const TOL = window.TOL;
  const { el, fmt, color, showTip, hideTip, widthOf, segmented, hashParams } = TOL.ui;
  const C = TOL.charts;
  const R = TOL.recipe, S = TOL.stats;

  const EXAMPLES = [["GLRA1", "poster example"], ["GRIN2B", "many labelled variants"], ["SHANK3", "long, largely disordered"],
                    ["NLGN4X", ""], ["RYR1", "test gene"], ["VHL", "test gene"]];
  const CLASS_NAME = { p: "pathogenic", b: "benign", v: "VUS", o: "other ClinVar" };
  const CLASS_VAR = { p: "--pathogenic", b: "--benign", v: "--vus", o: "--other" };
  const AA3 = { A: "Ala", R: "Arg", N: "Asn", D: "Asp", C: "Cys", Q: "Gln", E: "Glu", G: "Gly", H: "His", I: "Ile", L: "Leu",
                K: "Lys", M: "Met", F: "Phe", P: "Pro", S: "Ser", T: "Thr", W: "Trp", Y: "Tyr", V: "Val" };
  const ANN_ROWS = [
    { label: "Domain / repeat", keys: [["domain", "domain (InterPro)", "--accent"], ["repeat", "repeat", "--boundary"]] },
    { label: "Membrane / signal", keys: [["tm", "transmembrane segment", "--tol"], ["signal", "signal peptide", "--other"]] },
    { label: "Interface", keys: [["iface", "touches another chain (PDBe-KB)", "--vus"]] },
    { label: "Sites", keys: [["ligand", "ligand site", "--pathogenic"], ["active", "active site", "--pathogenic"],
                              ["disulfide", "disulfide", "--disorder-line"], ["ptm", "PTM site", "--esm"]] },
    { label: "Motifs / binding", keys: [["motif", "curated motif (ELM)", "--benign"], ["motifPred", "predicted motif (ELM)", "--benign"],
                                         ["lip", "curated binding peptide (MobiDB)", "--benign"], ["lipPred", "predicted binding peptide (MobiDB)", "--benign"],
                                         ["bindingRegion", "binding region (UniProt)", "--benign"]] },
  ];

  const state = {
    gene: null, data: null, info: null, opts: { ...R.DEFAULTS }, view: null, official: null, current: null,
    show: { gnomad: true, p: true, b: true, v: true, o: true }, marker: null, hover: null,
  };

  // ---------------------------------------------------------------- the recipe for the current settings
  function variantsFor(opts) {
    const d = state.data;
    const gn = d.gnomad.filter((v) => v.flag === "u" || (opts.gnomad === "ac1" && v.flag === "s"));
    const keep = opts.clinvar === "pathogenic" ? ["p"] : opts.clinvar === "noBenign" ? ["p", "v", "o"] : ["p", "b", "v", "o"];
    const cv = d.clinvar.filter((v) => keep.includes(v.flag));
    return { gn, cv };
  }
  const isOfficial = (o) => Object.keys(R.DEFAULTS).every((k) => o[k] === R.DEFAULTS[k]);
  function compute() {
    const d = state.data;
    const run = (opts) => { const { gn, cv } = variantsFor(opts); return { ...R.tolerance(gn, cv, d.dis, opts), gn, cv }; };
    state.official = run(R.DEFAULTS);
    state.current = isOfficial(state.opts) ? state.official : run(state.opts);
  }

  // ---------------------------------------------------------------- gene loading
  async function loadGene(name) {
    const idx = await TOL.data.index();
    const info = idx.byName.get(name);
    if (!info) { setStatus(`No gene called "${name}" in this set.`); return; }
    const card = document.getElementById("explorer-card");
    card.classList.add("loading");
    try {
      state.data = await TOL.data.gene(name);
    } catch (e) {
      setStatus(e.message); card.classList.remove("loading"); return;
    }
    state.gene = name; state.info = info; state.view = [1, state.data.L]; state.marker = null;
    document.getElementById("gene-input").value = name;
    compute();
    header();
    draw();
    variantPanel();
    recipeStatus();
    card.classList.remove("loading");
    try { history.replaceState(null, "", "#gene=" + encodeURIComponent(name)); } catch (e) { /* file:// in some browsers */ }
  }
  function setStatus(t) { document.getElementById("recipe-status").textContent = t; }

  function header() {
    const d = state.data, i = state.info;
    const SPLIT = { train: "training gene", val: "validation gene", test: "untouched test gene", excluded: "excluded (no usable score)" };
    const cvCounts = { p: 0, b: 0, v: 0, o: 0 };
    d.clinvar.forEach((v) => cvCounts[v.flag]++);
    const uniprot = d.uniprot ? el("a", { href: `https://www.uniprot.org/uniprotkb/${d.uniprot}/entry`, text: d.uniprot, target: "_blank", rel: "noopener" }) : null;
    document.getElementById("gene-head").replaceChildren(
      el("span", { class: "name", text: d.gene }),
      el("span", { class: "badge " + d.split, text: SPLIT[d.split] || d.split }),
      el("span", { class: "facts" }, `${fmt.int(d.L)} residues`, uniprot ? " · UniProt " : "", uniprot || "",
        ` · gnomAD: ${fmt.int(i.gnomad_used)} variants used · ClinVar: ${fmt.int(cvCounts.p)} pathogenic, ${fmt.int(cvCounts.b)} benign, `
        + `${fmt.int(cvCounts.v)} VUS, ${fmt.int(cvCounts.o)} other · labelled (≥ 1 star): ${i.labelled_pathogenic} pathogenic, ${i.labelled_benign} benign`));
  }

  function legend() {
    const host = document.getElementById("explorer-legend");
    const items = [["gnomad", "gnomAD (healthy population)", "--gnomad"], ["p", "ClinVar pathogenic", "--pathogenic"],
                   ["b", "ClinVar benign", "--benign"], ["v", "ClinVar VUS", "--vus"], ["o", "ClinVar other", "--other"]];
    host.replaceChildren(el("span", { class: "note", text: "Show:" }), ...items.map(([k, label, c]) => {
      const b = el("button", { class: "item" + (state.show[k] ? "" : " off"), type: "button", "aria-pressed": String(state.show[k]) },
        el("span", { class: "swatch box", style: `background:${color(c)}` }), label);
      b.addEventListener("click", () => { state.show[k] = !state.show[k]; legend(); draw(); });
      return b;
    }), el("span", { class: "item" }, el("span", { class: "swatch box", style: `background:${color("--disorder-band")};outline:1px solid ${color("--axis")}` }), "disordered region"));
  }

  // ---------------------------------------------------------------- drawing
  const TRACKS = [
    { key: "variants", label: ["Variants", "length = severity"], h: 112 },
    { key: "density", label: ["Smoothed", "density (0–1)"], h: 62 },
    { key: "tol", label: ["Tolerance", "gnomAD − ClinVar"], h: 100 },
    { key: "esm", label: ["Evolution", "ESM-1v mean LLR"], h: 52 },
    { key: "dis", label: ["Disorder", "AIUPred"], h: 44 },
    { key: "rsa", label: ["Burial", "RSA"], h: 44 },
    { key: "plddt", label: ["AlphaFold", "pLDDT"], h: 36 },
    { key: "ann", label: [], h: ANN_ROWS.length * 13 },
  ];
  const GAP = 12, LEFT = 112, RIGHT = 14, TOP = 6;

  function runsOf(mask, lo, hi) {
    const out = [];
    let s = null;
    for (let i = lo; i <= hi; i++) {
      const on = mask[i - 1];
      if (on && s === null) s = i;
      if ((!on || i === hi) && s !== null) { out.push([s, on && i === hi ? i : i - 1]); s = null; }
    }
    return out;
  }

  function draw() {
    if (!state.data) return;
    const host = document.getElementById("explorer-chart"), d = state.data, cur = state.current, off = state.official;
    const W = widthOf(host);
    let y = TOP;
    const pos = {};
    TRACKS.forEach((t) => { pos[t.key] = y; y += t.h + GAP; });
    const axisY = y - GAP + 4, ovY = axisY + 34, ovH = 34, H = ovY + ovH + 22;
    const svg = C.svg(host, W, H, `Tracks along ${d.gene}`);
    const [v0, v1] = state.view;
    const x = d3.scaleLinear().domain([v0 - 0.5, v1 + 0.5]).range([LEFT, W - RIGHT]);
    const step = (x(v0 + 1) - x(v0));
    const inView = (p) => p >= v0 && p <= v1;
    const clipId = "clip-explorer";
    svg.append("defs").append("clipPath").attr("id", clipId).append("rect").attr("x", LEFT).attr("y", 0).attr("width", W - LEFT - RIGHT).attr("height", axisY);
    const plot = svg.append("g").attr("clip-path", `url(#${clipId})`);

    // disordered regions, behind every track
    for (const [s, e] of runsOf(cur.disordered, v0, v1)) {
      plot.append("rect").attr("x", x(s - 0.5)).attr("width", Math.max(1, x(e + 0.5) - x(s - 0.5))).attr("y", TOP - 2)
        .attr("height", axisY - TOP - 4).attr("fill", color("--disorder-band"));
    }
    // track labels and baselines
    TRACKS.forEach((t) => {
      const ty = pos[t.key];
      const g = svg.append("g").attr("transform", `translate(0,${ty})`);
      if (t.label[0]) g.append("text").attr("class", "lab-strong").attr("x", 0).attr("y", 12).text(t.label[0]);
      if (t.label[1]) g.append("text").attr("class", "lab").attr("x", 0).attr("y", 27).style("font-size", "11px").text(t.label[1]);
    });
    const emptyNote = (yy, text) => svg.append("text").attr("class", "lab").attr("x", (LEFT + W - RIGHT) / 2).attr("y", yy)
      .attr("text-anchor", "middle").style("font-size", "11px").style("font-style", "italic").text(text);
    const line = (yScale, arr) => d3.line().defined((i) => Number.isFinite(arr[i - 1])).x((i) => x(i)).y((i) => yScale(arr[i - 1]));
    const idx = d3.range(v0, v1 + 1);

    // 1. variants: gnomAD up, ClinVar down; length = the severity weight the recipe gives the swap
    {
      const ty = pos.variants, mid = ty + TRACKS[0].h / 2, half = TRACKS[0].h / 2 - 5;
      plot.append("line").attr("x1", LEFT).attr("x2", W - RIGHT).attr("y1", mid).attr("y2", mid).attr("stroke", color("--axis"));
      const sev = (v) => (cur.disordered[v.pos - 1] ? R.granthamSeverity(v.ref, v.alt) : R.blosumSeverity(v.ref, v.alt));
      const r = step > 6 ? 2.6 : step > 2 ? 2 : 1.4;
      const draw1 = (list, up, col) => {                    // one path for the stems and one for the heads: fast for long proteins
        let stems = "", heads = "";
        for (const v of list) {
          if (!inView(v.pos)) continue;
          const len = 4 + sev(v) * (half - 4), yy = up ? mid - len : mid + len, xx = x(v.pos);
          stems += `M${xx.toFixed(1)},${mid}V${yy.toFixed(1)}`;
          heads += `M${(xx - r).toFixed(1)},${yy.toFixed(1)}a${r},${r} 0 1,0 ${2 * r},0a${r},${r} 0 1,0 ${-2 * r},0`;
        }
        if (!stems) return;
        plot.append("path").attr("d", stems).attr("stroke", col).attr("stroke-width", 1).attr("stroke-opacity", 0.75).attr("fill", "none");
        plot.append("path").attr("d", heads).attr("fill", col);
      };
      if (state.show.gnomad) draw1(cur.gn, true, color("--gnomad"));
      for (const k of ["o", "v", "b", "p"]) if (state.show[k]) draw1(cur.cv.filter((v) => v.flag === k), false, color(CLASS_VAR[k]));
      svg.append("text").attr("class", "lab").attr("x", W - RIGHT).attr("y", ty + 9).attr("text-anchor", "end").style("font-size", "10.5px").text("gnomAD ↑");
      svg.append("text").attr("class", "lab").attr("x", W - RIGHT).attr("y", ty + TRACKS[0].h - 1).attr("text-anchor", "end").style("font-size", "10.5px").text("ClinVar ↓");
    }
    // 2. smoothed densities
    {
      const ty = pos.density, yS = d3.scaleLinear().domain([0, 1]).range([ty + TRACKS[1].h, ty + 2]);
      plot.append("path").attr("d", line(yS, cur.gnomad)(idx)).attr("fill", "none").attr("stroke", color("--gnomad")).attr("stroke-width", 2);
      plot.append("path").attr("d", line(yS, cur.clinvar)(idx)).attr("fill", "none").attr("stroke", color("--pathogenic")).attr("stroke-width", 2);
      svg.append("line").attr("x1", LEFT).attr("x2", W - RIGHT).attr("y1", yS(0)).attr("y2", yS(0)).attr("stroke", color("--grid"));
    }
    // 3. tolerance: blue above 0 (tolerated), red below (intolerant); the official score in grey when the recipe differs
    {
      const ty = pos.tol, yS = d3.scaleLinear().domain([-1, 1]).range([ty + TRACKS[2].h, ty]);
      const area = (sign) => d3.area().x((i) => x(i)).y0(yS(0)).y1((i) => yS(sign > 0 ? Math.max(0, cur.tol[i - 1]) : Math.min(0, cur.tol[i - 1])));
      plot.append("path").attr("d", area(1)(idx)).attr("fill", color("--tolerated")).attr("fill-opacity", 0.28);
      plot.append("path").attr("d", area(-1)(idx)).attr("fill", color("--intolerant")).attr("fill-opacity", 0.28);
      if (cur !== off) plot.append("path").attr("d", line(yS, off.tol)(idx)).attr("fill", "none").attr("stroke", color("--muted")).attr("stroke-width", 1.2);
      plot.append("path").attr("d", line(yS, cur.tol)(idx)).attr("fill", "none").attr("stroke", color("--ink")).attr("stroke-width", 1.4);
      svg.append("line").attr("x1", LEFT).attr("x2", W - RIGHT).attr("y1", yS(0)).attr("y2", yS(0)).attr("stroke", color("--axis"));
      C.axis(svg.append("g").attr("transform", `translate(${LEFT},0)`), d3.axisLeft(yS).tickValues([-1, 0, 1]).tickSize(3), { domain: false });
      svg.append("text").attr("class", "lab").attr("x", W - RIGHT).attr("y", ty + 10).attr("text-anchor", "end").style("font-size", "10.5px").text("tolerated ↑");
      svg.append("text").attr("class", "lab").attr("x", W - RIGHT).attr("y", ty + TRACKS[2].h - 2).attr("text-anchor", "end").style("font-size", "10.5px").text("intolerant ↓");
    }
    // 4. ESM mean LLR (lower = evolution allows fewer swaps here)
    {
      const ty = pos.esm, ext = d3.extent(d.esm.filter(Number.isFinite));
      const yS = d3.scaleLinear().domain([Math.min(-1, ext[0] || -10), 0]).range([ty + TRACKS[3].h, ty + 2]);
      plot.append("path").attr("d", d3.area().defined((i) => Number.isFinite(d.esm[i - 1])).x((i) => x(i)).y0(yS(0)).y1((i) => yS(d.esm[i - 1]))(idx))
        .attr("fill", color("--esm")).attr("fill-opacity", 0.22);
      plot.append("path").attr("d", line(yS, d.esm)(idx)).attr("fill", "none").attr("stroke", color("--esm")).attr("stroke-width", 1.2);
      svg.append("text").attr("class", "lab").attr("x", W - RIGHT).attr("y", ty + 9).attr("text-anchor", "end").style("font-size", "10.5px").text("lower = fewer swaps allowed");
    }
    // 5. disorder with the current threshold
    {
      const ty = pos.dis, yS = d3.scaleLinear().domain([0, 1]).range([ty + TRACKS[4].h, ty + 2]);
      plot.append("line").attr("x1", LEFT).attr("x2", W - RIGHT).attr("y1", yS(state.opts.threshold)).attr("y2", yS(state.opts.threshold)).attr("stroke", color("--axis"));
      plot.append("path").attr("d", line(yS, d.dis)(idx)).attr("fill", "none").attr("stroke", color("--disorder-line")).attr("stroke-width", 1.4);
      plot.append("path").attr("d", line(yS, d.bind)(idx)).attr("fill", "none").attr("stroke", color("--vus")).attr("stroke-width", 1).attr("stroke-opacity", 0.8);
      svg.append("text").attr("class", "lab").attr("x", W - RIGHT).attr("y", ty + 9).attr("text-anchor", "end").style("font-size", "10.5px").text("disorder · binding");
    }
    // 6. RSA bars, coloured by class
    {
      const ty = pos.rsa, h = TRACKS[5].h, yS = d3.scaleLinear().domain([0, 1]).range([ty + h, ty + 2]);
      const bw = Math.max(1, step - (step > 4 ? 1 : 0)), paths = { "--core": "", "--boundary": "", "--surface-res": "" };
      for (let i = v0; i <= v1; i++) {
        const v = d.rsa[i - 1];
        if (!Number.isFinite(v)) continue;
        const c = v < 0.2 ? "--core" : v <= 0.4 ? "--boundary" : "--surface-res", top = yS(Math.min(1, v));
        paths[c] += `M${(x(i) - bw / 2).toFixed(1)},${yS(0).toFixed(1)}h${bw.toFixed(2)}V${Math.min(top, yS(0) - 0.5).toFixed(1)}h${(-bw).toFixed(2)}Z`;
      }
      for (const [c, dpath] of Object.entries(paths)) if (dpath) plot.append("path").attr("d", dpath).attr("fill", color(c));
      if (!d.rsa.some(Number.isFinite)) emptyNote(ty + h / 2 + 4, "no structure-based values were computed for this protein");
      svg.append("text").attr("class", "lab").attr("x", W - RIGHT).attr("y", ty + 9).attr("text-anchor", "end").style("font-size", "10.5px").text("core < 0.2 · surface > 0.4 · blank = no structure");
    }
    // 7. pLDDT
    {
      const ty = pos.plddt, yS = d3.scaleLinear().domain([0, 100]).range([ty + TRACKS[6].h, ty + 2]);
      plot.append("line").attr("x1", LEFT).attr("x2", W - RIGHT).attr("y1", yS(70)).attr("y2", yS(70)).attr("stroke", color("--grid"));
      plot.append("path").attr("d", line(yS, d.plddt)(idx)).attr("fill", "none").attr("stroke", color("--accent")).attr("stroke-width", 1.2);
      if (!d.plddt.some(Number.isFinite)) emptyNote(ty + TRACKS[6].h / 2 + 4, "no AlphaFold model was used for this protein");
    }
    // 8. annotations
    {
      const ty = pos.ann;
      ANN_ROWS.forEach((row, ri) => {
        const yy = ty + ri * 13;
        svg.append("text").attr("class", row.label ? "lab" : "lab").attr("x", 0).attr("y", yy + 9).style("font-size", "11px").text(row.label);
        for (const [key, , c] of row.keys) {
          for (const [s, e] of d.ann[key] || []) {
            if (e < v0 || s > v1) continue;
            plot.append("rect").attr("x", x(Math.max(s, v0) - 0.5)).attr("width", Math.max(1.5, x(Math.min(e, v1) + 0.5) - x(Math.max(s, v0) - 0.5)))
              .attr("y", yy + 2).attr("height", 8).attr("rx", 1.5).attr("fill", color(c)).attr("fill-opacity", key.endsWith("Pred") ? 0.5 : 0.9);
          }
        }
      });
    }
    // marker for a selected variant
    if (state.marker && inView(state.marker.pos)) {
      const mx = x(state.marker.pos);
      plot.append("line").attr("x1", mx).attr("x2", mx).attr("y1", TOP).attr("y2", axisY - 4).attr("stroke", color("--ink")).attr("stroke-width", 1.5);
      svg.append("text").attr("class", "lab-strong").attr("x", Math.min(mx + 4, W - RIGHT - 60)).attr("y", TOP + 10).text(state.marker.label);
    }
    // x axis
    C.axis(svg.append("g").attr("transform", `translate(0,${axisY})`), d3.axisBottom(x).ticks(Math.min(10, Math.round(W / 90))).tickFormat(d3.format("d")).tickSize(4));
    svg.append("text").attr("class", "lab").attr("x", W - RIGHT).attr("y", axisY + 30).attr("text-anchor", "end").text(`${d.gene} residue`);

    // crosshair layer
    const cross = svg.append("line").attr("y1", TOP).attr("y2", axisY).attr("stroke", color("--ink-2")).attr("stroke-width", 1).style("display", "none");
    const overlay = svg.append("rect").attr("class", "hit focusable").attr("x", LEFT).attr("y", TOP).attr("width", W - LEFT - RIGHT).attr("height", axisY - TOP)
      .attr("tabindex", 0).attr("role", "slider").attr("aria-label", `Residue position in ${d.gene}; use the left and right arrow keys`)
      .attr("aria-valuemin", v0).attr("aria-valuemax", v1);
    const show = (i, evt) => {
      i = Math.max(v0, Math.min(v1, i));
      state.hover = i;
      cross.style("display", null).attr("x1", x(i)).attr("x2", x(i));
      overlay.attr("aria-valuenow", i);
      showTip(evt || { clientX: overlay.node().getBoundingClientRect().left + x(i) - LEFT, clientY: overlay.node().getBoundingClientRect().top + 40 }, residueRows(i));
    };
    overlay.on("pointermove", (e) => show(Math.round(x.invert(d3.pointer(e)[0])), e))
      .on("pointerleave", () => { cross.style("display", "none"); hideTip(); })
      .on("keydown", (e) => {
        const k = e.key, cur0 = state.hover || Math.round((v0 + v1) / 2);
        if (k === "ArrowRight" || k === "ArrowLeft") { e.preventDefault(); show(cur0 + (k === "ArrowRight" ? 1 : -1) * (e.shiftKey ? 10 : 1)); }
      })
      .on("blur", () => { cross.style("display", "none"); hideTip(); });

    // overview strip with a brush to zoom
    {
      const xo = d3.scaleLinear().domain([0.5, d.L + 0.5]).range([LEFT, W - RIGHT]);
      const yo = d3.scaleLinear().domain([-1, 1]).range([ovY + ovH, ovY]);
      const all = d3.range(1, d.L + 1);
      svg.append("text").attr("class", "lab").attr("x", 0).attr("y", ovY + 12).text("Zoom");
      svg.append("text").attr("class", "lab").attr("x", 0).attr("y", ovY + 26).style("font-size", "10.5px").text("drag here");
      svg.append("rect").attr("x", LEFT).attr("y", ovY).attr("width", W - LEFT - RIGHT).attr("height", ovH).attr("fill", color("--surface-2"));
      svg.append("path").attr("d", d3.area().x((i) => xo(i)).y0(yo(0)).y1((i) => yo(Math.max(0, cur.tol[i - 1])))(all)).attr("fill", color("--tolerated")).attr("fill-opacity", 0.45);
      svg.append("path").attr("d", d3.area().x((i) => xo(i)).y0(yo(0)).y1((i) => yo(Math.min(0, cur.tol[i - 1])))(all)).attr("fill", color("--intolerant")).attr("fill-opacity", 0.45);
      const brush = d3.brushX().extent([[LEFT, ovY], [W - RIGHT, ovY + ovH]]).on("end", (e) => {
        if (!e.sourceEvent) return;
        if (!e.selection) { state.view = [1, d.L]; draw(); return; }
        const [a, b] = e.selection.map((p) => Math.round(xo.invert(p)));
        state.view = [Math.max(1, Math.min(a, b)), Math.min(d.L, Math.max(Math.max(a, b), Math.min(a, b) + 9))];
        draw();
      });
      const bg = svg.append("g").call(brush);
      if (v0 > 1 || v1 < d.L) bg.call(brush.move, [xo(v0 - 0.5), xo(v1 + 0.5)]);
      bg.selectAll(".selection").attr("fill", color("--accent")).attr("fill-opacity", 0.15).attr("stroke", color("--accent"));
      if (v0 > 1 || v1 < d.L) {
        const reset = svg.append("text").attr("class", "lab").attr("x", W - RIGHT).attr("y", ovY + ovH + 16).attr("text-anchor", "end")
          .style("cursor", "pointer").style("fill", color("--accent-ink")).attr("tabindex", 0).attr("role", "button")
          .text(`Showing ${v0}–${v1} · show the whole protein`);
        reset.on("click keydown", (e) => { if (e.type === "click" || e.key === "Enter") { state.view = [1, d.L]; draw(); } });
      }
    }
  }

  function residueRows(i) {
    const d = state.data, cur = state.current, off = state.official, aa = d.seq[i - 1];
    const rows = [{ head: `Residue ${i} · ${aa} (${AA3[aa] || aa}) · ${cur.disordered[i - 1] ? "disordered" : "folded"} region` }];
    rows.push({ key: color("--ink"), value: fmt.signed(cur.tol[i - 1]), label: cur === off ? "tolerance" : "tolerance, your recipe" });
    if (cur !== off) rows.push({ key: color("--muted"), value: fmt.signed(off.tol[i - 1]), label: "tolerance, official recipe" });
    rows.push({ key: color("--gnomad"), value: fmt.f(cur.gnomad[i - 1]), label: "gnomAD density" });
    rows.push({ key: color("--pathogenic"), value: fmt.f(cur.clinvar[i - 1]), label: "ClinVar density" });
    rows.push({ key: color("--esm"), value: fmt.f(d.esm[i - 1], 1), label: "ESM-1v mean LLR" });
    rows.push({ key: color("--disorder-line"), value: fmt.f(d.dis[i - 1]), label: "disorder (AIUPred)" });
    const rsa = d.rsa[i - 1];
    rows.push({ key: null, value: Number.isFinite(rsa) ? fmt.f(rsa) : "–", label: Number.isFinite(rsa) ? `RSA (${rsa < 0.2 ? "core" : rsa <= 0.4 ? "boundary" : "surface"})` : "RSA (no structure)" });
    rows.push({ key: null, value: Number.isFinite(d.plddt[i - 1]) ? String(d.plddt[i - 1]) : "–", label: "AlphaFold pLDDT" });
    const anns = [];
    for (const row of ANN_ROWS) for (const [key, label] of row.keys)
      if ((d.ann[key] || []).some(([s, e]) => i >= s && i <= e)) anns.push(label);
    if (anns.length) rows.push({ note: "Annotations: " + anns.join(", ") });
    const g = cur.gn.filter((v) => v.pos === i).map((v) => `${v.ref}${v.pos}${v.alt}`);
    const c = state.data.clinvar.filter((v) => v.pos === i).map((v) => `${v.ref}${v.pos}${v.alt} ${CLASS_NAME[v.flag]}${v.stars ? " ★" + v.stars : ""}`);
    if (g.length) rows.push({ note: "gnomAD: " + g.join(", ") });
    if (c.length) rows.push({ note: "ClinVar: " + c.join("; ") });
    return rows;
  }

  // ---------------------------------------------------------------- recipe controls
  function recipeControls() {
    const host = document.getElementById("recipe-controls");
    const slider = (label, key, min, max, stepv, fmtv) => {
      const out = el("output", { text: fmtv(state.opts[key]) });
      const input = el("input", { type: "range", min, max, step: stepv, value: state.opts[key], "aria-label": label });
      let raf = null;
      input.addEventListener("input", () => {
        state.opts[key] = +input.value; out.textContent = fmtv(+input.value);
        cancelAnimationFrame(raf);
        raf = requestAnimationFrame(() => { compute(); draw(); recipeStatus(); });
      });
      return { node: el("div", { class: "control" }, el("span", { text: label }), el("div", { class: "range-row" }, input, out)),
               set: (v) => { input.value = v; out.textContent = fmtv(v); } };
    };
    const span = (w) => `${w} (${2 * w + 1} res.)`;
    const sO = slider("Window, folded regions (w)", "wOrdered", 1, 61, 1, span);
    const sD = slider("Window, disordered regions (w)", "wDisordered", 1, 31, 1, span);
    const sT = slider("Disorder threshold", "threshold", 0.3, 0.8, 0.05, (v) => (+v).toFixed(2));
    const segHost = el("div", { class: "controls" });
    const mk = (label, key, options) => {
      const wrap = el("div", { class: "control" }, el("span", { text: label }));
      const seg = segmented(wrap, label, options, state.opts[key], (v) => { state.opts[key] = v; compute(); draw(); recipeStatus(); });
      segHost.appendChild(wrap);
      return seg;
    };
    const gA = mk("Variants at one residue", "aggregate", [{ value: "max", label: "most severe" }, { value: "sum", label: "sum" }, { value: "binary", label: "present / absent" }]);
    const gG = mk("gnomAD side", "gnomad", [{ value: "ac2", label: "seen ≥ 2 times" }, { value: "ac1", label: "include singletons" }]);
    const gC = mk("ClinVar side", "clinvar", [{ value: "all", label: "every entry" }, { value: "noBenign", label: "no benign" }, { value: "pathogenic", label: "pathogenic only" }]);
    const reset = el("button", { class: "btn small", type: "button", text: "Reset to the official recipe" });
    reset.addEventListener("click", () => {
      state.opts = { ...R.DEFAULTS };
      sO.set(state.opts.wOrdered); sD.set(state.opts.wDisordered); sT.set(state.opts.threshold);
      gA.set("max"); gG.set("ac2"); gC.set("all");
      compute(); draw(); recipeStatus();
    });
    host.replaceChildren(sO.node, sD.node, sT.node, segHost, el("div", { class: "control" }, el("span", { text: " " }), reset));
  }

  function recipeStatus() {
    const host = document.getElementById("recipe-status"), d = state.data, cur = state.current;
    if (!d) return;
    if (cur === state.official) {
      let sum = 0, sa = 0;
      for (const v of cur.tol) { sum += v; sa += Math.abs(v); }
      const ok = d.check && Math.abs(sum - d.check[0]) < 1e-4 && Math.abs(sa - d.check[1]) < 1e-4;
      host.replaceChildren(el("span", { text: "Official recipe. " }),
        ok ? el("span", { class: "ok", text: "✓ Rebuilt in your browser, it matches the saved score for this gene." })
           : el("span", { text: "Rebuilt in your browser (could not confirm against the saved score)." }));
    } else {
      const rho = S.spearman(cur.tol, state.official.tol);
      let diff = 0;
      for (let i = 0; i < cur.tol.length; i++) diff += Math.abs(cur.tol[i] - state.official.tol[i]);
      host.textContent = `Your recipe (black line) vs the official one (grey): within-gene Spearman ${fmt.f(rho)}, `
        + `mean change ${fmt.f(diff / cur.tol.length, 3)} per residue.`;
    }
  }

  // ---------------------------------------------------------------- labelled variants of this gene
  async function variantPanel() {
    const host = document.getElementById("variant-panel"), d = state.data;
    host.replaceChildren();
    if (d.split === "test") {
      if (!d.test.length) { host.appendChild(el("p", { class: "note", text: "No labelled test variants in this gene." })); return; }
      const rows = d.test.map(([v, benign, p, esm, tol]) => ({ v, benign, p, esm, tol, pos: +v.slice(1, -1) }));
      host.appendChild(el("p", { class: "fig-title", text: `Untouched test gene: the official predictions for its ${rows.length} labelled variants` }));
      host.appendChild(el("p", { class: "fig-sub", text: "P(benign) from the model fitted on the training and validation genes, before these variants were seen. Hover a row to mark the residue." }));
      const table = el("table", { class: "data-table" },
        el("thead", {}, el("tr", {}, ["Variant", "ClinVar", "Model P(benign)", "ESM-1v LLR", "Tolerance (own residue hidden)"].map((h) => el("th", { scope: "col", text: h })))),
        el("tbody", {}, rows.map((r) => {
          const tr = el("tr", {}, el("th", { scope: "row", text: r.v }), el("td", { text: r.benign ? "benign" : "pathogenic" }),
            el("td", { class: "num", text: fmt.f(r.p, 3) }), el("td", { class: "num", text: fmt.f(r.esm, 2) }), el("td", { class: "num", text: fmt.signed(r.tol, 3) }));
          tr.addEventListener("mouseenter", () => { state.marker = { pos: r.pos, label: r.v }; draw(); });
          return tr;
        })));
      host.appendChild(el("div", { class: "table-scroll" }, table));
      return;
    }
    const info = state.info;
    if (!info.labelled_pathogenic && !info.labelled_benign) {
      host.appendChild(el("p", { class: "note", text: "This gene has no labelled (≥ 1 star) pathogenic or benign ClinVar variants." }));
      return;
    }
    let model;
    try { model = await TOL.data.model(d.gene); } catch (e) { return; }
    if (state.gene !== d.gene) return;
    host.appendChild(el("p", { class: "fig-title", text: `Labelled variants in ${d.gene}: how the model scores them` }));
    host.appendChild(el("p", { class: "fig-sub", text: "Pick one to mark it on the tracks and see its features add up (the final model, which was trained with these variants)." }));
    const list = el("div", { class: "variant-list" }), chart = el("div", { class: "chart" });
    host.append(list, chart);
    const buttons = model.variants.map((v, k) => {
      const b = el("button", { class: "vbtn", type: "button" }, el("span", { class: "dot", style: `background:${color(v.benign ? "--benign" : "--pathogenic")}` }), v.v);
      b.addEventListener("click", () => select(k));
      list.appendChild(b);
      return b;
    });
    const select = (k) => {
      buttons.forEach((b, j) => b.classList.toggle("on", j === k));
      const v = model.variants[k];
      state.marker = { pos: +v.v.slice(1, -1), label: v.v };
      draw();
      TOL.waterfall(chart, d.gene, v);
    };
    const pick = model.variants.findIndex((v) => v.v === "A300P");
    select(pick >= 0 ? pick : 0);
  }

  // ---------------------------------------------------------------- set-up
  TOL.sections = TOL.sections || {};
  TOL.sections.explorer = async function () {
    const idx = await TOL.data.index();
    const list = document.getElementById("gene-list");
    list.replaceChildren(...idx.list.map((g) => el("option", { value: g.gene })));
    const input = document.getElementById("gene-input");
    const go = () => { const name = input.value.trim().toUpperCase(); if (idx.byName.has(name)) loadGene(name); };
    input.addEventListener("change", go);
    input.addEventListener("keydown", (e) => { if (e.key === "Enter") go(); });
    document.getElementById("gene-examples").replaceChildren(...EXAMPLES.filter(([g]) => idx.byName.has(g)).map(([g, why]) => {
      const b = el("button", { class: "chip", type: "button", title: why || null, text: g });
      b.addEventListener("click", () => loadGene(g));
      return b;
    }));
    legend();
    recipeControls();
    TOL.ui.onRedraw(() => { legend(); draw(); });
    const h = hashParams();
    await loadGene(h.gene && idx.byName.has(h.gene.toUpperCase()) ? h.gene.toUpperCase() : "GLRA1");
    if (h.gene) document.getElementById("explorer").scrollIntoView({ behavior: "instant", block: "start" });
  };
  TOL.explorer = { loadGene, state };
})();
