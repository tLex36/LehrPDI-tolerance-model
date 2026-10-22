/* Context: conservation across species (species tree + heatmap), as on the poster. */
(function () {
  "use strict";
  const TOL = window.TOL;
  const { el, fmt, color, showTip, hideTip, widthOf } = TOL.ui;
  const C = TOL.charts;
  const SHOW_CLADES = ["Mammalia", "Amniota", "Vertebrata", "Ecdysozoa"];

  /** Minimal Newick parser: "(A:1,(B:2,C:2):3);" -> {name, length, children} */
  function parseNewick(s) {
    let i = 0;
    function node() {
      const n = { children: [], name: "", length: 0 };
      if (s[i] === "(") {
        i++;
        n.children.push(node());
        while (s[i] === ",") { i++; n.children.push(node()); }
        i++;                                               // ")"
      }
      let name = "";
      while (i < s.length && !",):;".includes(s[i])) name += s[i++];
      n.name = name.trim().replace(/^'|'$/g, "");
      if (s[i] === ":") {
        i++;
        let num = "";
        while (i < s.length && !",);".includes(s[i])) num += s[i++];
        n.length = +num;
      }
      return n;
    }
    return node();
  }

  TOL.sections = TOL.sections || {};
  TOL.sections.conservation = async function () {
    const s = await TOL.data.summary(), cons = s.conservation;
    if (!cons) return;
    const host = document.getElementById("cons-chart");
    const tree = parseNewick(cons.newick);
    // ages: distance from the tips (the tree is ultrametric, branch lengths in million years)
    (function depth(n, d) { n.depth = d; n.children.forEach((c) => depth(c, d + c.length)); })(tree, 0);
    const tips = [];
    (function collect(n) { if (!n.children.length) tips.push(n); n.children.forEach(collect); })(tree);
    const total = d3.max(tips, (t) => t.depth);
    const tipIndex = new Map(cons.species.map((sp, i) => [sp, i]));
    (function place(n) {
      n.age = total - n.depth;
      if (!n.children.length) { n.row = tipIndex.has(n.name) ? tipIndex.get(n.name) : tips.indexOf(n); n.tips = [n.name]; return; }
      n.children.forEach(place);
      n.row = (d3.min(n.children, (c) => c.row) + d3.max(n.children, (c) => c.row)) / 2;
      n.tips = n.children.flatMap((c) => c.tips);
    })(tree);
    const labelOf = new Map();
    for (const [name, list] of Object.entries(cons.clades || {})) if (SHOW_CLADES.includes(name)) labelOf.set([...list].sort().join("|"), name);
    const ramp = d3.scaleSequential().domain([0, 100]).interpolator(d3.interpolateRgbBasis(["#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"]));
    const nProt = 3;

    document.getElementById("cons-legend").replaceChildren(
      el("span", { class: "item" }, "0%", el("span", { style: `display:inline-block;width:140px;height:10px;border-radius:3px;background:linear-gradient(90deg,${d3.range(0, 101, 20).map((v) => ramp(v)).join(",")})` }), "100% conserved vs human"),
      el("span", { class: "item" }, el("span", { class: "swatch box", style: `background:repeating-linear-gradient(45deg,${color("--axis")} 0 2px,transparent 2px 5px)` }), "too few orthologues, or no genome alignment"));

    const draw = () => {
      const W = widthOf(host), narrow = W < 680, rowH = 22, top = 44;
      const treeW = narrow ? Math.max(170, W * 0.36) : Math.min(380, W * 0.36), heatW = W - treeW - 8, cw = heatW / cons.columns.length;
      const Hh = top + cons.species.length * rowH + 44;
      const svg = C.svg(host, W, Hh, "Conservation across species");
      const nameW = narrow ? 92 : 120;
      const x = d3.scaleLinear().domain([total, 0]).range([10, treeW - nameW]);
      const y = (row) => top + row * rowH + rowH / 2;
      const g = svg.append("g");
      (function edges(n) {
        n.children.forEach((c) => {
          g.append("path").attr("d", `M${x(n.age)},${y(c.row)}H${x(c.age)}`).attr("stroke", color("--ink-2")).attr("fill", "none").attr("stroke-width", 1.3);
          edges(c);
        });
        if (n.children.length) g.append("line").attr("x1", x(n.age)).attr("x2", x(n.age)).attr("y1", y(d3.min(n.children, (c) => c.row))).attr("y2", y(d3.max(n.children, (c) => c.row))).attr("stroke", color("--ink-2")).attr("stroke-width", 1.3);
        const lab = n.children.length ? labelOf.get([...n.tips].sort().join("|")) : null;
        if (lab) {
          const fits = x(n.age) - 3 > lab.length * 6;             // otherwise write it to the right of the node
          g.append("text").attr("class", "lab").attr("x", fits ? x(n.age) - 3 : x(n.age) + 3).attr("y", y(n.row) - 4)
            .attr("text-anchor", fits ? "end" : "start").style("font-size", "10.5px").style("font-style", "italic").text(lab);
        }
      })(tree);
      tips.forEach((t) => g.append("text").attr("class", "lab").attr("x", x(0) + 4).attr("y", y(t.row) + 4).style("font-size", narrow ? "11px" : "12px")
        .style("font-style", /^[A-Z]\. /.test(t.name) || t.name === "Drosophila" || t.name === "Xenopus tropicalis" ? "italic" : null).text(t.name));
      C.axis(svg.append("g").attr("transform", `translate(0,${top + cons.species.length * rowH + 6})`), d3.axisBottom(x).ticks(narrow ? 3 : 4).tickSize(4));
      svg.append("text").attr("class", "lab").attr("x", (x(total) + x(0)) / 2).attr("y", Hh - 4).attr("text-anchor", "middle").style("font-size", "11px").text("million years ago");
      // heatmap
      const hx = treeW + 8;
      const defs = svg.append("defs");
      defs.append("pattern").attr("id", "hatch").attr("width", 6).attr("height", 6).attr("patternUnits", "userSpaceOnUse").attr("patternTransform", "rotate(45)")
        .append("line").attr("x1", 0).attr("y1", 0).attr("x2", 0).attr("y2", 6).attr("stroke", color("--axis")).attr("stroke-width", 2);
      const SHORT = { "Postsynaptic proteins": "Post- syn.", "Olfactory receptors": "Olf. rec.", "Ribosomal proteins": "Ribo- somal",
                      "Regulatory elements": "Reg. DNA", "Coding regions": "Cod. DNA" };
      svg.append("text").attr("class", "title-text").attr("x", hx + (cw * nProt) / 2).attr("y", 12).attr("text-anchor", "middle").text(narrow ? "Protein" : "Protein: % similarity");
      svg.append("text").attr("class", "title-text").attr("x", hx + cw * nProt + (cw * (cons.columns.length - nProt)) / 2).attr("y", 12).attr("text-anchor", "middle").text(narrow ? "DNA" : "DNA: % identical");
      cons.columns.forEach((col, j) => {
        const words = (narrow ? SHORT[col] || col : col).split(" ");
        words.forEach((w, k) => svg.append("text").attr("class", "lab").attr("x", hx + j * cw + cw / 2).attr("y", 26 + k * 11).attr("text-anchor", "middle").style("font-size", narrow ? "9px" : "10.5px").text(w));
      });
      cons.species.forEach((sp, i) => {
        cons.columns.forEach((col, j) => {
          const v = cons.values[i][j], cov = cons.coverage[i][j], prot = j < nProt;
          const hatched = v === null || (prot && cov !== null && cov < 10);
          const cx = hx + j * cw + 1, cy = top + i * rowH + 1, w = cw - 2, h = rowH - 2;
          const cell = svg.append("g").attr("tabindex", 0).attr("class", "focusable");
          cell.append("rect").attr("x", cx).attr("y", cy).attr("width", w).attr("height", h).attr("rx", 2)
            .attr("fill", hatched ? "url(#hatch)" : ramp(v)).attr("stroke", hatched ? color("--grid") : "none");
          const dark = !hatched && v > 55;
          if (hatched) {
            const full = v === null && !prot ? "no alignment" : cov === 0 ? "none" : `${fmt.f(cov, 0)}% have one`;
            cell.append("text").attr("class", "lab").attr("x", cx + w / 2).attr("y", cy + 14).attr("text-anchor", "middle").style("font-size", "9.5px")
              .text(cw < 70 ? (v === null && !prot ? "–" : cov === 0 ? "none" : `${fmt.f(cov, 0)}%`) : full);
          } else {
            const t = cell.append("text").attr("x", cx + w / 2 - (prot && sp !== "Human" && cw > 64 ? 8 : 0)).attr("y", cy + 14).attr("text-anchor", "middle")
              .style("font-size", "11.5px").style("font-weight", 650).style("fill", dark ? "#ffffff" : "#0b0b0b").text(fmt.f(v, 0));
            if (prot && sp !== "Human" && cw > 64) cell.append("text").attr("x", cx + w / 2 + 6).attr("y", cy + 14).style("font-size", "9px").style("fill", dark ? "#ffffff" : "#0b0b0b").text(`${fmt.f(cov, 0)}%`);
            void t;
          }
          cell.on("pointermove focus", (e) => showTip(e.type === "focus" ? e.target : e, [{ head: `${sp} · ${col}` },
            { key: hatched ? null : ramp(v ?? 0), value: v === null ? "no alignment" : `${fmt.f(v, 1)}%`, label: prot ? "median similarity" : "identical bases" },
            prot ? { key: null, value: cov === null ? "–" : `${fmt.f(cov, 0)}%`, label: "of genes have an orthologue" } : { note: "Mean over 80 regions" }]))
            .on("pointerleave blur", hideTip);
        });
      });
    };
    C.live(draw);
    TOL.ui.tableToggle(host.parentNode, () => [["Species", ...cons.columns.map((c, j) => (j < nProt ? `${c} (similarity; % with orthologue)` : `${c} (% identical)`))],
      ...cons.species.map((sp, i) => [sp, ...cons.columns.map((c, j) => {
        const v = cons.values[i][j], cov = cons.coverage[i][j];
        if (v === null) return "no alignment";
        return j < nProt ? `${fmt.f(v, 0)}% (${fmt.f(cov, 0)}%)` : `${fmt.f(v, 1)}%`;
      })])]);
  };
})();
