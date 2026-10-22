/* Experiments 1-3 and Design A: held-out variants, structural classes, model choice, feature groups. */
(function () {
  "use strict";
  const TOL = window.TOL;
  const { el, fmt, color, showTip, hideTip, tableToggle, widthOf, segmented } = TOL.ui;
  const C = TOL.charts, S = TOL.stats;

  // ================================================================ experiment 1: held-out variants
  const GROUPS = [
    { key: "gnomad", label: "gnomAD (healthy population)", c: "--gnomad" },
    { key: "benign", label: "ClinVar benign", c: "--benign" },
    { key: "vus", label: "ClinVar VUS", c: "--vus" },
    { key: "pathogenic", label: "ClinVar pathogenic", c: "--pathogenic" },
    { key: "other", label: "ClinVar other (conflicting, …)", c: "--other" },
    { key: "clinvar", label: "ClinVar, any class", c: "--ink-2" },
  ];
  const ORDER = ["gnomad", "benign", "vus", "pathogenic", "other"];

  function heldout(ms) {
    const H = ms.heldout, st = { score: "tolerance", subset: "all", a: "gnomad", b: "pathogenic" };
    const controls = document.getElementById("heldout-controls"), host = document.getElementById("heldout-chart");
    const wrap = (label, node) => el("div", { class: "control" }, el("span", { text: label }), node);
    const scoreHost = el("div"), subsetHost = el("div");
    segmented(scoreHost, "Score", [{ value: "tolerance", label: "Tolerance (from 90%)" }, { value: "esm", label: "ESM-1v LLR" }], st.score, (v) => { st.score = v; draw(); });
    segmented(subsetHost, "Variants", [{ value: "all", label: "All held-out" }, { value: "poster", label: "As on the poster" }], st.subset, (v) => { st.subset = v; draw(); });
    const select = (key) => {
      const s = el("select", { "aria-label": key === "a" ? "First group" : "Second group" },
        GROUPS.map((g) => el("option", { value: g.key, text: g.label, selected: g.key === st[key] || null })));
      s.addEventListener("change", () => { st[key] = s.value; draw(); });
      return s;
    };
    controls.replaceChildren(wrap("Score", scoreHost), wrap("Variants", subsetHost), wrap("Compare", select("a")), wrap("with", select("b")));

    const block = () => H[st.subset].scores[st.score];
    const countsOf = (blk, key) => key === "clinvar"
      ? ["benign", "vus", "pathogenic", "other"].map((k) => blk.groups[k]).filter(Boolean)
          .reduce((acc, g) => acc.map((v, i) => v + g.counts[i]), new Array(blk.groups.gnomad.counts.length).fill(0))
      : (blk.groups[key] ? blk.groups[key].counts : null);
    const meanOf = (blk, key) => {
      if (key !== "clinvar") return blk.groups[key] ? blk.groups[key].mean : NaN;
      const parts = ["benign", "vus", "pathogenic", "other"].map((k) => blk.groups[k]).filter(Boolean);
      return d3.sum(parts, (g) => g.mean * g.n) / d3.sum(parts, (g) => g.n);
    };
    const nOf = (blk, key) => d3.sum(countsOf(blk, key) || []);

    function draw() {
      const blk = block(), W = widthOf(host), narrow = W < 600, left = narrow ? 8 : 196, right = narrow ? 24 : 64, rowH = narrow ? 80 : 54;
      const rows = ORDER.concat(st.a === "clinvar" || st.b === "clinvar" ? ["clinvar"] : []).filter((k) => countsOf(blk, k));
      const Hh = rows.length * rowH + 52;
      const svg = C.svg(host, W, Hh, "Held-out variants: distribution of the score in each group");
      const nb = blk.groups.gnomad.counts.length, edges = d3.range(nb + 1).map((i) => blk.lo + i * blk.step);
      // show the range where any group has data
      let first = nb, last = 0;
      rows.forEach((k) => countsOf(blk, k).forEach((v, i) => { if (v) { first = Math.min(first, i); last = Math.max(last, i); } }));
      const x = d3.scaleLinear().domain([edges[Math.max(0, first - 1)], edges[Math.min(nb, last + 2)]]).range([left, W - right]);
      rows.forEach((k, ri) => {
        const g = GROUPS.find((q) => q.key === k), counts = countsOf(blk, k), n = d3.sum(counts);
        const dens = counts.map((v) => v / n / blk.step), y0 = 18 + ri * rowH + rowH - 8;
        const y = d3.scaleLinear().domain([0, d3.max(dens)]).range([y0, y0 - rowH + (narrow ? 40 : 14)]);
        const on = k === st.a || k === st.b;
        const area = d3.area().curve(d3.curveStepAfter).x((_, i) => x(edges[i])).y0(y0).y1((v) => y(v));
        svg.append("path").attr("d", area(dens.concat([0]))).attr("fill", color(g.c)).attr("fill-opacity", on ? 0.75 : 0.18);
        svg.append("line").attr("x1", left).attr("x2", W - right).attr("y1", y0).attr("y2", y0).attr("stroke", color("--axis"));
        const m = meanOf(blk, k);
        svg.append("line").attr("x1", x(m)).attr("x2", x(m)).attr("y1", y0 + 2).attr("y2", y0 - rowH + 12).attr("stroke", color("--ink")).attr("stroke-width", on ? 2 : 1).attr("stroke-opacity", on ? 1 : 0.5);
        svg.append("text").attr("class", on ? "lab-strong" : "lab").attr("x", narrow ? left : 0).attr("y", narrow ? y0 - rowH + 20 : y0 - 20).text(g.label);
        svg.append("text").attr("class", "lab").attr("x", narrow ? left : 0).attr("y", narrow ? y0 - rowH + 34 : y0 - 5)
          .style("font-size", "11px").text(`n = ${fmt.int(n)} · mean ${fmt.signed(m)}`);
        svg.append("rect").attr("x", left).attr("y", y0 - rowH + 10).attr("width", W - left - right).attr("height", rowH - 10).attr("class", "hit")
          .on("pointermove", (e) => {
            const i = Math.max(0, Math.min(nb - 1, Math.floor((x.invert(d3.pointer(e)[0]) - blk.lo) / blk.step)));
            showTip(e, [{ head: g.label }, { key: color(g.c), value: fmt.pct(counts[i] / n, 1), label: `of variants between ${fmt.f(edges[i], 2)} and ${fmt.f(edges[i + 1], 2)}` },
              { key: null, value: fmt.signed(m), label: "mean" }]);
          }).on("pointerleave", hideTip);
      });
      C.axis(svg.append("g").attr("transform", `translate(0,${Hh - 34})`), d3.axisBottom(x).ticks(8).tickSize(4));
      svg.append("text").attr("class", "lab").attr("x", W - right).attr("y", Hh - 4).attr("text-anchor", "end").style("font-size", "11px")
        .text(st.score === "tolerance" ? "tolerance (higher = more tolerant) →" : "ESM-1v LLR (higher = less surprising swap) →");
      // readout
      const ca = countsOf(blk, st.a), cb = countsOf(blk, st.b), A = GROUPS.find((q) => q.key === st.a), B = GROUPS.find((q) => q.key === st.b);
      const out = document.getElementById("heldout-readout");
      if (st.a === st.b) { out.textContent = "Pick two different groups to compare."; return; }
      const auc = S.binnedAuroc(ca, cb);
      const ia = ORDER.indexOf(st.a), ib = ORDER.indexOf(st.b);
      let exact = null;
      if (ia >= 0 && ib >= 0) exact = ia < ib ? blk.exact_auroc[`${st.a}|${st.b}`] : 1 - blk.exact_auroc[`${st.b}|${st.a}`];
      const what = st.score === "tolerance" ? "sits in the more tolerant region" : "gets the less surprising ESM-1v score";
      const better = auc >= 0.5 ? A : B, share = auc >= 0.5 ? auc : 1 - auc;
      out.replaceChildren(el("strong", { text: `AUROC ${fmt.f(auc, 3)}` }),
        exact !== null ? ` (exact ${fmt.f(exact, 3)})` : " (from 0.01-wide bins)",
        `: in ${fmt.pct(share)} of random pairs the ${better.label} variant ${what}. `,
        `Means ${fmt.signed(meanOf(blk, st.a))} (${A.label}, n = ${fmt.int(nOf(blk, st.a))}) vs ${fmt.signed(meanOf(blk, st.b))} (${B.label}, n = ${fmt.int(nOf(blk, st.b))}).`);
    }
    C.live(draw);
    tableToggle(document.getElementById("heldout-card"), () => {
      const out = [["Group", "Variants", "Tolerance: mean", "Tolerance: median", "ESM-1v LLR: mean", "ESM-1v LLR: median"]];
      for (const k of ORDER) {
        const t = H[st.subset].scores.tolerance.groups[k], e = H[st.subset].scores.esm.groups[k];
        if (t) out.push([GROUPS.find((q) => q.key === k).label, fmt.int(t.n), fmt.signed(t.mean, 3), fmt.signed(t.median, 3), fmt.f(e.mean, 2), fmt.f(e.median, 2)]);
      }
      return out;
    });
  }

  // ================================================================ experiment 2: structural classes
  function structure(s) {
    const host = document.getElementById("struct-chart");
    const rows = s.struct_classes.map(([key, label, ben, path]) => ({ key, label, ben, path, n: ben + path, p: path / (ben + path) }))
      .sort((a, b) => b.p - a.p);
    const draw = () => {
      const W = widthOf(host), narrow = W < 600, left = narrow ? 0 : Math.min(270, W * 0.42), right = 70, rowH = narrow ? 48 : 34, Hh = rows.length * rowH + 44;
      const svg = C.svg(host, W, Hh, "Share of held-out ClinVar variants that are pathogenic, by structural class");
      const x = d3.scaleLinear().domain([0, 1]).range([left, W - right]);
      const lg = svg.append("g").attr("transform", `translate(${left},10)`);
      [["pathogenic", "--pathogenic"], ["benign", "--benign"]].forEach(([lab, c], i) => {
        lg.append("rect").attr("x", i * 110).attr("y", -6).attr("width", 10).attr("height", 10).attr("rx", 2).attr("fill", color(c));
        lg.append("text").attr("class", "lab").attr("x", i * 110 + 14).attr("y", 3).text(lab);
      });
      rows.forEach((r, i) => {
        const y = 26 + i * rowH + (narrow ? 16 : 0);
        svg.append("text").attr("class", "lab").attr("x", 0).attr("y", narrow ? y - 4 : y + 13).text(r.label);
        const g = svg.append("g").attr("tabindex", 0).attr("class", "focusable")
          .on("pointermove focus", (e) => showTip(e.type === "focus" ? e.target : e, [{ head: r.label },
            { key: color("--pathogenic"), value: `${fmt.pct(r.p)} (${r.path})`, label: "pathogenic" },
            { key: color("--benign"), value: `${fmt.pct(1 - r.p)} (${r.ben})`, label: "benign" }]))
          .on("pointerleave blur", hideTip);
        g.append("path").attr("d", `M${x(0)},${y}H${x(r.p) - 1}V${y + 18}H${x(0)}Z`).attr("fill", color("--pathogenic"));
        g.append("path").attr("d", C.hBar(x(r.p) + 1, x(1), y, 18)).attr("fill", color("--benign")).attr("fill-opacity", 0.55);
        svg.append("text").attr("class", "lab-strong").attr("x", W - right + 8).attr("y", y + 13).text(fmt.pct(r.p));
        svg.append("text").attr("class", "lab").attr("x", W - right + 44).attr("y", y + 13).style("font-size", "11px").text(`n ${r.n}`);
      });
      C.axis(svg.append("g").attr("transform", `translate(0,${Hh - 16})`), d3.axisBottom(x).ticks(5).tickFormat(d3.format(".0%")).tickSize(4));
    };
    C.live(draw);
    tableToggle(host.parentNode, () => [["Class", "Pathogenic", "Benign", "% pathogenic"], ...rows.map((r) => [r.label, r.path, r.ben, fmt.pct(r.p, 1)])]);
  }

  // ================================================================ experiment 3: models (cross-validation and test)
  const TEST_KEY = ["ESM-1v LLR", "tolerance (own residue hidden)", "logistic: ESM LLR + tolerance", "logistic: all features",
                    "gradient-boosted trees: all features"];
  function models(s, test) {
    const host = document.getElementById("cv-chart"), res = test.results.scores;
    const rows = s.cv.map(([name, inputs, all, within, unsure, testAll], i) => {
      const r = TEST_KEY[i] ? res[TEST_KEY[i]] : null;
      return { name, inputs, chosen: i === 3, cv: { all, within, unsure },
               test: r ? { all: r.AUROC, within: r["within-gene AUROC"], unsure: r["AUROC where ESM is unsure"] } : null };
    });
    const st = { metric: "all" };
    const LABEL = { all: "All variants", within: "Within genes", unsure: "Where ESM-1v is unsure" };
    const ctl = document.getElementById("cv-controls");
    const m = el("div", { class: "control" }, el("span", { text: "AUROC measured on" }));
    segmented(m, "Metric", Object.entries(LABEL).map(([value, label]) => ({ value, label })), st.metric, (v) => { st.metric = v; draw(); });
    ctl.replaceChildren(m, el("div", { class: "legend" },
      el("span", { class: "item" }, el("span", { class: "swatch box", style: `background:${color("--model")};border-radius:50%` }), "cross-validation (training + validation genes)"),
      el("span", { class: "item" }, el("span", { class: "swatch box", style: `border:2px solid ${color("--ink")};border-radius:50%;background:transparent` }), "untouched test genes")));
    const draw = () => {
      const W = widthOf(host), narrow = W < 620, left = narrow ? 150 : 300, right = 104, rowH = narrow ? 40 : 34, Hh = rows.length * rowH + 40;
      const svg = C.svg(host, W, Hh, "AUROC of each model");
      const lo = st.metric === "unsure" ? 0.5 : 0.6;
      const x = d3.scaleLinear().domain([lo, 1]).range([left, W - right]);
      const nt = narrow ? 2 : 5;
      C.axis(svg.append("g").attr("transform", `translate(0,${Hh - 20})`), d3.axisBottom(x).ticks(nt).tickSize(-(Hh - 26)), { domain: false, grid: true });
      C.axis(svg.append("g").attr("transform", `translate(0,${Hh - 20})`), d3.axisBottom(x).ticks(nt).tickSize(4));
      rows.forEach((r, i) => {
        const y = 10 + i * rowH + rowH / 2;
        if (r.chosen) svg.append("rect").attr("x", 0).attr("y", y - rowH / 2 + 2).attr("width", W).attr("height", rowH - 4).attr("rx", 6).attr("fill", color("--tint-blue"));
        svg.append("text").attr("class", r.chosen ? "lab-strong" : "lab").attr("x", 6).attr("y", narrow ? y - 3 : y + 4).text(r.name);
        svg.append("text").attr("class", "lab").attr("x", narrow ? 6 : 176).attr("y", narrow ? y + 12 : y + 4).style("font-size", "11px").text(r.inputs);
        const v = r.cv[st.metric], t = r.test ? r.test[st.metric] : null;
        if (t !== null && Number.isFinite(t)) {
          svg.append("line").attr("x1", x(v)).attr("x2", x(Math.max(lo, t))).attr("y1", y).attr("y2", y).attr("stroke", color("--axis")).attr("stroke-width", 2);
          svg.append("circle").attr("cx", x(Math.max(lo, t))).attr("cy", y).attr("r", 6).attr("fill", color("--surface")).attr("stroke", color("--ink")).attr("stroke-width", 2);
        }
        svg.append("circle").attr("cx", x(v)).attr("cy", y).attr("r", 6).attr("fill", color("--model")).attr("stroke", color("--surface")).attr("stroke-width", 2);
        svg.append("text").attr("class", "lab").attr("x", W - right + 10).attr("y", y + 4).style("font-variant-numeric", "tabular-nums")
          .text(`${fmt.f(v, 3)}${t !== null && Number.isFinite(t) ? " · " + fmt.f(t, 3) : " · –"}`);
        svg.append("rect").attr("class", "hit").attr("x", 0).attr("y", y - rowH / 2).attr("width", W).attr("height", rowH)
          .on("pointermove", (e) => showTip(e, [{ head: `${r.name}: ${r.inputs}` },
            { key: color("--model"), value: fmt.f(v, 3), label: "cross-validation" },
            { key: color("--ink"), value: t !== null && Number.isFinite(t) ? fmt.f(t, 3) : "not tested", label: "untouched test genes" },
            r.chosen ? { note: "Chosen in advance as the one model to test." } : (t === null ? { note: "Left out of the test on purpose: it did not beat logistic regression." } : { note: "" })]))
          .on("pointerleave", hideTip);
      });
      svg.append("text").attr("class", "lab").attr("x", W - right + 10).attr("y", 8).style("font-size", "11px").text("CV · test");
    };
    C.live(draw);
    tableToggle(document.getElementById("cv-card"), () => [["Model", "Inputs", "CV: all", "CV: within genes", "CV: ESM unsure", "Test: all", "Test: within genes", "Test: ESM unsure"],
      ...rows.map((r) => [r.name, r.inputs, fmt.f(r.cv.all, 3), fmt.f(r.cv.within, 3), fmt.f(r.cv.unsure, 3),
        r.test ? fmt.f(r.test.all, 3) : "not tested", r.test ? fmt.f(r.test.within, 3) : "", r.test ? fmt.f(r.test.unsure, 3) : ""])]);
  }

  // ================================================================ feature groups
  const GROUP_NAME = { "ESM scores": "ESM scores (evolution)", "human variants nearby": "Human variants nearby",
    "gene-wide variant rates": "Gene-wide variant rates", "substitution type": "Substitution type",
    "structure (§7)": "Structure", "annotations (§8)": "Annotations", "disorder + conservation": "Disorder + conservation" };
  function groups(s) {
    const host = document.getElementById("groups-chart"), st = { view: "added" };
    const entries = Object.entries(s.feature_groups).map(([k, v]) => ({ key: k, name: GROUP_NAME[k] || k, ...v }));
    const ctl = document.getElementById("groups-controls");
    const m = el("div", { class: "control" }, el("span", { text: "Show" }));
    segmented(m, "View", [{ value: "added", label: "Added to the ESM scores" }, { value: "alone", label: "Each group alone" }], st.view, (v) => { st.view = v; draw(); });
    ctl.replaceChildren(m, el("div", { class: "legend" },
      el("span", { class: "item" }, el("span", { class: "swatch box", style: `background:${color("--model")}` }), "all variants"),
      el("span", { class: "item" }, el("span", { class: "swatch box", style: `background:${color("--boundary")};opacity:.45` }), "within genes")));
    const draw = () => {
      const rows = (st.view === "added" ? entries.filter((e) => Number.isFinite(e.gain)).sort((a, b) => b.gain - a.gain)
                                         : entries.slice().sort((a, b) => b.alone - a.alone));
      const W = widthOf(host), narrow = W < 600, left = narrow ? 4 : Math.min(200, W * 0.36), right = 60, rowH = narrow ? 50 : 36, Hh = rows.length * rowH + 36;
      const svg = C.svg(host, W, Hh, "AUROC by feature group");
      const x = st.view === "added" ? d3.scaleLinear().domain([-0.015, 0.08]).range([left, W - right])
                                    : d3.scaleLinear().domain([0.5, 1]).range([left, W - right]);
      C.axis(svg.append("g").attr("transform", `translate(0,${Hh - 20})`), d3.axisBottom(x).ticks(5).tickSize(-(Hh - 24)), { domain: false, grid: true });
      C.axis(svg.append("g").attr("transform", `translate(0,${Hh - 20})`), d3.axisBottom(x).ticks(5).tickFormat(st.view === "added" ? d3.format("+.2f") : d3.format(".2f")).tickSize(4));
      const base = st.view === "added" ? x(0) : x(0.5);
      if (st.view === "added") svg.append("line").attr("x1", x(0)).attr("x2", x(0)).attr("y1", 4).attr("y2", Hh - 20).attr("stroke", color("--axis"));
      rows.forEach((r, i) => {
        const y = 8 + i * rowH + (narrow ? 14 : 0), main = st.view === "added" ? r.gain : r.alone, within = st.view === "added" ? r.gain_within : r.alone_within;
        svg.append("text").attr("class", "lab").attr("x", 0).attr("y", narrow ? y - 2 : y + 15).text(r.name);
        svg.append("path").attr("d", C.hBar(base, x(main), y + 2, 14)).attr("fill", color("--model"));
        svg.append("path").attr("d", C.hBar(base, x(within), y + 18, 7, 3)).attr("fill", color("--boundary")).attr("fill-opacity", 0.45);
        if (st.view === "added") {
          svg.append("line").attr("x1", x(r.gain_lo)).attr("x2", x(r.gain_hi)).attr("y1", y + 9).attr("y2", y + 9).attr("stroke", color("--ink")).attr("stroke-width", 1.5);
          [r.gain_lo, r.gain_hi].forEach((v) => svg.append("line").attr("x1", x(v)).attr("x2", x(v)).attr("y1", y + 5).attr("y2", y + 13).attr("stroke", color("--ink")).attr("stroke-width", 1.5));
        }
        svg.append("text").attr("class", "lab-strong").attr("x", Math.max(x(main), st.view === "added" ? x(r.gain_hi) : 0) + 8).attr("y", y + 14)
          .text(st.view === "added" ? fmt.signed(main, 3) : fmt.f(main, 3));
        svg.append("rect").attr("class", "hit").attr("x", 0).attr("y", y).attr("width", W).attr("height", rowH)
          .on("pointermove", (e) => showTip(e, st.view === "added"
            ? [{ head: r.name + ", added to the ESM scores" }, { key: color("--model"), value: fmt.signed(r.gain, 3), label: `AUROC gained (95% range ${fmt.signed(r.gain_lo, 3)} to ${fmt.signed(r.gain_hi, 3)})` },
               { key: color("--boundary"), value: fmt.signed(r.gain_within, 3), label: "gained within genes" }]
            : [{ head: r.name + ", alone" }, { key: color("--model"), value: fmt.f(r.alone, 3), label: "AUROC, all variants" },
               { key: color("--boundary"), value: fmt.f(r.alone_within, 3), label: "AUROC within genes" }]))
          .on("pointerleave", hideTip);
      });
    };
    C.live(draw);
    tableToggle(document.getElementById("groups-card"), () => [["Group", "Alone", "Alone, within genes", "Gain added to ESM", "95% range", "Gain within genes"],
      ...entries.map((r) => [r.name, fmt.f(r.alone, 3), fmt.f(r.alone_within, 3), fmt.signed(r.gain, 3),
        Number.isFinite(r.gain_lo) ? `${fmt.signed(r.gain_lo, 3)} to ${fmt.signed(r.gain_hi, 3)}` : "–", fmt.signed(r.gain_within, 3)])]);
  }

  // ================================================================ Design A
  function designA(s) {
    const a = s.design_a, host = document.getElementById("designa-chart"), host2 = document.getElementById("recipe-chart");
    const rows = [
      { label: "ESM-1v scores, no training", lo: a.esm_zero_shot, hi: a.esm_zero_shot },
      { label: "Networks trained from scratch on sequence", lo: a.from_scratch[0], hi: a.from_scratch[1] },
      { label: "Small networks on frozen ESM-1v (best: BiLSTM head)", lo: a.esm_heads[0], hi: a.esm_heads[1] },
      { label: "Trees with all 97 per-residue inputs, incl. structure", lo: a.trees_all_inputs, hi: a.trees_all_inputs },
    ];
    const draw = () => {
      const W = widthOf(host), narrow = W < 640, left = narrow ? 4 : Math.min(330, W * 0.46), right = 60, rowH = narrow ? 48 : 34, Hh = rows.length * rowH + 50;
      const svg = C.svg(host, W, Hh, "Design A results");
      const x = d3.scaleLinear().domain([0, 0.45]).range([left, W - right]);
      C.axis(svg.append("g").attr("transform", `translate(0,${Hh - 20})`), d3.axisBottom(x).ticks(5).tickSize(4));
      svg.append("line").attr("x1", x(a.ceiling)).attr("x2", x(a.ceiling)).attr("y1", 8).attr("y2", Hh - 20).attr("stroke", color("--pathogenic")).attr("stroke-width", 2);
      svg.append("text").attr("class", "lab").attr("x", x(a.ceiling) - 6).attr("y", 14).attr("text-anchor", "end").text(`noise ceiling ${fmt.f(a.ceiling)}`);
      rows.forEach((r, i) => {
        const y = 22 + i * rowH + (narrow ? 14 : 0);
        svg.append("text").attr("class", "lab").attr("x", 0).attr("y", narrow ? y - 2 : y + 13).text(r.label);
        svg.append("path").attr("d", C.hBar(x(0), x(r.hi), y + 2, 14)).attr("fill", color("--model")).attr("fill-opacity", r.lo === r.hi ? 1 : 0.35);
        if (r.lo !== r.hi) svg.append("path").attr("d", C.hBar(x(0), x(r.lo), y + 2, 14)).attr("fill", color("--model"));
        const t = (v) => String(+v.toFixed(3));
        svg.append("text").attr("class", "lab-strong").attr("x", x(r.hi) + 8).attr("y", y + 13).text(r.lo === r.hi ? t(r.hi) : `${t(r.lo)}–${t(r.hi)}`);
      });
    };
    C.live(draw);
    const p = a.pathogenic_only;
    const pairs = [
      { label: "How reproducible the target is (noise ceiling)", before: p.ceiling_before, after: p.ceiling_after },
      { label: "What a simple model learns (linear probe, Spearman)", before: p.probe_before, after: p.probe_after },
      { label: "The track itself on held-out ClinVar variants (AUROC)", before: p.auroc_before, after: p.auroc_after },
    ];
    const draw2 = () => {
      const W = widthOf(host2), narrow = W < 640, left = narrow ? 10 : Math.min(360, W * 0.5), right = 90, rowH = narrow ? 50 : 36, Hh = pairs.length * rowH + 56;
      const svg = C.svg(host2, W, Hh, "Effect of a pathogenic-only ClinVar side");
      const x = d3.scaleLinear().domain([0, 1]).range([left, W - right]);
      C.axis(svg.append("g").attr("transform", `translate(0,${Hh - 20})`), d3.axisBottom(x).ticks(5).tickSize(4));
      const lg = svg.append("g").attr("transform", `translate(${left},10)`);
      [["every ClinVar entry (official)", "--esm"], ["pathogenic only", "--model"]].forEach(([lab, c], i) => {
        const dx = narrow ? 0 : i * 190, dy = narrow ? i * 16 - 4 : 0;
        lg.append("circle").attr("cx", dx + 5).attr("cy", dy).attr("r", 5).attr("fill", color(c));
        lg.append("text").attr("class", "lab").attr("x", dx + 14).attr("y", dy + 4).text(lab);
      });
      pairs.forEach((r, i) => {
        const y = 36 + i * rowH + (narrow ? 16 : 0);
        svg.append("text").attr("class", "lab").attr("x", 0).attr("y", narrow ? y - 12 : y + 4).text(r.label);
        svg.append("line").attr("x1", x(r.before)).attr("x2", x(r.after)).attr("y1", y).attr("y2", y).attr("stroke", color("--axis")).attr("stroke-width", 3);
        svg.append("circle").attr("cx", x(r.before)).attr("cy", y).attr("r", 6).attr("fill", color("--esm")).attr("stroke", color("--surface")).attr("stroke-width", 2);
        svg.append("circle").attr("cx", x(r.after)).attr("cy", y).attr("r", 6).attr("fill", color("--model")).attr("stroke", color("--surface")).attr("stroke-width", 2);
        svg.append("text").attr("class", "lab-strong").attr("x", x(r.after) + 10).attr("y", y + 4).text(`${fmt.f(r.before)} → ${fmt.f(r.after)}`);
      });
    };
    C.live(draw2);
  }

  TOL.sections = TOL.sections || {};
  TOL.sections.experiments = async function () {
    const [s, ms, test] = await Promise.all([TOL.data.summary(), TOL.data.modelSummary(), TOL.data.test()]);
    heldout(ms);
    structure(s);
    models(s, test);
    groups(s);
    designA(s);
  };
})();
