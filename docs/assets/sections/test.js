/* Experiment 4: the pre-registered test on untouched genes. */
(function () {
  "use strict";
  const TOL = window.TOL;
  const { el, fmt, color, showTip, hideTip, widthOf, segmented } = TOL.ui;
  const C = TOL.charts, S = TOL.stats;

  const CURVES = [
    { key: "model", label: "Logistic regression, 42 features (the model)", c: "--model", on: true },
    { key: "tolerance", label: "Tolerance alone", c: "--tol", on: true },
    { key: "esm", label: "ESM-1v alone", c: "--esm", on: true },
    { key: "esm_tolerance", label: "Logistic regression, ESM + tolerance", c: "--esmtol", on: false },
    { key: "trees", label: "Gradient-boosted trees, 42 features", c: "--trees", on: false, dash: "5 3" },
  ];
  const BIG3 = ["RYR1", "ABL1", "VHL"];

  TOL.sections = TOL.sections || {};
  TOL.sections.test = async function () {
    const test = await TOL.data.test();
    const rows = test.list, thr = test.esm_unsure_threshold;
    const st = { subset: "all", genes: { RYR1: true, ABL1: true, VHL: true, other: true }, curves: Object.fromEntries(CURVES.map((c) => [c.key, c.on])), cutoff: 0.5 };
    const nOther = new Set(rows.filter((r) => !BIG3.includes(r.gene)).map((r) => r.gene)).size;

    // ---------------------------------------------------------------- controls (one row; they scope everything below)
    const ctl = document.getElementById("test-controls");
    const subsetC = el("div", { class: "control" }, el("span", { text: "Variants" }));
    segmented(subsetC, "Variants", [{ value: "all", label: "All" }, { value: "unsure", label: `Where ESM-1v is unsure (LLR > ${fmt.f(thr, 2)})` }], st.subset, (v) => { st.subset = v; update(); });
    const geneC = el("div", { class: "control" }, el("span", { text: "Genes included" }));
    const chips = el("div", { class: "chips" });
    const count = (g) => rows.filter((r) => (g === "other" ? !BIG3.includes(r.gene) : r.gene === g)).length;
    [...BIG3, "other"].forEach((g) => {
      const b = el("button", { class: "chip on", type: "button", "aria-pressed": "true",
        text: g === "other" ? `other ${nOther} genes (${count(g)})` : `${g} (${count(g)})` });
      b.addEventListener("click", () => {
        st.genes[g] = !st.genes[g]; b.classList.toggle("on", st.genes[g]); b.setAttribute("aria-pressed", String(st.genes[g]));
        update();
      });
      chips.appendChild(b);
    });
    geneC.appendChild(chips);
    ctl.replaceChildren(subsetC, geneC);

    const selected = () => rows.filter((r) => (BIG3.includes(r.gene) ? st.genes[r.gene] : st.genes.other) && (st.subset === "all" || r.esm > thr));

    // ---------------------------------------------------------------- ROC
    const rocHost = document.getElementById("roc-chart"), legendHost = document.getElementById("roc-legend");
    function drawRoc() {
      const sel = selected(), ben = sel.filter((r) => r.benign), path = sel.filter((r) => !r.benign);
      legendHost.replaceChildren(...CURVES.map((cv) => {
        const auc = S.auroc(ben.map((r) => r[cv.key]), path.map((r) => r[cv.key]));
        const b = el("label", { class: "check" }, el("input", { type: "checkbox", checked: st.curves[cv.key] || null }),
          el("span", { class: "swatch", style: `background:${color(cv.c)};height:3px;${cv.dash ? "background:repeating-linear-gradient(90deg," + color(cv.c) + " 0 5px,transparent 5px 8px)" : ""}` }),
          `${cv.label}: `, el("strong", { text: Number.isFinite(auc) ? fmt.f(auc, 3) : "–" }));
        b.querySelector("input").addEventListener("change", (e) => { st.curves[cv.key] = e.target.checked; drawRoc(); });
        return b;
      }));
      legendHost.style.flexDirection = "column"; legendHost.style.alignItems = "flex-start";
      const W = widthOf(rocHost), size = Math.min(W, 420), m = { l: 46, r: 12, t: 8, b: 46 };
      const svg = C.svg(rocHost, size, size, "ROC curves on the untouched test genes");
      const x = d3.scaleLinear().domain([0, 1]).range([m.l, size - m.r]), y = d3.scaleLinear().domain([0, 1]).range([size - m.b, m.t]);
      C.axis(svg.append("g").attr("transform", `translate(0,${size - m.b})`), d3.axisBottom(x).ticks(5).tickSize(-(size - m.b - m.t)), { domain: false, grid: true });
      C.axis(svg.append("g").attr("transform", `translate(${m.l},0)`), d3.axisLeft(y).ticks(5).tickSize(-(size - m.l - m.r)), { domain: false, grid: true });
      C.axis(svg.append("g").attr("transform", `translate(0,${size - m.b})`), d3.axisBottom(x).ticks(5).tickSize(4));
      C.axis(svg.append("g").attr("transform", `translate(${m.l},0)`), d3.axisLeft(y).ticks(5).tickSize(4));
      svg.append("line").attr("x1", x(0)).attr("y1", y(0)).attr("x2", x(1)).attr("y2", y(1)).attr("stroke", color("--axis"));
      svg.append("text").attr("class", "lab").attr("x", (m.l + size - m.r) / 2).attr("y", size - 10).attr("text-anchor", "middle").text("pathogenic variants called benign (share)");
      svg.append("text").attr("class", "lab").attr("transform", `translate(12,${(m.t + size - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("benign variants called benign (share)");
      if (!ben.length || !path.length) {
        svg.append("text").attr("class", "lab").attr("x", size / 2).attr("y", size / 2).attr("text-anchor", "middle").text("Needs both benign and pathogenic variants");
        return;
      }
      const lineGen = d3.line().curve(d3.curveStepAfter).x((p) => x(p[0])).y((p) => y(p[1]));
      CURVES.slice().reverse().forEach((cv) => {
        if (!st.curves[cv.key]) return;
        const pts = S.rocCurve(ben.map((r) => r[cv.key]), path.map((r) => r[cv.key]));
        svg.append("path").attr("d", lineGen(pts)).attr("fill", "none").attr("stroke", color(cv.c)).attr("stroke-width", cv.key === "model" ? 2.5 : 2)
          .attr("stroke-dasharray", cv.dash || null).attr("stroke-linejoin", "round");
      });
    }

    // ---------------------------------------------------------------- the S-curve and a cut-off
    const sHost = document.getElementById("scurve-chart"), cutCtl = document.getElementById("cutoff-controls");
    const out = el("output", { text: st.cutoff.toFixed(2) });
    const slider = el("input", { type: "range", min: 0.05, max: 0.95, step: 0.01, value: st.cutoff, "aria-label": "Cut-off for P(benign)" });
    slider.addEventListener("input", () => { st.cutoff = +slider.value; out.textContent = st.cutoff.toFixed(2); drawS(); });
    cutCtl.replaceChildren(el("div", { class: "control" }, el("span", { text: "Call a variant benign when P(benign) ≥" }), el("div", { class: "range-row" }, slider, out)),
      el("div", { class: "legend" },
        el("span", { class: "item" }, el("span", { class: "swatch box", style: `background:${color("--pathogenic")};border-radius:50%` }), "ClinVar pathogenic"),
        el("span", { class: "item" }, el("span", { class: "swatch box", style: `background:${color("--benign")};border-radius:50%` }), "ClinVar benign")));
    const logit = (p) => Math.log(Math.max(1e-6, p) / Math.max(1e-6, 1 - p));
    const jitter = (s) => { let h = 0; for (const ch of s) h = (h * 31 + ch.charCodeAt(0)) >>> 0; return ((h % 1000) / 1000 - 0.5) * 0.06; };
    function drawS() {
      const sel = selected(), W = widthOf(sHost), Hh = Math.min(W, 380), m = { l: 46, r: 12, t: 10, b: 42 };
      const svg = C.svg(sHost, W, Hh, "Test variants on the model's S-curve");
      const x = d3.scaleLinear().domain([-9, 9]).range([m.l, W - m.r]), y = d3.scaleLinear().domain([-0.05, 1.05]).range([Hh - m.b, m.t]);
      C.axis(svg.append("g").attr("transform", `translate(0,${Hh - m.b})`), d3.axisBottom(x).ticks(7).tickSize(4));
      C.axis(svg.append("g").attr("transform", `translate(${m.l},0)`), d3.axisLeft(y).tickValues([0, 0.25, 0.5, 0.75, 1]).tickFormat(d3.format(".2~f")).tickSize(4));
      svg.append("text").attr("class", "lab").attr("x", (m.l + W - m.r) / 2).attr("y", Hh - 8).attr("text-anchor", "middle").text("z = weighted sum of the 42 features");
      svg.append("text").attr("class", "lab").attr("transform", `translate(12,${(m.t + Hh - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("P(benign)");
      svg.append("rect").attr("x", m.l).attr("width", W - m.l - m.r).attr("y", y(1.05)).attr("height", y(st.cutoff) - y(1.05)).attr("fill", color("--tint-blue"));
      svg.append("path").attr("d", d3.line().x((z) => x(z)).y((z) => y(1 / (1 + Math.exp(-z))))(d3.range(-9, 9.01, 0.1))).attr("fill", "none").attr("stroke", color("--ink-2")).attr("stroke-width", 2);
      svg.append("line").attr("x1", m.l).attr("x2", W - m.r).attr("y1", y(st.cutoff)).attr("y2", y(st.cutoff)).attr("stroke", color("--ink")).attr("stroke-width", 1.5);
      svg.append("text").attr("class", "lab").attr("x", m.l + 6).attr("y", y(st.cutoff) - 6).text("called benign ↑");
      svg.append("text").attr("class", "lab").attr("x", m.l + 6).attr("y", y(st.cutoff) + 14).text("called pathogenic ↓");
      const pts = sel.map((r) => ({ r, cx: x(Math.max(-9, Math.min(9, logit(r.model)))), cy: y(r.model + jitter(r.gene + r.variant)) }));
      const order = pts.slice().sort((a, b) => b.r.benign - a.r.benign);      // pathogenic drawn last, on top
      svg.append("g").selectAll("circle").data(order).join("circle").attr("cx", (p) => p.cx).attr("cy", (p) => p.cy).attr("r", 3)
        .attr("fill", (p) => color(p.r.benign ? "--benign" : "--pathogenic")).attr("fill-opacity", 0.6);
      if (pts.length) {
        const dl = d3.Delaunay.from(pts, (p) => p.cx, (p) => p.cy);
        const ring = svg.append("circle").attr("r", 6).attr("fill", "none").attr("stroke", color("--ink")).attr("stroke-width", 1.5).style("display", "none");
        svg.append("rect").attr("class", "hit").attr("x", m.l).attr("y", m.t).attr("width", W - m.l - m.r).attr("height", Hh - m.t - m.b)
          .on("pointermove", (e) => {
            const [px, py] = d3.pointer(e), k = dl.find(px, py), p = pts[k];
            if (!p || Math.hypot(p.cx - px, p.cy - py) > 24) { ring.style("display", "none"); hideTip(); return; }
            ring.style("display", null).attr("cx", p.cx).attr("cy", p.cy);
            showTip(e, [{ head: `${p.r.gene} ${p.r.variant} · ClinVar ${p.r.benign ? "benign" : "pathogenic"}` },
              { key: color("--model"), value: fmt.f(p.r.model, 3), label: "model P(benign)" },
              { key: color("--esm"), value: fmt.f(p.r.esm, 2), label: "ESM-1v LLR" },
              { key: color("--tol"), value: fmt.signed(p.r.tolerance, 3), label: "tolerance (own residue hidden)" }]);
          }).on("pointerleave", () => { ring.style("display", "none"); hideTip(); });
      }
      const ben = sel.filter((r) => r.benign), path = sel.filter((r) => !r.benign);
      const cleared = ben.filter((r) => r.model >= st.cutoff).length, caught = path.filter((r) => r.model < st.cutoff).length;
      document.getElementById("cutoff-readout").replaceChildren(
        `At P(benign) ≥ ${st.cutoff.toFixed(2)}: `, el("strong", { text: `${fmt.pct(cleared / ben.length)} of benign variants cleared` }),
        ` (${cleared} of ${ben.length}) and `, el("strong", { text: `${fmt.pct(caught / path.length)} of pathogenic variants caught` }),
        ` (${caught} of ${path.length}); accuracy ${fmt.pct((cleared + caught) / sel.length)}.`);
    }

    // ---------------------------------------------------------------- tables
    function questions() {
      const q = test.results.questions;
      const NAMES = { "Q1: primary beats ESM-1v on all variants": "Q1 · The model beats ESM-1v on all variants",
        "Q2: primary beats ESM-1v within genes": "Q2 · The model beats ESM-1v within genes (8 genes)",
        "Q3: primary beats ESM-1v where ESM is unsure": "Q3 · The model beats ESM-1v where ESM-1v is unsure",
        "Q4: human variant data adds to ESM": "Q4 · Human-variant data add to ESM-1v (ESM + tolerance vs ESM)" };
      const table = el("table", { class: "data-table" },
        el("thead", {}, el("tr", {}, ["Question", "Difference in AUROC", "95% range", "Verdict"].map((h) => el("th", { scope: "col", text: h })))),
        el("tbody", {}, Object.entries(q).map(([k, v]) => el("tr", {}, el("th", { scope: "row", text: NAMES[k] || k }),
          el("td", { class: "num", text: fmt.signed(v.difference, 3) }),
          el("td", { class: "num", text: `${fmt.signed(v.lo, 3)} to ${fmt.signed(v.hi, 3)}` }),
          el("td", {}, v.passes ? el("strong", { text: "✓ passes" }) : "range includes 0: not shown")))));
      document.getElementById("questions-table").replaceChildren(table);
    }
    function geneTable() {
      const byGene = d3.groups(rows, (r) => r.gene).map(([g, rs]) => {
        const b = rs.filter((r) => r.benign), p = rs.filter((r) => !r.benign), ok = b.length >= 3 && p.length >= 3;
        return { g, n: rs.length, b: b.length, p: p.length,
                 model: ok ? S.auroc(b.map((r) => r.model), p.map((r) => r.model)) : NaN,
                 esm: ok ? S.auroc(b.map((r) => r.esm), p.map((r) => r.esm)) : NaN };
      }).sort((a, b) => b.n - a.n);
      const qualifying = byGene.filter((r) => Number.isFinite(r.model));
      const table = el("table", { class: "data-table" },
        el("thead", {}, el("tr", {}, ["Gene", "Variants", "Share of test", "Pathogenic", "Benign", "Within-gene AUROC: model", "ESM-1v"].map((h) => el("th", { scope: "col", text: h })))),
        el("tbody", {}, byGene.map((r) => el("tr", { class: BIG3.includes(r.g) ? "hi" : null },
          el("th", { scope: "row" }, el("a", { href: "#gene=" + r.g, text: r.g, onclick: (e) => { e.preventDefault(); TOL.explorer.loadGene(r.g); document.getElementById("explorer").scrollIntoView(); } })),
          el("td", { class: "num", text: String(r.n) }), el("td", { class: "num", text: fmt.pct(r.n / rows.length, 1) }),
          el("td", { class: "num", text: String(r.p) }), el("td", { class: "num", text: String(r.b) }),
          el("td", { class: "num", text: Number.isFinite(r.model) ? fmt.f(r.model, 3) : "" }),
          el("td", { class: "num", text: Number.isFinite(r.esm) ? fmt.f(r.esm, 3) : "" })))),
        el("tfoot", {}, el("tr", {}, el("th", { scope: "row", text: `Mean over the ${qualifying.length} qualifying genes` }), el("td"), el("td"), el("td"), el("td"),
          el("td", { class: "num", text: fmt.f(d3.mean(qualifying, (r) => r.model), 3) }), el("td", { class: "num", text: fmt.f(d3.mean(qualifying, (r) => r.esm), 3) }))));
      document.getElementById("test-genes-table").replaceChildren(table);
    }

    function update() { drawRoc(); drawS(); }
    C.live(update);
    questions();
    geneTable();
  };
})();
