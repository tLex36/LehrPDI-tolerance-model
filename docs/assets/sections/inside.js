/* Inside the model: the weights, and one variant's prediction feature by feature. */
(function () {
  "use strict";
  const TOL = window.TOL;
  const { el, fmt, color, showTip, hideTip, tableToggle, widthOf, segmented } = TOL.ui;
  const C = TOL.charts;
  let MS = null;                                   // model summary (weights, intercept)
  const GROUP_NAME = { "ESM scores": "evolution (ESM)", "human variants nearby": "human variants nearby",
    "gene-wide variant rates": "gene-wide rates", "substitution type": "substitution type", "structure (§7)": "structure",
    "annotations (§8)": "annotations", "disorder + conservation": "disorder + conservation", "missing-value flags": "missing-value flag" };

  // ---------------------------------------------------------------- weights
  function weights() {
    const host = document.getElementById("weights-chart"), st = { sort: "size", all: false, flags: false };
    const ctl = document.getElementById("weights-controls");
    const sortC = el("div", { class: "control" }, el("span", { text: "Order" }));
    segmented(sortC, "Order", [{ value: "size", label: "largest first" }, { value: "group", label: "by kind of data" }], st.sort, (v) => { st.sort = v; draw(); });
    const allC = el("label", { class: "check" }, el("input", { type: "checkbox" }), "show all features, not just the 15 largest");
    allC.querySelector("input").addEventListener("change", (e) => { st.all = e.target.checked; draw(); });
    const flagC = el("label", { class: "check" }, el("input", { type: "checkbox" }), "include the 4 \"value missing\" flags");
    flagC.querySelector("input").addEventListener("change", (e) => { st.flags = e.target.checked; draw(); });
    ctl.replaceChildren(sortC, el("div", { class: "control" }, el("span", { text: " " }), allC), el("div", { class: "control" }, el("span", { text: " " }), flagC));
    const draw = () => {
      let rows = MS.features.filter((f) => st.flags || !f.missing_flag);
      rows = rows.slice().sort((a, b) => Math.abs(b.weight) - Math.abs(a.weight));
      if (!st.all) rows = rows.slice(0, 15);
      if (st.sort === "group") rows.sort((a, b) => MS.groups.indexOf(a.group) - MS.groups.indexOf(b.group) || Math.abs(b.weight) - Math.abs(a.weight));
      const W = widthOf(host), narrow = W < 560, left = narrow ? 4 : 300, right = narrow ? 40 : 54, rowH = narrow ? 34 : 24;
      let extra = 0;
      if (st.sort === "group") extra = new Set(rows.map((r) => r.group)).size * 20;
      const Hh = rows.length * rowH + extra + (narrow ? 70 : 54);
      const svg = C.svg(host, W, Hh, "Weights of the final model");
      const mx = d3.max(MS.features, (f) => Math.abs(f.weight));
      const x = d3.scaleLinear().domain([-mx, mx]).nice().range([left, W - right]);
      const base = Hh - (narrow ? 50 : 34), nt = narrow ? 4 : 6;
      C.axis(svg.append("g").attr("transform", `translate(0,${base})`), d3.axisBottom(x).ticks(nt).tickSize(-(base - 10)), { domain: false, grid: true });
      C.axis(svg.append("g").attr("transform", `translate(0,${base})`), d3.axisBottom(x).ticks(nt).tickFormat(d3.format("+.1f")).tickSize(4));
      svg.append("line").attr("x1", x(0)).attr("x2", x(0)).attr("y1", 4).attr("y2", base).attr("stroke", color("--axis"));
      svg.append("text").attr("class", "lab").attr("x", W - right).attr("y", narrow ? Hh - 18 : Hh - 4).attr("text-anchor", "end").text("pushes towards benign →");
      svg.append("text").attr("class", "lab").attr("x", left).attr("y", Hh - 4).text("← pushes towards pathogenic");
      let y = 6, lastGroup = null;
      rows.forEach((f) => {
        if (st.sort === "group" && f.group !== lastGroup) {
          y += 6;
          svg.append("text").attr("class", "lab-strong").attr("x", 0).attr("y", y + 10).text(GROUP_NAME[f.group] || f.group);
          y += 14; lastGroup = f.group;
        }
        const c = color(f.weight >= 0 ? "--tolerated" : "--intolerant");
        const by = narrow ? y + 14 : y;
        svg.append("text").attr("class", "lab").attr("x", st.sort === "group" ? 12 : 0).attr("y", narrow ? y + 10 : y + 13).style("font-size", narrow ? "11px" : null).text(f.label);
        svg.append("path").attr("d", C.hBar(x(0), x(f.weight), by + 3, narrow ? 12 : 14)).attr("fill", c);
        svg.append("text").attr("class", "lab").attr("x", f.weight >= 0 ? x(f.weight) + 5 : x(f.weight) - 5).attr("y", by + 13)
          .attr("text-anchor", f.weight >= 0 ? "start" : "end").style("font-variant-numeric", "tabular-nums").text(fmt.signed(f.weight, 2));
        svg.append("rect").attr("class", "hit").attr("x", 0).attr("y", y).attr("width", W).attr("height", rowH)
          .on("pointermove", (e) => showTip(e, [{ head: f.label }, { key: c, value: fmt.signed(f.weight, 3), label: "weight per standard deviation" },
            { note: `Kind of data: ${GROUP_NAME[f.group] || f.group}. ${f.weight >= 0 ? "Higher values push towards benign." : "Higher values push towards pathogenic."}` }]))
          .on("pointerleave", hideTip);
        y += rowH;
      });
    };
    C.live(draw);
    tableToggle(document.getElementById("weights-card"), () => [["Feature", "Kind of data", "Weight"],
      ...MS.features.slice().sort((a, b) => Math.abs(b.weight) - Math.abs(a.weight)).map((f) => [f.label, GROUP_NAME[f.group] || f.group, fmt.signed(f.weight, 3)])]);
  }

  // ---------------------------------------------------------------- one variant, feature by feature
  /** Draw the waterfall for one labelled variant (from data/model/<GENE>.js) into host. */
  async function waterfall(host, gene, v) {
    if (!MS) MS = await TOL.data.modelSummary();
    host._wf = { gene, v };
    const draw = () => {
      const { gene, v } = host._wf;
      const steps = [{ label: "starting point (intercept)", value: MS.intercept, raw: null }]
        .concat(v.top.map(([j, c, raw]) => ({ label: MS.features[j].label, value: c, raw, key: MS.features[j].key, flag: MS.features[j].missing_flag })))
        .concat([{ label: `the other ${MS.features.length - v.top.length} features together`, value: v.rest, raw: null }]);
      let acc = 0;
      steps.forEach((s) => { s.from = acc; acc += s.value; s.to = acc; });
      const z = acc, p = 1 / (1 + Math.exp(-z));
      const W = widthOf(host), narrow = W < 560, left = narrow ? 4 : 320, right = narrow ? 44 : 70, rowH = narrow ? 36 : 24, Hh = (steps.length + 1) * rowH + (narrow ? 84 : 70);
      const svg = C.svg(host, W, Hh, `How the model scores ${gene} ${v.v}`);
      const ext = d3.extent(steps.flatMap((s) => [s.from, s.to]).concat([0]));
      const x = d3.scaleLinear().domain([Math.min(-1, ext[0]), Math.max(1, ext[1])]).nice().range([left, W - right]);
      svg.append("text").attr("class", "title-text").attr("x", 0).attr("y", 14)
        .text(`${gene} p.${v.v} · ClinVar ${v.benign ? "benign" : "pathogenic"}${v.stars ? " (" + "★".repeat(Math.min(4, v.stars)) + ")" : ""}`);
      const verdict = `the model leans ${p >= 0.5 ? "benign" : "pathogenic"}${(p >= 0.5) === !!v.benign ? ", as ClinVar says" : ", unlike ClinVar"}.`;
      svg.append("text").attr("class", "lab").attr("x", 0).attr("y", 32).text(`z = ${fmt.signed(z, 2)} → P(benign) = ${fmt.f(p, 3)}${narrow ? ":" : ": " + verdict}`);
      if (narrow) svg.append("text").attr("class", "lab").attr("x", 0).attr("y", 48).text(verdict);
      const top = narrow ? 60 : 44;
      const nt = narrow ? 4 : 6;
      C.axis(svg.append("g").attr("transform", `translate(0,${Hh - 22})`), d3.axisBottom(x).ticks(nt).tickSize(-(Hh - 22 - top)), { domain: false, grid: true });
      C.axis(svg.append("g").attr("transform", `translate(0,${Hh - 22})`), d3.axisBottom(x).ticks(nt).tickFormat(d3.format("+.1~f")).tickSize(4));
      svg.append("line").attr("x1", x(0)).attr("x2", x(0)).attr("y1", top).attr("y2", Hh - 22).attr("stroke", color("--axis"));
      steps.forEach((s, i) => {
        const y = top + i * rowH + (narrow ? 12 : 0), c = color(s.value >= 0 ? "--tolerated" : "--intolerant");
        const rawTxt = s.raw === null || s.raw === undefined ? "" : s.flag ? `? ${s.raw ? "yes" : "no"}`
          : ` = ${Math.abs(s.raw) >= 100 ? fmt.f(s.raw, 0) : +(+s.raw).toPrecision(3)}`;
        const label = s.label + rawTxt;
        svg.append("text").attr("class", "lab").attr("x", 0).attr("y", narrow ? y : y + 15).style("font-size", narrow ? "11px" : null).text(label);
        const a = Math.min(s.from, s.to), b = Math.max(s.from, s.to);
        svg.append("rect").attr("x", x(a)).attr("width", Math.max(1.5, x(b) - x(a))).attr("y", y + 4).attr("height", 14).attr("rx", 2).attr("fill", c);
        if (i < steps.length - 1) svg.append("line").attr("x1", x(s.to)).attr("x2", x(s.to)).attr("y1", y + 18).attr("y2", y + rowH + 4).attr("stroke", color("--ink-2"));
        svg.append("text").attr("class", "lab").attr("x", x(b) + 5).attr("y", y + 15).style("font-variant-numeric", "tabular-nums").text(fmt.signed(s.value, 2));
        svg.append("rect").attr("class", "hit").attr("x", 0).attr("y", y).attr("width", W).attr("height", rowH)
          .on("pointermove", (e) => showTip(e, [{ head: s.label }, { key: c, value: fmt.signed(s.value, 3), label: "added to z" },
            s.raw !== null && s.raw !== undefined ? { key: null, value: s.flag ? (s.raw ? "yes" : "no") : String(s.raw), label: "value of the feature for this variant" }
              : { note: i === 0 ? "The model's log-odds for a variant whose every feature is at the training average." : "" },
            { key: null, value: fmt.signed(s.to, 2), label: "running total" }]))
          .on("pointerleave", hideTip);
      });
      const yz = top + steps.length * rowH + (narrow ? 12 : 0);
      svg.append("text").attr("class", "lab-strong").attr("x", 0).attr("y", yz + 15).text("total: z");
      svg.append("circle").attr("cx", x(z)).attr("cy", yz + 11).attr("r", 6).attr("fill", color(p >= 0.5 ? "--tolerated" : "--intolerant")).attr("stroke", color("--surface")).attr("stroke-width", 2);
      svg.append("text").attr("class", "lab-strong").attr("x", x(z) + 10).attr("y", yz + 15).text(`${fmt.signed(z, 2)} → P(benign) ${fmt.f(p, 3)}`);
    };
    if (!host._wfLive) { host._wfLive = true; TOL.ui.onRedraw(() => host._wf && draw()); }
    draw();
  }

  async function waterfallSection() {
    const idx = await TOL.data.index();
    const genes = idx.list.filter((g) => (g.split === "train" || g.split === "val") && (g.labelled_pathogenic + g.labelled_benign) > 0);
    const ctl = document.getElementById("waterfall-controls"), list = document.getElementById("waterfall-variants"), chart = document.getElementById("waterfall-chart");
    const input = el("input", { type: "search", list: "waterfall-genes", value: "GLRA1", "aria-label": "Gene", autocomplete: "off", spellcheck: "false" });
    const dl = el("datalist", { id: "waterfall-genes" }, genes.map((g) => el("option", { value: g.gene })));
    ctl.replaceChildren(el("div", { class: "control" }, el("label", { text: "Gene (training or validation)" }), input, dl),
      el("div", { class: "control" }, el("span", { text: " " }), el("span", { class: "note", text: `${genes.length} genes have labelled variants` })));
    const show = async (name, pick) => {
      if (!genes.some((g) => g.gene === name)) return;
      const m = await TOL.data.model(name);
      const buttons = m.variants.map((v, k) => {
        const b = el("button", { class: "vbtn", type: "button" }, el("span", { class: "dot", style: `background:${color(v.benign ? "--benign" : "--pathogenic")}` }), v.v);
        b.addEventListener("click", () => select(k));
        return b;
      });
      list.replaceChildren(...buttons);
      const select = (k) => { buttons.forEach((b, j) => b.classList.toggle("on", j === k)); waterfall(chart, name, m.variants[k]); };
      const k = m.variants.findIndex((v) => v.v === pick);
      select(k >= 0 ? k : 0);
    };
    input.addEventListener("change", () => show(input.value.trim().toUpperCase()));
    await show("GLRA1", "A300P");
  }

  TOL.waterfall = waterfall;
  TOL.sections = TOL.sections || {};
  TOL.sections.inside = async function () {
    MS = await TOL.data.modelSummary();
    weights();
    await waterfallSection();
  };
  TOL.insideReady = () => MS !== null;
})();
