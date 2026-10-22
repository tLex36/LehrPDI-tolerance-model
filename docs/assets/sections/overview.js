/* "The data": the pipeline, the sources and the gene split. */
(function () {
  "use strict";
  const TOL = window.TOL;
  const { el, fmt, color, showTip, hideTip, tableToggle, widthOf } = TOL.ui;
  const C = TOL.charts;

  const USES = { tolerance: "tolerance score", features: "model features", labels: "labels", context: "context" };
  const SOURCES = [
    { name: "gnomAD v4.1", what: "Variants in about 800,000 people, with those known to have severe childhood disease excluded: a healthy reference population.",
      took: (s) => `${fmt.int(s.variants.gnomad_total)} missense variants in these genes; ${fmt.int(s.variants.gnomad_used)} used (seen at least twice, plus a rule for X-linked genes).`,
      uses: ["tolerance", "features"] },
    { name: "ClinVar", what: "Variants reported by clinical laboratories, each classified (pathogenic, benign, uncertain = VUS) with a review status in stars.",
      took: (s) => `${fmt.int(s.variants.clinvar_total)} missense entries of any class for the tolerance score; ${fmt.int(s.variants.labelled_trainval + s.variants.labelled_test)} pathogenic or benign with at least one star as labels.`,
      uses: ["tolerance", "features", "labels"] },
    { name: "SynGO", what: "An expert-curated knowledge base of synapse genes.",
      took: () => "1,082 postsynaptic proteins, one MANE Select transcript each.", uses: ["context"] },
    { name: "ESM-1v", what: "A protein language model (650 million weights) trained to predict hidden amino acids in natural protein sequences.",
      took: () => "The probability of every amino acid at every residue, turned into 4 features (the variant's LLR and 3 per-residue summaries).",
      uses: ["features"] },
    { name: "AIUPred", what: "Predicts intrinsically disordered (unstructured) regions, and disordered regions likely to bind a partner, from sequence.",
      took: () => "A disorder and a binding score per residue. Disorder ≥ 0.5 switches the tolerance recipe to its disordered-region settings.",
      uses: ["tolerance", "features"] },
    { name: "PDB (via PDBe) and AlphaFold DB", what: "Experimental 3D structures, and predicted structures with a per-residue confidence (pLDDT).",
      took: () => "Solvent accessibility (RSA) from the best experimental chain, else from AlphaFold where pLDDT ≥ 70; pLDDT itself.",
      uses: ["features"] },
    { name: "PDBe-KB", what: "Interaction interfaces pooled across all PDB entries.",
      took: () => "Whether a residue touches another chain: any, a copy of the same protein, or a different protein.", uses: ["features"] },
    { name: "InterPro and UniProt", what: "Protein domains and families, and curated sites.",
      took: () => "Domains, repeats, transmembrane segments, signal peptides, ligand and active sites, PTMs, disulfides. UniProt's natural variants and mutagenesis were left out so clinical labels can't leak in.",
      uses: ["features"] },
    { name: "MobiDB", what: "Curated and predicted disorder, and linear interacting peptides (short binding stretches in disordered regions).",
      took: () => "Curated disorder, binding peptides and a conservation score per residue.", uses: ["features"] },
    { name: "ELM", what: "The Eukaryotic Linear Motif resource: short motifs recognised by partner domains.",
      took: () => "Curated and predicted motif flags.", uses: ["features"] },
    { name: "Ensembl 116", what: "Orthologues and pairwise genome alignments across species.",
      took: () => "Conservation across 17 species, for context on the poster (not a model input).", uses: ["context"] },
  ];

  function pipeline(s) {
    const steps = [
      ["1", "Variants", "gnomAD and ClinVar missense variants in 1,082 postsynaptic proteins", "1.15M", "#explorer"],
      ["2", "Tolerance score", "One value per residue, from where healthy and clinical variants fall", fmt.int(780750) + " residues", "#explorer"],
      ["3", "Features", "42 numbers per variant: evolution, human variation, structure, annotations", "42", "#inside"],
      ["4", "Labels and splits", "ClinVar pathogenic or benign; gene families kept together", "15,339 + 1,357", "#split-card"],
      ["5", "Models and test", "8 models compared; the chosen one tested once on untouched genes", "0.89 AUROC", "#test"],
    ];
    const host = document.getElementById("pipeline");
    host.replaceChildren(...steps.map(([n, t, d, big, href]) =>
      el("a", { class: "step", href }, el("div", { class: "n", text: "Step " + n }), el("div", { class: "t", text: t }),
        el("div", { class: "d", text: d }), el("div", { class: "big", text: big }))));
  }

  function sources(s) {
    const grid = document.getElementById("sources"), chips = document.getElementById("source-filter");
    const cards = SOURCES.map((src) => el("div", { class: "source", "data-uses": src.uses.join(" ") },
      el("h4", { text: src.name }), el("p", { text: src.what }), el("div", { class: "took", text: src.took(s) }),
      el("div", { class: "uses" }, src.uses.map((u) => el("span", { class: "use", text: USES[u] })))));
    grid.replaceChildren(...cards);
    let active = "all";
    const options = [["all", "everything"], ...Object.entries(USES)];
    const buttons = options.map(([key, label]) => {
      const b = el("button", { class: "chip" + (key === "all" ? " on" : ""), type: "button", "aria-pressed": String(key === "all"), text: label });
      b.addEventListener("click", () => {
        active = key;
        buttons.forEach((x, i) => { const on = options[i][0] === key; x.classList.toggle("on", on); x.setAttribute("aria-pressed", String(on)); });
        cards.forEach((c) => c.classList.toggle("dim", key !== "all" && !c.dataset.uses.split(" ").includes(key)));
      });
      return b;
    });
    chips.replaceChildren(...buttons);
  }

  async function splitChart(s) {
    const idx = await TOL.data.index();
    const NAMES = { train: "Training", val: "Validation", test: "Test (untouched)" };
    const rows = ["train", "val", "test"].map((k) => {
      const genes = idx.list.filter((g) => g.split === k);
      return { key: k, name: NAMES[k], genes: s.split[k].genes, residues: s.split[k].residues,
               path: genes.reduce((a, g) => a + g.labelled_pathogenic, 0), ben: genes.reduce((a, g) => a + g.labelled_benign, 0) };
    });
    const host = document.getElementById("split-chart");
    const panels = [
      { title: "Genes", parts: [["genes", "--gnomad", "genes"]] },
      { title: "Residues", parts: [["residues", "--gnomad", "residues"]] },
      { title: "Labelled variants", parts: [["path", "--pathogenic", "pathogenic"], ["ben", "--benign", "benign"]] },
    ];
    host.before(el("div", { class: "legend" },
      el("span", { class: "item" }, el("span", { class: "swatch box", style: `background:${color("--pathogenic")}` }), "labelled pathogenic"),
      el("span", { class: "item" }, el("span", { class: "swatch box", style: `background:${color("--benign")}` }), "labelled benign")));
    const draw = () => {
      const W = widthOf(host), narrow = W < 640, labelW = 118;
      const pw = narrow ? W - labelW - 10 : (W - labelW - 10) / 3 - 16, rowH = 30, ph = rows.length * rowH + 34;
      const H = narrow ? panels.length * (ph + 10) + 20 : ph + 24;
      const svg = C.svg(host, W, H, "Genes, residues and labelled variants per split");
      panels.forEach((p, pi) => {
        const ox = narrow ? labelW : labelW + pi * (pw + 16), oy = narrow ? pi * (ph + 10) + 22 : 22;
        const g = svg.append("g").attr("transform", `translate(${ox},${oy})`);
        const max = d3.max(rows, (r) => d3.sum(p.parts, (q) => r[q[0]]));
        const x = d3.scaleLinear().domain([0, max]).range([0, pw - 56]);
        g.append("text").attr("class", "title-text").attr("x", 0).attr("y", -8).text(p.title);
        rows.forEach((r, ri) => {
          const y = ri * rowH + 6;
          if (pi === 0 || narrow) svg.append("text").attr("class", "lab").attr("x", 0).attr("y", oy + y + 13).text(r.name);
          let acc = 0;
          p.parts.forEach(([k, c, lab], qi) => {
            const v = r[k], x0 = x(acc) + (qi ? 1 : 0), x1 = x(acc + v) - (qi < p.parts.length - 1 ? 1 : 0);
            const last = qi === p.parts.length - 1;
            g.append("path").attr("d", last ? C.hBar(x0, x1, y, 18) : `M${x0},${y}H${x1}V${y + 18}H${x0}Z`)
              .attr("fill", color(c)).attr("tabindex", 0).attr("class", "focusable")
              .on("pointermove focus", (e) => showTip(e.type === "focus" ? e.target : e, [{ head: `${r.name}: ${p.title.toLowerCase()}` },
                { key: color(c), value: fmt.int(v), label: lab }]))
              .on("pointerleave blur", hideTip);
            acc += v;
          });
          g.append("text").attr("class", "lab").attr("x", x(acc) + 6).attr("y", y + 13).text(fmt.int(acc));
        });
      });
    };
    C.live(draw);
    tableToggle(document.getElementById("split-card"), () => [["Split", "Genes", "Residues", "Labelled pathogenic", "Labelled benign"],
      ...rows.map((r) => [r.name, fmt.int(r.genes), fmt.int(r.residues), fmt.int(r.path), fmt.int(r.ben)])]);
  }

  TOL.sections = TOL.sections || {};
  TOL.sections.overview = async function () {
    const s = await TOL.data.summary();
    pipeline(s);
    sources(s);
    await splitChart(s);
  };
})();
