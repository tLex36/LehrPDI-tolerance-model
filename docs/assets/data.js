/* Loading the data files in docs/data/.
 *
 * Each data file is JSON wrapped in a function call (TOL.gene({...}); TOL.index({...}); ...). Loading it as a
 * <script> calls that function, which hands the data to whoever asked for it. Unlike fetch(), this also works
 * when index.html is opened straight from disk (file://). The files are written by tools/explorer/build_data.py
 * and tools/explorer/build_model_data.py.
 */
(function () {
  "use strict";
  const TOL = (window.TOL = window.TOL || {});
  const waiting = new Map();      // key -> {promise, resolve, reject}
  const cache = new Map();        // key -> data

  function deliver(key, data) {
    cache.set(key, data);
    const w = waiting.get(key);
    if (w) { w.resolve(data); waiting.delete(key); }
  }
  // the calls the data files make
  TOL.gene = (d) => deliver("gene:" + d.gene, d);
  TOL.model = (d) => deliver("model:" + d.gene, d);
  TOL.index = (d) => deliver("index", d);
  TOL.test = (d) => deliver("test", d);
  TOL.summary = (d) => deliver("summary", d);
  TOL.modelSummary = (d) => deliver("modelSummary", d);

  /** Load data/<path>.js once; resolves with what the file delivers under `key`. */
  function load(key, path) {
    if (cache.has(key)) return Promise.resolve(cache.get(key));
    if (waiting.has(key)) return waiting.get(key).promise;
    let resolve, reject;
    const promise = new Promise((res, rej) => { resolve = res; reject = rej; });
    waiting.set(key, { promise, resolve, reject });
    const s = document.createElement("script");
    s.src = "data/" + path + ".js";
    s.async = true;
    s.onerror = () => { waiting.delete(key); reject(new Error("could not load data/" + path + ".js")); };
    document.head.appendChild(s);
    return promise;
  }

  // per-residue values are stored as comma-separated integers: value = integer / scale; empty = missing
  const SCALES = { dis: 100, bind: 100, rsa: 100, plddt: 1, esm: 10, cons: 10 };
  function decodeInts(text, scale) {
    const parts = text.split(","), out = new Float64Array(parts.length);
    for (let i = 0; i < parts.length; i++) out[i] = parts[i] === "" ? NaN : +parts[i] / scale;
    return out;
  }

  const TOKEN = /^([A-Z])(\d+)([A-Z*])([usxpbvo])(\d?)$/;
  /** "E154Du" / "E154Dv1" -> {ref: "E", pos: 154, alt: "D", flag: "u" | class, stars} */
  function decodeVariants(text) {
    if (!text) return [];
    return text.split(" ").map((t) => {
      const m = TOKEN.exec(t);
      return { ref: m[1], pos: +m[2], alt: m[3], flag: m[4], stars: m[5] === "" ? null : +m[5] };
    });
  }

  /** A gene file, decoded into typed arrays and variant lists. */
  async function gene(name) {
    const key = "gene:" + name;
    if (cache.has(key) && cache.get(key).decoded) return cache.get(key).decoded;
    const raw = await load(key, "genes/" + encodeURIComponent(name));
    const d = {
      gene: raw.gene, split: raw.split, uniprot: raw.uniprot, seq: raw.seq, L: raw.seq.length,
      ann: raw.ann || {}, check: raw.check, test: raw.test || [],
      gnomad: decodeVariants(raw.gn), clinvar: decodeVariants(raw.cv),
    };
    for (const k of Object.keys(SCALES)) d[k] = decodeInts(raw[k], SCALES[k]);
    raw.decoded = d;
    return d;
  }

  async function index() {
    const raw = await load("index", "genes");
    if (!raw.byName) {
      raw.list = raw.rows.map((r) => Object.fromEntries(raw.columns.map((c, i) => [c, r[i]])));
      raw.byName = new Map(raw.list.map((g) => [g.gene, g]));
    }
    return raw;
  }

  async function test() {
    const raw = await load("test", "test");
    if (!raw.list) raw.list = raw.rows.map((r) => Object.fromEntries(raw.columns.map((c, i) => [c, r[i]])));
    return raw;
  }

  TOL.data = {
    load, gene, index, test, SCALES,
    summary: () => load("summary", "summary"),
    modelSummary: () => load("modelSummary", "model"),
    model: (name) => load("model:" + name, "model/" + encodeURIComponent(name)),
  };
})();
