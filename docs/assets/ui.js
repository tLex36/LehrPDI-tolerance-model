/* Shared pieces of the page: theme, colours, tooltip, small DOM helpers, table views, resize handling. */
(function () {
  "use strict";
  const TOL = (window.TOL = window.TOL || {});

  // ---------------------------------------------------------------- DOM helpers
  /** el("div", {class: "x", text: "label"}, child, ...) - text always goes in through textContent. */
  function el(tag, attrs, ...children) {
    const n = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs || {})) {
      if (v === undefined || v === null || v === false) continue;
      if (k === "text") n.textContent = v;
      else if (k === "class") n.className = v;
      else if (k.startsWith("on")) n.addEventListener(k.slice(2), v);
      else n.setAttribute(k, v === true ? "" : v);
    }
    for (const c of children.flat()) {
      if (c === null || c === undefined || c === false) continue;
      n.appendChild(typeof c === "string" || typeof c === "number" ? document.createTextNode(String(c)) : c);
    }
    return n;
  }
  const $ = (sel, root) => (root || document).querySelector(sel);
  const $$ = (sel, root) => Array.from((root || document).querySelectorAll(sel));

  // ---------------------------------------------------------------- number formats
  const fmt = {
    int: (x) => (Number.isFinite(x) ? Math.round(x).toLocaleString("en-US") : "–"),
    f: (x, d = 2) => (Number.isFinite(x) ? x.toFixed(d).replace("-", "−") : "–"),
    signed: (x, d = 2) => (Number.isFinite(x) ? (x > 0 ? "+" : x < 0 ? "−" : "") + Math.abs(x).toFixed(d) : "–"),
    pct: (x, d = 0) => (Number.isFinite(x) ? (100 * x).toFixed(d) + "%" : "–"),
    compact: (x) => (Math.abs(x) >= 1e6 ? (x / 1e6).toFixed(2).replace(/\.?0+$/, "") + "M"
      : Math.abs(x) >= 1e3 ? (x / 1e3).toFixed(1).replace(/\.0$/, "") + "K" : String(x)),
  };

  // ---------------------------------------------------------------- theme
  const THEME_KEY = "tol-theme";
  function storedTheme() { try { return localStorage.getItem(THEME_KEY); } catch (e) { return null; } }
  function applyTheme(t) {
    if (t === "light" || t === "dark") document.documentElement.setAttribute("data-theme", t);
    else document.documentElement.removeAttribute("data-theme");
  }
  const hashTheme = /(?:^#|&)theme=(light|dark)/.exec(location.hash);   // e.g. index.html#theme=light
  applyTheme(hashTheme ? hashTheme[1] : storedTheme());
  function isDark() {
    const t = document.documentElement.getAttribute("data-theme");
    if (t) return t === "dark";
    return window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches;
  }
  const redrawers = new Set();
  /** Charts register a redraw function; it runs on theme change and (debounced) on resize. */
  function onRedraw(fn) { redrawers.add(fn); return fn; }
  function redrawAll() { for (const fn of redrawers) { try { fn(); } catch (e) { console.error(e); } } }
  function setTheme(t) {
    try { if (t) localStorage.setItem(THEME_KEY, t); else localStorage.removeItem(THEME_KEY); } catch (e) { /* private mode */ }
    applyTheme(t);
    redrawAll();
  }
  if (window.matchMedia) {
    const mq = window.matchMedia("(prefers-color-scheme: dark)");
    (mq.addEventListener ? mq.addEventListener.bind(mq, "change") : mq.addListener.bind(mq))(() => redrawAll());
  }
  let resizeTimer = null, lastWidth = window.innerWidth;
  window.addEventListener("resize", () => {
    if (window.innerWidth === lastWidth) return;          // ignore mobile address-bar height changes
    lastWidth = window.innerWidth;
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(redrawAll, 150);
  });

  /** A colour token from the stylesheet (e.g. "--gnomad"), resolved for the current theme. */
  function color(name) {
    return getComputedStyle(document.documentElement).getPropertyValue(name).trim() || "#888";
  }

  // ---------------------------------------------------------------- tooltip (one for the page)
  let tip = null;
  function tooltip() {
    if (!tip) {
      tip = el("div", { class: "tooltip", role: "status", "aria-live": "polite" });
      document.body.appendChild(tip);
    }
    return tip;
  }
  /** rows: [{key: colour or null, value: "0.83", label: "model"} | {head: "title"} | {note: "text"}] */
  function showTip(evtOrRect, rows) {
    const t = tooltip();
    t.replaceChildren(...rows.map((r) => {
      if (r.head) return el("div", { class: "tip-head", text: r.head });
      if (r.note) return el("div", { class: "tip-note", text: r.note });
      return el("div", { class: "tip-row" },
        r.key ? el("span", { class: "tip-key", style: `background:${r.key}` }) : el("span", { class: "tip-key blank" }),
        el("span", { class: "tip-value", text: r.value }),
        el("span", { class: "tip-label", text: r.label || "" }));
    }));
    t.classList.add("on");
    let x, y;
    if (evtOrRect && evtOrRect.clientX !== undefined) { x = evtOrRect.clientX; y = evtOrRect.clientY; }
    else if (evtOrRect && evtOrRect.getBoundingClientRect) {
      const b = evtOrRect.getBoundingClientRect(); x = b.left + b.width / 2; y = b.top;
    } else { x = 0; y = 0; }
    const w = t.offsetWidth, h = t.offsetHeight, pad = 14;
    let left = x + pad, top = y + pad;
    if (left + w > window.innerWidth - 8) left = x - w - pad;
    if (top + h > window.innerHeight - 8) top = y - h - pad;
    t.style.left = Math.max(8, left) + "px";
    t.style.top = Math.max(8, top) + "px";
  }
  function hideTip() { if (tip) tip.classList.remove("on"); }

  // ---------------------------------------------------------------- figure frame with a table view
  /** Wrap a chart host: adds a "Show table" toggle that renders tableFn() (an array of rows, first = header). */
  function tableToggle(figure, tableFn, label) {
    const box = el("div", { class: "table-view", hidden: true });
    const btn = el("button", { class: "btn small ghost", type: "button", "aria-expanded": "false",
      text: label || "Show the numbers as a table" });
    btn.addEventListener("click", () => {
      const open = box.hidden;
      box.hidden = !open;
      btn.setAttribute("aria-expanded", String(open));
      btn.textContent = open ? "Hide the table" : (label || "Show the numbers as a table");
      if (open) fillTable(box, tableFn());
    });
    figure.appendChild(el("div", { class: "fig-foot" }, btn));
    figure.appendChild(box);
    return { refresh: () => { if (!box.hidden) fillTable(box, tableFn()); } };
  }
  function fillTable(box, rows) {
    const [head, ...body] = rows;
    const table = el("table", { class: "data-table" },
      el("thead", {}, el("tr", {}, head.map((h) => el("th", { scope: "col", text: h })))),
      el("tbody", {}, body.map((r) => el("tr", {}, r.map((c, i) => el(i ? "td" : "th", i ? { text: c } : { scope: "row", text: c }))))));
    box.replaceChildren(el("div", { class: "table-scroll" }, table));
  }

  /** Segmented control: options [{value, label}], calls onChange(value). Returns {get, set}. */
  function segmented(host, name, options, value, onChange) {
    const wrap = el("div", { class: "seg", role: "radiogroup", "aria-label": name });
    let current = value;
    const buttons = options.map((o) => {
      const b = el("button", { type: "button", role: "radio", "aria-checked": String(o.value === value),
        class: o.value === value ? "on" : "", text: o.label, title: o.title || null });
      b.addEventListener("click", () => set(o.value, true));
      wrap.appendChild(b);
      return b;
    });
    function set(v, fire) {
      current = v;
      buttons.forEach((b, i) => {
        const on = options[i].value === v;
        b.classList.toggle("on", on); b.setAttribute("aria-checked", String(on));
      });
      if (fire) onChange(v);
    }
    host.appendChild(wrap);
    return { get: () => current, set };
  }

  /** Width available to a chart host. */
  const widthOf = (node) => Math.max(280, Math.floor(node.getBoundingClientRect().width));

  /** Read "#key=value&..." from the address bar (used for shareable links and for testing). */
  function hashParams() {
    const out = {};
    for (const part of location.hash.replace(/^#/, "").split("&")) {
      if (!part.includes("=")) continue;
      const [k, v] = part.split("=");
      out[decodeURIComponent(k)] = decodeURIComponent(v);
    }
    return out;
  }

  TOL.ui = { el, $, $$, fmt, color, isDark, setTheme, storedTheme, onRedraw, redrawAll, showTip, hideTip,
             tableToggle, segmented, widthOf, hashParams };
})();
