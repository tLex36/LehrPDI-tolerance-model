/* Start-up: theme button, section navigation, then each section in turn. */
(function () {
  "use strict";
  const TOL = window.TOL;
  const { el } = TOL.ui;

  function fail(message) {
    const main = document.getElementById("main");
    main.prepend(el("div", { class: "callout warn", role: "alert" }, el("strong", { text: "Something did not load" }), el("p", { text: message })));
    document.documentElement.setAttribute("data-js-error", message);
  }
  window.addEventListener("error", (e) => document.documentElement.setAttribute("data-js-error", String(e.message)));
  window.addEventListener("unhandledrejection", (e) => document.documentElement.setAttribute("data-js-error", String(e.reason && e.reason.message || e.reason)));

  // theme button: cycles to the opposite of what is showing
  const btn = document.getElementById("theme-btn");
  const label = () => { btn.textContent = TOL.ui.isDark() ? "Light" : "Dark"; };
  btn.addEventListener("click", () => { TOL.ui.setTheme(TOL.ui.isDark() ? "light" : "dark"); label(); });
  label();

  // highlight the section in view
  const links = Array.from(document.querySelectorAll("nav.toc a"));
  if ("IntersectionObserver" in window) {
    const byId = new Map(links.map((a) => [a.getAttribute("href").slice(1), a]));
    const io = new IntersectionObserver((entries) => {
      entries.forEach((en) => {
        if (!en.isIntersecting) return;
        links.forEach((a) => a.classList.remove("active"));
        const a = byId.get(en.target.id);
        if (a) {
          a.classList.add("active");
          const ol = a.closest("ol");                         // keep the active link visible in the scrolling nav
          if (ol && (a.offsetLeft < ol.scrollLeft || a.offsetLeft + a.offsetWidth > ol.scrollLeft + ol.clientWidth)) ol.scrollLeft = a.offsetLeft - 16;
        }
      });
    }, { rootMargin: "-45% 0px -50% 0px" });
    document.querySelectorAll("main section[id]").forEach((s) => io.observe(s));
  }

  async function start() {
    if (typeof d3 === "undefined") {
      fail("The charting library (D3) could not be loaded. It comes from cdn.jsdelivr.net, so this page needs an internet connection.");
      return;
    }
    const order = ["overview", "explorer", "experiments", "test", "inside", "conservation"];
    for (const name of order) {
      try { await TOL.sections[name](); }
      catch (e) { console.error(e); fail(`The "${name}" section failed: ${e.message}`); }
    }
    document.documentElement.setAttribute("data-ready", "true");
  }
  start();
})();
