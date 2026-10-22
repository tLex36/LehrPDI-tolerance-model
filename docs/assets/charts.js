/* Small chart helpers shared by the sections (built on D3 v7). */
(function () {
  "use strict";
  const TOL = (window.TOL = window.TOL || {});

  /** Path for a horizontal bar from x0 to x1 (either direction) with the data end rounded (r px) and the
   *  baseline end square. */
  function hBar(x0, x1, y, h, r = 4) {
    const dir = x1 >= x0 ? 1 : -1, len = Math.abs(x1 - x0), rr = Math.min(r, len, h / 2);
    if (len < 0.5) return `M${x0},${y}h0v${h}h0Z`;
    const xe = x1 - dir * rr;
    return dir > 0
      ? `M${x0},${y}H${xe}A${rr},${rr} 0 0 1 ${x1},${y + rr}V${y + h - rr}A${rr},${rr} 0 0 1 ${xe},${y + h}H${x0}Z`
      : `M${x0},${y}H${xe}A${rr},${rr} 0 0 0 ${x1},${y + rr}V${y + h - rr}A${rr},${rr} 0 0 0 ${xe},${y + h}H${x0}Z`;
  }

  /** Path for a vertical bar from baseline y0 up/down to y1, data end rounded. */
  function vBar(x, w, y0, y1, r = 4) {
    const dir = y1 <= y0 ? -1 : 1, len = Math.abs(y1 - y0), rr = Math.min(r, len, w / 2);
    if (len < 0.5) return `M${x},${y0}h${w}v0h${-w}Z`;
    const ye = y1 - dir * rr;
    return dir < 0
      ? `M${x},${y0}V${ye}A${rr},${rr} 0 0 1 ${x + rr},${y1}H${x + w - rr}A${rr},${rr} 0 0 1 ${x + w},${ye}V${y0}Z`
      : `M${x},${y0}V${ye}A${rr},${rr} 0 0 0 ${x + rr},${y1}H${x + w - rr}A${rr},${rr} 0 0 0 ${x + w},${ye}V${y0}Z`;
  }

  /** A fresh responsive SVG in host (viewBox = pixel size, so text stays the same size). */
  function svg(host, width, height, label) {
    d3.select(host).selectAll("svg").remove();
    return d3.select(host).append("svg")
      .attr("viewBox", `0 0 ${width} ${height}`).attr("width", width).attr("height", height)
      .attr("role", "img").attr("aria-label", label || "");
  }

  /** Style an axis group and drop its domain line when asked. */
  function axis(g, gen, { domain = true, grid = false, size = 0 } = {}) {
    g.attr("class", grid ? "grid" : "axis").call(gen);
    if (!domain) g.select(".domain").remove();
    if (grid) g.selectAll(".tick text").remove();
    return g;
  }

  /** Register a draw function that re-runs on resize and theme change; returns it after a first draw. */
  function live(draw) {
    TOL.ui.onRedraw(draw);
    draw();
    return draw;
  }

  TOL.charts = { hBar, vBar, svg, axis, live };
})();
