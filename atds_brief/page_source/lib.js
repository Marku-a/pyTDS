const NS = "http://www.w3.org/2000/svg";
function el(tag, attrs, parent) { const e = document.createElementNS(NS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; }
function scale(d0, d1, r0, r1, log) {
  if (log) { const a = Math.log(d0), b = Math.log(d1); return v => r0 + (Math.log(v) - a) / (b - a) * (r1 - r0); }
  return v => r0 + (v - d0) / (d1 - d0) * (r1 - r0);
}
function fmt(v) { return typeof v === "number" ? (Math.abs(v) >= 100 ? v.toLocaleString("en-US") : String(v)) : v; }

// Builds one SVG plot with axes and a nearest-point tooltip layer.
function plot(container, o) {
  const W = o.w || 680, H = o.h || 300, m = Object.assign({ l: 58, r: 20, t: o.title ? 28 : 14, b: 48 }, o.m || {});
  const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": o.aria || "" }, container);
  const x = scale(o.x.min, o.x.max, m.l, W - m.r, o.x.log), y = scale(o.y.min, o.y.max, H - m.b, m.t, o.y.log);
  const ax = el("g", { class: "axis" }, svg);
  if (o.title) el("text", { x: m.l, y: 16, class: "panel-title" }, svg).textContent = o.title;
  o.y.ticks.forEach(t => {
    el("line", { x1: m.l, x2: W - m.r, y1: y(t), y2: y(t), class: "gridline" }, ax);
    el("text", { x: m.l - 8, y: y(t) + 4, "text-anchor": "end" }, ax).textContent = (o.y.fmt || fmt)(t);
  });
  el("line", { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, class: "baseline" }, ax);
  o.x.ticks.forEach(t => {
    el("line", { x1: x(t), x2: x(t), y1: H - m.b, y2: H - m.b + 4, class: "baseline" }, ax);
    el("text", { x: x(t), y: H - m.b + 18, "text-anchor": "middle" }, ax).textContent = (o.x.fmt || fmt)(t);
  });
  el("text", { x: (m.l + W - m.r) / 2, y: H - 8, "text-anchor": "middle", class: "axis-label" }, svg).textContent = o.x.label;
  const yl = el("text", { x: 0, y: 0, "text-anchor": "middle", class: "axis-label", transform: `translate(14 ${(m.t + H - m.b) / 2}) rotate(-90)` }, svg);
  yl.textContent = o.y.label;
  const layer = el("g", {}, svg);
  const pts = [];
  return { svg, x, y, layer, pts, W, H, m, container };
}
function attachTip(p) {
  const tip = document.createElement("div"); tip.className = "tip"; tip.hidden = true; p.container.appendChild(tip);
  const ring = el("circle", { r: 7, class: "hover-ring", visibility: "hidden" }, p.svg);
  function move(ev) {
    const r = p.svg.getBoundingClientRect(), sx = p.W / r.width;
    const mx = (ev.clientX - r.left) * sx, my = (ev.clientY - r.top) * sx;
    let best = null, bd = 1e9;
    p.pts.forEach(q => { const d = Math.hypot(q.px - mx, q.py - my); if (d < bd) { bd = d; best = q; } });
    if (!best || bd > 40) { tip.hidden = true; ring.setAttribute("visibility", "hidden"); return; }
    ring.setAttribute("cx", best.px); ring.setAttribute("cy", best.py); ring.setAttribute("visibility", "visible");
    const cr = p.container.getBoundingClientRect();
    tip.innerHTML = best.html; tip.hidden = false;
    let left = r.left - cr.left + best.px / sx; left = Math.max(80, Math.min(cr.width - 80, left));
    tip.style.left = left + "px"; tip.style.top = (r.top - cr.top + best.py / sx) + "px";
  }
  p.svg.addEventListener("pointermove", move);
  p.svg.addEventListener("pointerleave", () => { tip.hidden = true; ring.setAttribute("visibility", "hidden"); });
}
function line(p, data, color, opts = {}) {
  const d = data.map((v, i) => (i ? "L" : "M") + p.x(v[0]).toFixed(1) + " " + p.y(v[1]).toFixed(1)).join(" ");
  el("path", { d, fill: "none", stroke: color, "stroke-width": 2, "stroke-linejoin": "round", "stroke-dasharray": opts.dash || "none" }, p.layer);
  if (opts.dots !== false) data.forEach(v => el("circle", { cx: p.x(v[0]), cy: p.y(v[1]), r: 4, fill: color, stroke: "var(--surface)", "stroke-width": 2 }, p.layer));
  data.forEach(v => p.pts.push({ px: p.x(v[0]), py: p.y(v[1]), html: opts.tip(v) }));
}
function legend(container, items) {
  const L = document.createElement("div"); L.className = "legend";
  L.innerHTML = items.map(i => `<span><i class="sw ${i.k || ""}" style="background:${i.c};border-color:${i.c}"></i>${i.n}</span>`).join("");
  container.prepend(L);
}
function table(id, head, rows) {
  document.getElementById(id).innerHTML = "<table><thead><tr>" + head.map(h => `<th>${h}</th>`).join("") + "</tr></thead><tbody>" +
    rows.map(r => "<tr>" + r.map(c => `<td>${c}</td>`).join("") + "</tr>").join("") + "</tbody></table>";
}

