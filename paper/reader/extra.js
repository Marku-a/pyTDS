/* Reader-only figures. Uses frame/path/el/hover/hitDot/legend/css/fmt from the paper's chart code. */
const ORDER_B = ["LS", "awake", "REM", "DS"];   // Bashan: most → least connected
function meanCI(a) { const x = a.filter(v => v != null && isFinite(v)); const n = x.length; if (!n) return [NaN, NaN, NaN, 0];
  const m = x.reduce((s, v) => s + v, 0) / n, sd = Math.sqrt(x.reduce((s, v) => s + (v - m) ** 2, 0) / Math.max(1, n - 1));
  const h = 1.96 * sd / Math.sqrt(n); return [m, m - h, m + h, n]; }

function fig12() {
  [["published", "fig12a", "--s1", "TDS"], ["v2_calibrated", "fig12b", "--s2", "aTDS"]].forEach(([v, id, col, name]) => {
    const C = D.connectivity[v].strict; const all = ORDER_B.flatMap(s => C[s].filter(x => x != null));
    const hi = Math.ceil(Math.max(...all) / 5) * 5, lo = Math.min(0, Math.floor(Math.min(...all)));
    const f = frame(document.getElementById(id), { xband: true, y: [lo, hi], height: 270, ylabel: "connectivity (pts above chance)" });
    const bw = (f.W - f.m.l - f.m.r) / ORDER_B.length;
    ORDER_B.forEach((s, k) => { const cx = f.m.l + k * bw + bw / 2, [m, a, b, n] = meanCI(C[s]);
      const t = el("text", { x: cx, y: f.H - f.m.b + 16, "text-anchor": "middle" }, f.svg); t.textContent = STAGE_NAME[s];
      C[s].forEach((x, i) => { if (x == null) return; const jx = cx - bw * 0.22 + ((i * 37) % 100) / 100 * bw * 0.16;
        el("circle", { cx: jx, cy: f.y(x), r: 2.2, fill: css(col), "fill-opacity": .35 }, f.g); });
      el("line", { x1: cx + bw * 0.08, x2: cx + bw * 0.08, y1: f.y(a), y2: f.y(b), stroke: css(col), "stroke-width": 2.5 }, f.g);
      el("circle", { cx: cx + bw * 0.08, cy: f.y(m), r: 6, fill: css(col), stroke: css("--paper"), "stroke-width": 2 }, f.g);
      hitDot(f.g, cx + bw * 0.08, f.y(m), 8, () => `<b>${name}, ${STAGE_NAME[s]}</b><br>mean ${fmt(m, 1)} pts [${fmt(a, 1)}, ${fmt(b, 1)}]<br>${n} of 35 people scored`); });
    if (lo < 0) el("line", { x1: f.m.l, x2: f.W - f.m.r, y1: f.y(0), y2: f.y(0), stroke: css("--ink-3") }, f.g); });
  legend("leg12", [{ label: "TDS", color: "--s1", kind: "dot" }, { label: "aTDS", color: "--s2", kind: "dot" }, { label: "large dot = mean, bar = 95 % CI, small dots = people", color: "--ink-3", kind: "dot" }]);
}

function linePanel(id, xs, series, o) {
  const f = frame(document.getElementById(id), Object.assign({ height: 250 }, o));
  if (o.hline != null) { el("line", { x1: f.m.l, x2: f.W - f.m.r, y1: f.y(o.hline), y2: f.y(o.hline), stroke: css("--ink-3"), "stroke-dasharray": "2 4" }, f.g);
    const t = el("text", { x: f.W - f.m.r - 2, y: f.y(o.hline) - 4, "text-anchor": "end" }, f.svg); t.textContent = o.hlabel || ""; }
  series.forEach(s => { const pts = xs.map((x, i) => [f.x(x), s.y[i] == null ? null : f.y(s.y[i])]);
    path(f.g, pts, css(s.col), s.dash ? 1.8 : 2.2, s.dash ? "5 4" : "");
    pts.forEach((p, i) => { if (p[1] == null) return; el("circle", { cx: p[0], cy: p[1], r: s.dash ? 3 : 4.5, fill: s.dash ? css("--paper") : css(s.col), stroke: css(s.col), "stroke-width": 1.6 }, f.g);
      hitDot(f.g, p[0], p[1], 5, () => `<b>${s.name}</b><br>${o.xname} ${xs[i]}: ${o.fmtv(s.y[i])}${s.extra ? "<br>" + s.extra(i) : ""}`); }); });
  return f;
}
const KEY3 = [{ label: "TDS (fixed default)", color: "--s1" }, { label: "aTDS", color: "--s2" }, { label: "oracle (best, found with the answer)", color: "--null" }];

function fig21() {
  const E = R.synthetic.E1_weak_coupling_vs_slowness, xs = E.fixed_default.r;
  const mk = (key, scale) => [["fixed_default", "TDS", "--s1"], ["v2_calibrated", "aTDS", "--s2"], ["oracle", "oracle", "--null", true]].map(([v, name, col, dash]) =>
    ({ name, col, dash, y: E[v][key].map(x => x * scale), extra: i => `window ${E[v].window[i]}, tolerance ±${E[v].tolerance[i]}` }));
  linePanel("fig21a", xs, mk("auc", 1), { x: [0.8, 20], xlog: true, xticks: xs, y: [0, 1], yticks: [0, 0.25, 0.5, 0.75, 1], hline: 0.5, hlabel: "chance", xlabel: "signals slower by ×", ylabel: "AUC", xname: "×", fmtv: v => fmt(v, 2) });
  linePanel("fig21b", xs, mk("hit_rate", 100), { x: [0.8, 20], xlog: true, xticks: xs, y: [0, 100], xlabel: "signals slower by ×", ylabel: "% correct delay", xname: "×", fmtv: v => fmt(v, 0) + " %" });
  legend("leg21", KEY3);
}
function fig22() {
  const E = R.synthetic.E2_rulkov_vs_delay, xs = E.fixed_default.delay;
  const mk = (key, scale) => [["fixed_default", "TDS", "--s1"], ["v2_calibrated", "aTDS", "--s2"], ["oracle", "oracle", "--null", true]].map(([v, name, col, dash]) =>
    ({ name, col, dash, y: E[v][key].map(x => x * scale), extra: i => `window ${E[v].window[i]}, tolerance ±${E[v].tolerance[i]}` }));
  linePanel("fig22a", xs, mk("auc", 1), { x: [6, 40], xlog: true, xticks: xs, y: [0, 1], yticks: [0, 0.25, 0.5, 0.75, 1], hline: 0.5, hlabel: "chance", xlabel: "true delay (steps)", ylabel: "AUC", xname: "delay", fmtv: v => fmt(v, 2) });
  linePanel("fig22b", xs, mk("hit_rate", 100), { x: [6, 40], xlog: true, xticks: xs, y: [0, 100], xlabel: "true delay (steps)", ylabel: "% correct delay", xname: "delay", fmtv: v => fmt(v, 0) + " %" });
  legend("leg22", KEY3);
}
function tab23() {
  const S = R.synthetic, g = (arr, cond, v, k) => (arr.find(z => z.condition === cond && z.variant === v) || {})[k];
  const conds = S.E1c_mixed_speeds.map(z => z.condition).filter((c, i, a) => a.indexOf(c) === i);
  const worst = v => Math.min(...conds.map(c => g(S.E1c_mixed_speeds, c, v, "auc")));
  const rows = [["Signals of different speeds, worst of 4 combinations (AUC)", worst("fixed_default"), worst("v2_calibrated"), worst("oracle"), "aTDS"],
    ["Coupling switching on/off, fast signals (balanced accuracy)", g(S.E1d_on_off, "r=1", "fixed_default", "balanced_acc_mean"), g(S.E1d_on_off, "r=1", "v2_calibrated", "balanced_acc_mean"), null, "tie"],
    ["Coupling switching on/off, 4× slower (balanced accuracy)", g(S.E1d_on_off, "r=4", "fixed_default", "balanced_acc_mean"), g(S.E1d_on_off, "r=4", "v2_calibrated", "balanced_acc_mean"), null, "aTDS"]];
  document.getElementById("t23").innerHTML = rows.map(r => `<tr><td class="l">${r[0]}</td><td>${fmt(r[1])}</td><td>${fmt(r[2])}</td><td>${r[3] == null ? "–" : fmt(r[3])}</td><td class="l">${pill(r[4])}</td></tr>`).join("");
}
function fig24() {
  const E = R.synthetic.R4_sampling_rate, xs = E.fixed_default.factor;
  linePanel("fig24", xs, [["fixed_default", "TDS", "--s1"], ["v2_calibrated", "aTDS (no cap)", "--s2"]].map(([v, name, col]) =>
    ({ name, col, y: E[v].score, extra: i => `window ${E[v].window[i]} samples, tolerance ±${E[v].tolerance[i]}` })),
    { x: [0.4, 5], xlog: true, xticks: xs, y: [0, 100], xlabel: "sampling rate (× original 1 Hz)", ylabel: "TDS score (%)", xname: "rate ×", fmtv: v => fmt(v, 1) + " %", height: 260 });
  legend("leg24", [{ label: "TDS (60 samples / ±1 sample)", color: "--s1" }, { label: "aTDS, no window cap (saturates at 100 %)", color: "--s2" }]);
}
function forest(id, rows, o) {
  const g = frame(document.getElementById(id), { x: o.x, y: [0, rows.length], height: 34 + rows.length * 24 + 30, xlabel: o.xlabel, margin: { l: o.ml || 110 }, yticks: [] });
  rows.forEach((r, i) => { const yy = g.y(rows.length - i - 0.5);
    if (r.label) { const t = el("text", { x: g.m.l - 6, y: yy + 4, "text-anchor": "end" }, g.svg); t.textContent = r.label; }
    el("line", { x1: g.x(Math.max(o.x[0], r.lo)), x2: g.x(Math.min(o.x[1], r.hi)), y1: yy, y2: yy, stroke: css(r.col), "stroke-width": 2.2 }, g.g);
    el("circle", { cx: g.x(r.est), cy: yy, r: 5, fill: css(r.col) }, g.g);
    hitDot(g.g, g.x(r.est), yy, 6, () => r.tip); });
  if (o.x[0] < 0) el("line", { x1: g.x(0), x2: g.x(0), y1: g.m.t, y2: g.H - g.m.b, stroke: css("--ink-3") }, g.g);
}
function fig25() {
  const E = R.groundtruth.effects.filter(z => z.dataset === "sleep_edf"), effs = ["LS>DS dz", "REM>DS dz", "W>DS dz"], nm = { "LS>DS dz": "light > deep", "REM>DS dz": "REM > deep", "W>DS dz": "wake > deep" };
  const rows = []; effs.forEach(e => [["published", "TDS", "--s1"], ["v2_calibrated", "aTDS", "--s2"]].forEach(([v, name, col], j) => {
    const z = E.find(q => q.effect === e && q.variant === v);
    rows.push({ label: j ? "" : nm[e], est: z.estimate, lo: z.ci_lo, hi: z.ci_hi, col, tip: `<b>${name}, ${nm[e]}</b><br>dz ${fmt(z.estimate)} [${fmt(z.ci_lo)}, ${fmt(z.ci_hi)}]${v !== "published" ? `<br>gain ${fmt(z.diff_vs_published)} [${fmt(z.diff_ci_lo)}, ${fmt(z.diff_ci_hi)}]` : ""}` }); }));
  forest("fig25", rows, { x: [-0.2, 6.2], xlabel: "dz (paired effect size)" });
  legend("leg25", [{ label: "TDS (published)", color: "--s1", kind: "dot" }, { label: "aTDS", color: "--s2", kind: "dot" }]);
}
function fig26() {
  const E = R.groundtruth.effects.filter(z => z.dataset === "fantasia");
  const rows = [["published", "TDS", "--s1"], ["v2_calibrated", "aTDS", "--s2"]].map(([v, name, col]) => { const z = E.find(q => q.variant === v);
    return { label: name, est: z.estimate, lo: z.ci_lo, hi: z.ci_hi, col, tip: `<b>${name}</b><br>d ${fmt(z.estimate)} [${fmt(z.ci_lo)}, ${fmt(z.ci_hi)}]` }; });
  forest("fig26a", rows, { x: [0, 2.6], xlabel: "d (young > elderly)", ml: 60 });
  legend("leg26", [{ label: "TDS", color: "--s1", kind: "dot" }, { label: "aTDS", color: "--s2", kind: "dot" }]);
  const F = R.groundtruth.fantasia; const f = frame(document.getElementById("fig26b"), { x: [0, 10], y: [0, 2], height: 130, xlabel: "false links (of 40 cross-person pairs)", margin: { l: 60 }, yticks: [], xticks: [0, 2, 4, 6, 8, 10] });
  [["published", "TDS", "--s1"], ["v2_calibrated", "aTDS", "--s2"]].forEach(([v, name, col], i) => { const z = F.find(q => q.variant === v), n = Math.round(z.xsubj_fpr * 40), yy = f.y(2 - i - 0.5);
    const t = el("text", { x: f.m.l - 6, y: yy + 4, "text-anchor": "end" }, f.svg); t.textContent = name;
    el("rect", { x: f.x(0), y: yy - 9, width: Math.max(2, f.x(n) - f.x(0)), height: 18, rx: 3, fill: css(col) }, f.g);
    const vt = el("text", { x: f.x(n) + 5, y: yy + 4 }, f.svg); vt.textContent = `${n}/40`; });
}
const FIELD_NAME = { hydrology_hourly: "Rain → river flow (hourly)", hydrology_daily: "Rain → river flow (daily)", energy_vic_summer: "Heat → electricity demand", beijing_pm25: "Wind → PM2.5 pollution",
  climate_nino34_gistemp: "El Niño → global temperature", gas_furnace: "Gas feed → CO₂ (296 samples)", gasoline_weekly_diff: "Crude → petrol price", jena_T_rh: "Temperature ↔ humidity (control)",
  sleep_eeg_delta_sigma: "Sleep EEG δ ↔ σ", neurokit_hr_rsp: "Heart rate ↔ breathing" };
const FIELD_VERDICT = { hydrology_hourly: "aTDS", hydrology_daily: "tie", energy_vic_summer: "aTDS", beijing_pm25: "aTDS", climate_nino34_gistemp: "mixed", gas_furnace: "TDS",
  gasoline_weekly_diff: "tie", jena_T_rh: "tie", sleep_eeg_delta_sigma: "aTDS", neurokit_hr_rsp: "tie" };
function fig31() {
  const F = R.fields, ds = Object.keys(FIELD_NAME);
  const f = frame(document.getElementById("fig31"), { x: [-10, 100], y: [0, ds.length], height: 40 + ds.length * 30 + 30, xlabel: "score minus chance (pts)", margin: { l: 196, r: 90 }, yticks: [], xticks: [0, 20, 40, 60, 80, 100] });
  el("line", { x1: f.x(0), x2: f.x(0), y1: f.m.t, y2: f.H - f.m.b, stroke: css("--ink-3") }, f.g);
  ds.forEach((d, i) => { const yy = f.y(ds.length - i - 0.5);
    const t = el("text", { x: f.m.l - 6, y: yy + 4, "text-anchor": "end" }, f.svg); t.textContent = FIELD_NAME[d];
    const a = F.find(z => z.dataset === d && z.variant === "fixed_default"), b = F.find(z => z.dataset === d && z.variant === "adaptive");
    const ea = a.score - a.null95, eb = b.score - b.null95, mk = z => z.match === true ? " ✓" : z.match === false ? " ✗" : "";
    const right = Math.max(ea, eb) < 75, dl = el("text", { x: right ? f.x(Math.max(ea, eb)) + 10 : f.x(Math.min(ea, eb)) - 10, y: yy + 4, "text-anchor": right ? "start" : "end" }, f.svg);
    dl.textContent = `${a.tau_phys}${mk(a)} / ${b.tau_phys}${mk(b)}`; dl.style.fontSize = "10.5px";
    el("line", { x1: f.x(ea), x2: f.x(eb), y1: yy, y2: yy, stroke: css("--rule"), "stroke-width": 3 }, f.g);
    [[a, ea, "--s1", "TDS"], [b, eb, "--s2", "aTDS"]].forEach(([z, e, col, name]) => { el("circle", { cx: f.x(e), cy: yy, r: 6, fill: css(col), stroke: css("--paper"), "stroke-width": 1.5 }, f.g);
      hitDot(f.g, f.x(e), yy, 7, () => `<b>${name}: ${FIELD_NAME[d]}</b><br>score ${fmt(z.score, 1)} % · chance ${fmt(z.null95, 1)} %<br>delay ${z.tau_phys}${z.match === true ? " ✓" : z.match === false ? " ✗" : ""}<br>window ${z.window_phys}, tolerance ±${z.tolerance_phys}${z.calibration_failed === true ? "<br>calibration warning" : ""}`); });
    const v = el("text", { x: f.W - f.m.r + 8, y: yy + 4 }, f.svg); v.textContent = { aTDS: "aTDS better", TDS: "TDS better", tie: "tie", mixed: "mixed" }[FIELD_VERDICT[d]];
    v.style.fill = css(FIELD_VERDICT[d] === "aTDS" ? "--s2" : FIELD_VERDICT[d] === "TDS" ? "--s1" : "--ink-3"); v.style.fontWeight = 600; });
  legend("leg31", [{ label: "TDS (fixed default)", color: "--s1", kind: "dot" }, { label: "aTDS", color: "--s2", kind: "dot" }]);
}
function pill(w) { const map = { aTDS: ["aTDS", "--s2"], TDS: ["TDS", "--s1"], tie: ["tie", "--ink-3"], mixed: ["mixed", "--ink-3"], none: ["neither", "--ink-3"], na: ["no comparison", "--ink-3"] };
  const [t, c] = map[w] || [w, "--ink-3"]; return `<span class="status" style="color:var(${c})">${t}</span>`; }
const SCORE = [
  ["Bashan-2012 cohort, 35 people (Part 1)", "yes (stage ordering)", "known stage differences stand out more", "tie", "same ordering; dz not sharper (inconclusive)"],
  ["Bashan-2012 cohort, 35 people (Part 1)", "yes", "real coupling separated from chance (AUC)", "aTDS", "0.778 vs 0.736; loses wake in 25/35 people at 300 s"],
  ["Synthetic, fast signals (2.1)", "yes", "AUC, correct delay", "tie", "both perfect"],
  ["Synthetic, 4–16× slower (2.1)", "yes", "AUC, correct delay", "aTDS", "8×: 0.91 / 95 % vs 0.35 / 0 %; oracle 1.00"],
  ["Mixed speeds, on/off coupling (2.3)", "yes", "AUC, accuracy", "aTDS", "worst case 1.00 vs 0.76"],
  ["Rulkov neurons, delay 32 (2.2)", "yes", "correct delay", "aTDS", "100 % vs 0 % (fixed search too short)"],
  ["Rulkov neurons, delay 16 (2.2)", "yes", "correct delay", "TDS", "20 % vs 0 %; delay near half a burst period"],
  ["Same EEG at 4 sampling rates (2.4)", "yes (should not change)", "same answer at every rate", "aTDS", "fixed 75 → 27 %; aTDS 100 % (saturated, no cap)"],
  ["Sleep-EDF, 18 nights (2.5)", "yes", "known stage differences", "aTDS", "LS>DS +1.30, REM>DS +0.72; W>DS −0.23 (n.s.)"],
  ["Fantasia, 40 people (2.6)", "yes", "young vs elderly; false links", "tie", "d 1.54 vs 1.25 (n.s.); false links 1/40 vs 4/40"],
  ["ECG → finger pulse (2.7)", "yes (430 ms)", "correct delay", "none", "both report −80 ms: half-beat limit"],
  ["Heart rate ↔ breathing (2.7)", "yes", "score vs chance", "tie", "aTDS chose ~the default by itself"],
  ["Treadmill CPET (2.7)", "yes", "speed → HR detection", "na", "design flaw: uninformative"],
  ["Nine fields (3.1)", "roughly", "above chance, delay", "aTDS", "4 better, 4 tie/mixed, 1 worse (296 samples)"],
  ["Garmin, 30 + 4 runs (3.2)", "roughly (HR follows speed)", "speed → HR detected", "aTDS", "8 vs 4 (p = 0.29); held-out 2/4 vs 0/4"],
  ["Thesis apnea cohort (Appendix)", "–", "–", "na", "inputs unavailable; only TDS outputs"]];
function tabScore() { document.getElementById("tscore").innerHTML = SCORE.map(r => `<tr><td class="l">${r[0]}</td><td class="l">${r[1]}</td><td class="l">${r[2]}</td><td class="l">${pill(r[3])}</td><td class="l">${r[4]}</td></tr>`).join(""); }
