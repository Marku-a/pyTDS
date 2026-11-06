// Looping step animations for the algorithm section
(() => {
  const root = document.documentElement;
  const css = n => getComputedStyle(root).getPropertyValue(n).trim();
  const reduce = matchMedia("(prefers-reduced-motion: reduce)").matches;
  let seed = 11;
  const rnd = () => { seed |= 0; seed = seed + 0x6D2B79F5 | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; };
  const gauss = () => Math.sqrt(-2 * Math.log(rnd() || 1e-9)) * Math.cos(2 * Math.PI * rnd());
  function ar1(n, tc) { const phi = Math.exp(-1 / tc); let x = 0; const o = []; for (let i = 0; i < n + 300; i++) { x = phi * x + gauss(); if (i >= 300) o.push(x); } const m = o.reduce((a, b) => a + b) / n, s = Math.sqrt(o.reduce((a, b) => a + (b - m) ** 2, 0) / n); return o.map(v => (v - m) / s); }
  const N = 12000, DEL = 20, X = ar1(N, 8), E = ar1(N, 2);
  const Y = X.map((_, i) => 0.6 * X[(i - DEL + N) % N] + 0.8 * E[i]);
  function bart(S_) { const r = []; for (let k = 0; k <= 400; k++) { let s = 0; for (let t = 0; t + k < N; t++) s += S_[t] * S_[t + k]; r.push(s); }
    const r0 = r[0]; for (let k = 0; k < r.length; k++) r[k] /= r0;
    let K_ = r.findIndex(v => v <= 0); if (K_ < 0) K_ = r.length;
    const Bc_ = [1]; for (let k = 1; k < r.length; k++) Bc_.push(Bc_[k - 1] + (k < K_ ? 2 * r[k] ** 2 : 0));
    return { r, K: K_, Bc: Bc_, B: Bc_[Bc_.length - 1] }; }
  const AX = bart(X), AY = bart(Y);
  const B = Math.max(AX.B, AY.B), L = Math.round(30 * B), LAG = Math.floor((L - 1) / 2), STEP = Math.floor(L / 2);
  const font = (s, w) => `${w || 400} ${s}px "IBM Plex Sans", system-ui, sans-serif`;
  const mono = s => `${s}px "IBM Plex Mono", ui-monospace, monospace`;
  const clamp = (v, a, b) => Math.max(a, Math.min(b, v));
  const ease = t => t < 0.5 ? 2 * t * t : 1 - (-2 * t + 2) ** 2 / 2;
  function trace(ctx, arr, x0, x1, yMid, amp, col, lw, xs) {
    ctx.beginPath(); ctx.strokeStyle = col; ctx.lineWidth = lw || 1.5;
    arr.forEach((v, i) => { const x = xs ? xs(i) : x0 + i / (arr.length - 1) * (x1 - x0), y = yMid - v * amp; if (x == null) { ctx.stroke(); ctx.beginPath(); return; } i && xs == null ? ctx.lineTo(x, y) : ctx.moveTo(x, y); if (xs && i) ctx.lineTo(x, y); });
    ctx.stroke();
  }
  function label(ctx, txt, x, y, col, f, align) { ctx.fillStyle = col; ctx.font = f || font(12); ctx.textAlign = align || "left"; ctx.fillText(txt, x, y); }

  function anim(id, dur, draw, HH) {
    const cv = document.getElementById(id); if (!cv) return;
    const ctx = cv.getContext("2d"), W = 680, H = HH || 240;
    const d = window.devicePixelRatio || 1; cv.width = W * d; cv.height = H * d; ctx.setTransform(d, 0, 0, d, 0, 0);
    let paused = reduce, t = reduce ? dur * 0.999 : 0, prev = null, vis = true;
    const btn = document.querySelector(`[data-anim="${id}"]`);
    const setBtn = () => { if (btn) { btn.textContent = paused ? "Play" : "Pause"; btn.setAttribute("aria-pressed", String(paused)); } };
    if (btn) btn.addEventListener("click", () => { paused = !paused; prev = null; setBtn(); });
    setBtn();
    if ("IntersectionObserver" in window) new IntersectionObserver(es => { vis = es[0].isIntersecting; }).observe(cv);
    const render = () => { ctx.clearRect(0, 0, W, H); ctx.fillStyle = css("--surface"); ctx.fillRect(0, 0, W, H); draw(ctx, t, W, H); };
    render();
    function frame(now) {
      if (!paused && vis) { if (prev != null) t = (t + (now - prev) / 1000) % dur; prev = now; render(); } else prev = null;
      requestAnimationFrame(frame);
    }
    requestAnimationFrame(frame);
    matchMedia("(prefers-color-scheme: dark)").addEventListener("change", render);
  }

  // Step 1: autocorrelation of BOTH signals -> Bartlett factors -> slowest signal sets the window
  anim("a1", 11, (ctx, t, W, H) => {
    const ink = css("--ink"), ink2 = css("--ink-2"), ink3 = css("--ink-3"), acc = css("--accent"), or = css("--atds"), grid = css("--grid");
    const SH = 600, x0 = 40, x1 = 640, sx = i => x0 + i / (SH - 1) * (x1 - x0);
    label(ctx, "Two signals of one system (first 600 samples)", x0, 16, ink2, font(12, 600));
    trace(ctx, X.slice(0, SH), x0, x1, 48, 8, ink, 1.2); label(ctx, "x", x1 + 8, 52, ink, font(12, 600));
    trace(ctx, Y.slice(0, SH), x0, x1, 84, 8, acc, 1.2); label(ctx, "y", x1 + 8, 88, acc, font(12, 600));
    const KS = Math.max(AX.K, AY.K) + 8, p1 = clamp(t / 4.5, 0, 1), kNow = Math.floor(p1 * KS);
    const bx0 = 40, bx1 = 380, bw = (bx1 - bx0) / KS, bh = 34;
    label(ctx, "Step 0 · autocorrelation r(k) of each signal", bx0, 126, ink2, font(12, 600));
    [[AX, 170, ink, "x"], [AY, 222, acc, "y"]].forEach(([A, y0, col, nm]) => {
      ctx.strokeStyle = ink3; ctx.beginPath(); ctx.moveTo(bx0, y0); ctx.lineTo(bx1, y0); ctx.stroke();
      label(ctx, nm, bx0 - 14, y0 - 4, col, font(12, 600));
      for (let k = 0; k <= Math.min(kNow, KS - 1); k++) { const h = A.r[k] * bh; ctx.fillStyle = col; ctx.globalAlpha = k < A.K ? 0.9 : 0.25; ctx.fillRect(bx0 + k * bw, h >= 0 ? y0 - h : y0, Math.max(1, bw - 0.6), Math.abs(h)); }
      ctx.globalAlpha = 1;
      if (kNow >= A.K) { const xk = bx0 + A.K * bw; ctx.setLineDash([3, 3]); ctx.strokeStyle = ink3; ctx.beginPath(); ctx.moveTo(xk, y0 - bh); ctx.lineTo(xk, y0 + 4); ctx.stroke(); ctx.setLineDash([]); label(ctx, `K = ${A.K}`, xk + 4, y0 - bh + 9, ink3, font(10)); }
    });
    label(ctx, "lag k", bx1, 238, ink3, font(10), "right");
    const bNow = A => A.Bc[Math.max(0, Math.min(kNow, A.K - 1))];
    label(ctx, "Step 1 · B = 1 + 2 Σ r(k)²", 410, 136, ink2, mono(12));
    label(ctx, `B of x = ${bNow(AX).toFixed(1)}`, 410, 160, ink, mono(13));
    label(ctx, `B of y = ${bNow(AY).toFixed(1)}`, 410, 180, acc, mono(13));
    if (t > 4.7) label(ctx, `slowest: x, use B = ${B.toFixed(1)}`, 410, 202, ink, mono(12, 600));
    const p3 = clamp((t - 5.4) / 4, 0, 1);
    if (p3 > 0) {
      const n = Math.round(p3 * 30), a0 = 30, Lc = n * B;
      ctx.fillStyle = or; ctx.globalAlpha = 0.1; ctx.fillRect(sx(a0), 32, sx(a0 + Lc) - sx(a0), 68); ctx.globalAlpha = 1;
      ctx.strokeStyle = or; ctx.lineWidth = 2; ctx.strokeRect(sx(a0), 32, sx(a0 + Lc) - sx(a0), 68); ctx.lineWidth = 1;
      ctx.fillStyle = or; for (let j = 0; j < n; j++) { ctx.beginPath(); ctx.arc(sx(a0 + (j + 0.5) * B), 106, 2.3, 0, 7); ctx.fill(); }
      label(ctx, `${n} independent samples (dots)`, 410, 222, or, mono(12));
      label(ctx, n < 30 ? `${n} × B = ${n} × ${B.toFixed(1)} = ${Math.round(Lc)}` : `L = 30 × ${B.toFixed(1)} = ${L} samples`, 410, 242, or, mono(12, 600));
    }
  }, 250);

  // Step 2: lag sweep inside one window
  const s0 = 800, xw = X.slice(s0, s0 + L), yw = Y.slice(s0, s0 + L);
  const z = a => { const m = a.reduce((p, q) => p + q) / a.length, s = Math.sqrt(a.reduce((p, q) => p + (q - m) ** 2, 0) / a.length); return a.map(v => (v - m) / s); };
  const xz = z(xw), yz = z(yw), C = {};
  for (let k = -LAG; k <= LAG; k++) { let s = 0; for (let i = 0; i < L; i++) s += xz[((i + k) % L + L) % L] * yz[i]; C[k] = s / L; }
  let kPk = 0; for (let k = -LAG; k <= LAG; k++) if (Math.abs(C[k]) > Math.abs(C[kPk])) kPk = k;
  anim("a2", 9, (ctx, t, W, H) => {
    const ink = css("--ink"), ink2 = css("--ink-2"), ink3 = css("--ink-3"), acc = css("--accent"), or = css("--atds"), grid = css("--grid"), tint = css("--tint");
    const x0 = 40, x1 = 660, p = clamp(t / 6.5, 0, 1), kNow = Math.round(-LAG + p * 2 * LAG);
    const kShow = t > 6.5 ? kPk : kNow;
    label(ctx, `One window (L = ${L}): x fixed, y shifted by k`, x0, 16, ink2, font(12, 600));
    const sx = i => x0 + i / (L - 1) * (x1 - x0);
    trace(ctx, xz, x0, x1, 48, 9, ink, 1.3);
    ctx.beginPath(); ctx.strokeStyle = acc; ctx.lineWidth = 1.3; let prevX = null;
    yz.forEach((v, i) => { const j = ((i + kShow) % L + L) % L, X_ = sx(j), Y_ = 98 - v * 9; if (prevX == null || X_ < prevX) ctx.moveTo(X_, Y_); else ctx.lineTo(X_, Y_); prevX = X_; });
    ctx.stroke(); ctx.lineWidth = 1;
    label(ctx, "x", x1 + 6, 52, ink2, font(12)); label(ctx, "y", x1 + 6, 102, acc, font(12));
    const cy0 = 182, ch = 42, cx = k => x0 + (k + LAG) / (2 * LAG) * (x1 - x0);
    ctx.fillStyle = tint; ctx.fillRect(cx(-30), cy0 - ch - 4, cx(30) - cx(-30), 2 * ch + 8);
    label(ctx, "fixed ±30", cx(0), cy0 + ch + 16, ink3, font(11), "center");
    ctx.strokeStyle = grid; ctx.beginPath(); ctx.moveTo(x0, cy0); ctx.lineTo(x1, cy0); ctx.stroke();
    ctx.beginPath(); ctx.strokeStyle = or; ctx.lineWidth = 2;
    for (let k = -LAG; k <= (t > 6.5 ? LAG : kNow); k++) { const X_ = cx(k), Y_ = cy0 - C[k] * ch / 0.7; k === -LAG ? ctx.moveTo(X_, Y_) : ctx.lineTo(X_, Y_); }
    ctx.stroke(); ctx.lineWidth = 1;
    label(ctx, `−${LAG}`, x0, cy0 + ch + 16, ink3, mono(11)); label(ctx, `+${LAG}`, x1, cy0 + ch + 16, ink3, mono(11), "right");
    label(ctx, "C(k), aTDS lag range ±(L−1)/2", x0, 132, ink2, font(12, 600));
    if (t <= 6.5) { ctx.strokeStyle = ink3; ctx.beginPath(); ctx.moveTo(cx(kNow), cy0 - ch - 4); ctx.lineTo(cx(kNow), cy0 + ch + 4); ctx.stroke(); label(ctx, `k = ${kNow}`, x1, 132, ink, mono(13), "right"); }
    else { const X_ = cx(kPk), Y_ = cy0 - C[kPk] * ch / 0.7; ctx.fillStyle = or; ctx.beginPath(); ctx.arc(X_, Y_, 5, 0, 7); ctx.fill(); label(ctx, `τ₀ = ${kPk}: largest |C|`, x1, 132, or, mono(13), "right"); }
  });

  // Step 3: rotated copies -> tau0 per window -> stable windows at tol -> copy scores -> C(tol)
  function tau0(s1, s2) {
    const z_ = a_ => { const m = a_.reduce((p, q) => p + q) / a_.length, s = Math.sqrt(a_.reduce((p, q) => p + (q - m) ** 2, 0) / a_.length) || 1; return a_.map(v => (v - m) / s); };
    const out = [];
    for (let s = 0; s + L <= N; s += STEP) {
      const u = z_(s1.slice(s, s + L)), v = z_(s2.slice(s, s + L)); let best = 0, bk = 0;
      for (let k = -LAG; k <= LAG; k++) { let c = 0; for (let q = 0; q < L; q++) c += u[((q + k) % L + L) % L] * v[q]; if (Math.abs(c) > best) { best = Math.abs(c); bk = k; } }
      out.push(bk);
    }
    return out;
  }
  function stable(tau, tol) {
    const lab = tau.map(() => false);
    for (let i = 0; i + 5 <= tau.length; i++) { const w = tau.slice(i, i + 5);
      for (const c of new Set(w)) { const near = w.map(v => Math.abs(v - c) <= tol); if (near.filter(Boolean).length >= 4) near.forEach((ok, q) => { if (ok) lab[i + q] = true; }); } }
    return lab;
  }
  const NC = 8, aSh = Math.max(L, Math.floor(N / 10)), shifts = [], nullT = [];
  for (let j = 0; j < NC; j++) { const h = aSh + Math.floor(rnd() * (N - 2 * aSh)); shifts.push(h); const yr = Y.map((_, i) => Y[((i - h) % N + N) % N]); nullT.push(tau0(X, yr)); }
  const TM = Math.max(1, Math.floor(LAG / 4)), SC = [], CT = [], UB = [];
  for (let tol = 0; tol <= TM; tol++) { const sc = nullT.map(tt => { const lb = stable(tt, tol); return 100 * lb.filter(Boolean).length / lb.length; }); SC.push(sc); const m_ = sc.reduce((p, q) => p + q) / NC; CT.push(m_); UB.push(m_ + 2 * Math.sqrt(sc.reduce((p, q) => p + (q - m_) ** 2, 0) / (NC - 1)) / Math.sqrt(NC)); }
  let pick = 0, stop = TM; for (let tol = 0; tol <= TM; tol++) { if (UB[tol] <= 5) pick = tol; else { stop = tol; break; } }
  function vlabel(ctx, txt, x, y, col) { ctx.save(); ctx.translate(x, y); ctx.rotate(-Math.PI / 2); label(ctx, txt, 0, 0, col, font(11), "center"); ctx.restore(); }
  function axisX(ctx, x0, x1, y, ticks, fx, col, title, tcol) {
    ctx.strokeStyle = col; ctx.beginPath(); ctx.moveTo(x0, y); ctx.lineTo(x1, y); ctx.stroke();
    ticks.forEach(v => { const X_ = fx(v); ctx.beginPath(); ctx.moveTo(X_, y); ctx.lineTo(X_, y + 4); ctx.stroke(); label(ctx, String(v), X_, y + 15, tcol, mono(10), "center"); });
    if (title) label(ctx, title, (x0 + x1) / 2, y + 28, tcol, font(11), "center");
  }
  function axisY(ctx, x, y0, y1, ticks, fy, col, title, tcol, fmt_) {
    ctx.strokeStyle = col; ctx.beginPath(); ctx.moveTo(x, y0); ctx.lineTo(x, y1); ctx.stroke();
    ticks.forEach(v => { const Y_ = fy(v); ctx.beginPath(); ctx.moveTo(x - 4, Y_); ctx.lineTo(x, Y_); ctx.stroke(); label(ctx, fmt_ ? fmt_(v) : String(v), x - 6, Y_ + 3, tcol, mono(10), "right"); });
    if (title) vlabel(ctx, title, x - 34, (y0 + y1) / 2, tcol);
  }
  anim("a3", 13, (ctx, t, W, H) => {
    const ink = css("--ink"), ink2 = css("--ink-2"), ink3 = css("--ink-3"), acc = css("--accent"), or = css("--atds"), bl = css("--fixed"), grid = css("--grid");
    // Panel A: the two signals, y rotated
    const x0 = 70, x1 = 650, SH = 600, h = Math.round(ease(clamp(t / 2.5, 0, 1)) * shifts[0]);
    label(ctx, h ? `A · Rotate y by ${h} samples: each signal unchanged, coupling destroyed` : "A · x and y, coupled", 20, 14, ink2, font(12, 600));
    trace(ctx, X.slice(0, SH), x0, x1, 38, 7, ink, 1.1);
    const yr = []; for (let i = 0; i < SH; i++) yr.push(Y[((i - h) % N + N) % N]);
    trace(ctx, yr, x0, x1, 70, 7, acc, 1.1);
    label(ctx, "x", x0 - 10, 42, ink, font(12, 600), "right"); label(ctx, "y", x0 - 10, 74, acc, font(12, 600), "right");
    axisX(ctx, x0, x1, 90, [0, 200, 400, 600], v => x0 + v / (SH - 1) * (x1 - x0), ink3, "time (samples)", ink3);
    if (t < 2.6) return;
    const steps = stop + 1, p = clamp((t - 2.8) / 7, 0, 1), tol = t > 10 ? pick : Math.min(stop, Math.floor(p * steps));
    // Panel B: tau0 of copy 1, stable windows filled
    const tt = nullT[0], lb = stable(tt, tol), mx0 = 70, mx1 = 410, my0 = 150, my1 = 232, nW = tt.length;
    label(ctx, `B · Copy 1: τ₀ in each window (filled = stable at ±${tol})`, 20, 138, ink2, font(12, 600));
    const fyB = v => my1 - (v + LAG) / (2 * LAG) * (my1 - my0), fxB = i => mx0 + i / (nW - 1) * (mx1 - mx0);
    ctx.strokeStyle = grid; ctx.beginPath(); ctx.moveTo(mx0, fyB(0)); ctx.lineTo(mx1, fyB(0)); ctx.stroke();
    axisY(ctx, mx0, my0, my1, [-LAG, 0, LAG], fyB, ink3, "τ₀ (samples)", ink3);
    axisX(ctx, mx0, mx1, my1 + 6, [1, Math.round(nW / 2), nW], v => fxB(v - 1), ink3, "window number", ink3);
    tt.forEach((v, i) => { ctx.beginPath(); ctx.arc(fxB(i), fyB(v), 2.6, 0, 7); if (lb[i]) { ctx.fillStyle = bl; ctx.fill(); } else { ctx.strokeStyle = ink3; ctx.stroke(); } });
    label(ctx, `score of copy 1 = ${lb.filter(Boolean).length} / ${nW} windows = ${SC[tol][0].toFixed(1)} %`, mx0, 286, ink, mono(11));
    // Panel C: scores of all copies and their mean
    const rx0 = 500, rx1 = 660, ry1 = 232, ry0 = 150, smax = 20, bw = (rx1 - rx0) / NC;
    const fyC = v => ry1 - Math.min(v, smax) / smax * (ry1 - ry0);
    label(ctx, `C · Score of each copy`, 450, 138, ink2, font(12, 600));
    axisY(ctx, rx0, ry0, ry1, [0, 5, 10, 15, 20], fyC, ink3, "score (%)", ink3);
    SC[tol].forEach((v, j) => { const Y_ = fyC(v); ctx.fillStyle = j === 0 ? bl : ink3; ctx.globalAlpha = j === 0 ? 1 : 0.55; ctx.fillRect(rx0 + j * bw + 3, Y_, bw - 6, ry1 - Y_); });
    ctx.globalAlpha = 1;
    axisX(ctx, rx0, rx1, ry1 + 6, [], null, ink3, "copy number 1 … " + NC, ink3);
    for (let j = 0; j < NC; j++) label(ctx, String(j + 1), rx0 + (j + 0.5) * bw, ry1 + 21, ink3, mono(10), "center");
    const yM = fyC(CT[tol]); ctx.strokeStyle = or; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(rx0, yM); ctx.lineTo(rx1, yM); ctx.stroke(); ctx.lineWidth = 1;
    label(ctx, `mean (orange) = C(±${tol}) = ${CT[tol].toFixed(1)} %`, 450, 286, or, mono(11, 600));
    // Panel D: C(tol) curve
    const gx0 = 70, gx1 = 400, gy1 = 392, gy0 = 322, tmx = Math.min(TM, stop + 4), ymax = 10;
    const fxD = v => gx0 + v / tmx * (gx1 - gx0), fyD = v => gy1 - Math.min(v, ymax) / ymax * (gy1 - gy0);
    label(ctx, "D · C(tol) (dot) and its upper bound C + 2 SE (bar top)", 20, 310, ink2, font(12, 600));
    axisY(ctx, gx0, gy0, gy1, [0, 5, 10], fyD, ink3, "C(tol) (%)", ink3);
    const tks = []; for (let q = 0; q <= tmx; q += Math.max(1, Math.round(tmx / 5))) tks.push(q);
    axisX(ctx, gx0, gx1, gy1, tks, fxD, ink3, "tolerance (± samples)", ink3);
    ctx.setLineDash([4, 4]); ctx.strokeStyle = ink3; ctx.beginPath(); ctx.moveTo(gx0, fyD(5)); ctx.lineTo(gx1, fyD(5)); ctx.stroke(); ctx.setLineDash([]);
    label(ctx, "5 % budget", gx1 + 4, fyD(5) + 4, ink3, mono(10));
    const upto = t > 10 ? stop : tol;
    for (let q = 0; q <= upto; q++) { ctx.strokeStyle = UB[q] > 5 ? ink : bl; ctx.beginPath(); ctx.moveTo(fxD(q), fyD(CT[q])); ctx.lineTo(fxD(q), fyD(UB[q])); ctx.stroke(); ctx.beginPath(); ctx.moveTo(fxD(q) - 3, fyD(UB[q])); ctx.lineTo(fxD(q) + 3, fyD(UB[q])); ctx.stroke(); ctx.beginPath(); ctx.arc(fxD(q), fyD(CT[q]), 3.5, 0, 7); ctx.fillStyle = UB[q] > 5 ? ink : bl; ctx.fill(); }
    const over = UB[t > 10 ? stop : tol] > 5;
    label(ctx, `tol ±${tol}: C ${CT[tol].toFixed(1)}, C+2SE ${UB[tol].toFixed(1)} %`, 470, 340, ink, mono(11));
    label(ctx, UB[tol] <= 5 ? "C+2SE ≤ 5 %: keep, try wider" : "C+2SE > 5 %: stop", 470, 360, ink2, mono(11));
    if (t > 10 || (over && tol === stop)) label(ctx, `chosen: ±${pick} (last one ≤ 5 %)`, 470, 384, or, mono(12, 600));
  }, 430);
})();
