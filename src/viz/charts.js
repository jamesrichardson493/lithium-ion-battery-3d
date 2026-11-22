/**
 * Crisp HiDPI Canvas2D Scientific Plotting Suite (Zero External Dependencies).
 *
 * Plots:
 *   1. Cell Voltage V_cell(t) & Peak/Min Temperature T(t) over discharge time
 *   2. Through-thickness 1D profile across the 5 shaded battery regions
 *   3. Live Solver Convergence Race (log10 relative residual vs iteration count)
 *   4. Condition Number Scaling vs Conductivity Contrast & Geometric Aspect Ratio
 *   5. Report Fig. 29 Diagnostic (Left-Preconditioned vs True Residual Floor)
 */

function setupCanvas(canvas) {
  const rect = canvas.getBoundingClientRect();
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const w = Math.max(260, Math.floor(rect.width || canvas.clientWidth || 420));
  const h = Math.max(160, Math.floor(rect.height || canvas.clientHeight || 220));
  if (canvas.width !== Math.floor(w * dpr) || canvas.height !== Math.floor(h * dpr)) {
    canvas.width = Math.floor(w * dpr);
    canvas.height = Math.floor(h * dpr);
  }
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  return { ctx, w, h };
}

/**
 * 1. Dual-axis Voltage (V) & Max Temperature (K) Discharge Curve.
 */
export function drawDischargeChart(canvas, history) {
  const { ctx, w, h } = setupCanvas(canvas);
  const pad = { l: 48, r: 50, t: 18, b: 30 };
  const pw = w - pad.l - pad.r;
  const ph = h - pad.t - pad.b;

  // Grid
  ctx.strokeStyle = 'rgba(148, 163, 184, 0.14)';
  ctx.lineWidth = 1;
  for (let i = 0; i <= 4; i++) {
    const y = pad.t + (i / 4) * ph;
    ctx.beginPath();
    ctx.moveTo(pad.l, y);
    ctx.lineTo(pad.l + pw, y);
    ctx.stroke();
  }

  if (!history || history.length === 0) {
    ctx.fillStyle = '#64748b';
    ctx.font = '12px Inter, sans-serif';
    ctx.fillText('Press Play or Step to generate discharge telemetry...', pad.l + 12, pad.t + ph / 2);
    return;
  }

  const tMax = Math.max(60, history[history.length - 1].t);
  const vMin = 2.75;
  const vMax = 4.25;
  let tempMin = 297.5;
  let tempMax = 300.0;
  for (const pt of history) {
    if (pt.Tmax + 0.8 > tempMax) tempMax = pt.Tmax + 0.8;
    if (pt.Tmin - 0.5 < tempMin) tempMin = pt.Tmin - 0.5;
  }
  tempMax = Math.ceil(tempMax);

  const xOf = (t) => pad.l + (t / tMax) * pw;
  const yV = (v) => pad.t + ph - ((v - vMin) / (vMax - vMin)) * ph;
  const yT = (T) => pad.t + ph - ((T - tempMin) / Math.max(0.5, tempMax - tempMin)) * ph;

  // Axis labels
  ctx.fillStyle = '#94a3b8';
  ctx.font = '10px JetBrains Mono, monospace';
  ctx.textAlign = 'right';
  for (let i = 0; i <= 4; i++) {
    const val = vMax - (i / 4) * (vMax - vMin);
    ctx.fillText(val.toFixed(2), pad.l - 6, pad.t + (i / 4) * ph + 3);
  }
  ctx.textAlign = 'left';
  ctx.fillStyle = '#fb923c';
  for (let i = 0; i <= 4; i++) {
    const val = tempMax - (i / 4) * (tempMax - tempMin);
    ctx.fillText(val.toFixed(1), pad.l + pw + 6, pad.t + (i / 4) * ph + 3);
  }
  ctx.textAlign = 'center';
  ctx.fillStyle = '#94a3b8';
  for (let i = 0; i <= 4; i++) {
    const tVal = Math.round((i / 4) * tMax);
    ctx.fillText(`${tVal}s`, pad.l + (i / 4) * pw, h - 8);
  }

  // Cutoff line 2.80 V
  ctx.strokeStyle = 'rgba(244, 63, 94, 0.35)';
  ctx.setLineDash([4, 4]);
  const yCut = yV(2.8);
  ctx.beginPath();
  ctx.moveTo(pad.l, yCut);
  ctx.lineTo(pad.l + pw, yCut);
  ctx.stroke();
  ctx.setLineDash([]);

  // Temperature curve
  ctx.strokeStyle = '#fb923c';
  ctx.lineWidth = 2;
  ctx.beginPath();
  history.forEach((pt, idx) => {
    const x = xOf(pt.t);
    const y = yT(pt.Tmax);
    if (idx === 0) ctx.moveTo(x, y);
    else ctx.lineTo(x, y);
  });
  ctx.stroke();

  // Voltage curve
  ctx.strokeStyle = '#38bdf8';
  ctx.lineWidth = 2.4;
  ctx.beginPath();
  history.forEach((pt, idx) => {
    const x = xOf(pt.t);
    const y = yV(pt.V);
    if (idx === 0) ctx.moveTo(x, y);
    else ctx.lineTo(x, y);
  });
  ctx.stroke();

  // Points
  const last = history[history.length - 1];
  ctx.fillStyle = '#38bdf8';
  ctx.beginPath();
  ctx.arc(xOf(last.t), yV(last.V), 3.5, 0, Math.PI * 2);
  ctx.fill();

  ctx.fillStyle = '#fb923c';
  ctx.beginPath();
  ctx.arc(xOf(last.t), yT(last.Tmax), 3.5, 0, Math.PI * 2);
  ctx.fill();

  // Mini legend
  ctx.textAlign = 'left';
  ctx.font = '10px Inter, sans-serif';
  ctx.fillStyle = '#38bdf8';
  ctx.fillText('● V_cell (V)', pad.l + 6, pad.t + 11);
  ctx.fillStyle = '#fb923c';
  ctx.fillText('● T_max (K)', pad.l + 84, pad.t + 11);
}

/**
 * 2. Through-Thickness 1D Profile across the 5 shaded battery regions.
 */
export function drawProfileChart(canvas, profile, activeField) {
  const { ctx, w, h } = setupCanvas(canvas);
  const pad = { l: 52, r: 16, t: 18, b: 30 };
  const pw = w - pad.l - pad.r;
  const ph = h - pad.t - pad.b;

  if (!profile) return;
  const { xUm, reg } = profile;
  const xMax = 346.0;
  const xOf = (x) => pad.l + (x / xMax) * pw;

  // Shaded 5-region bands (matching thickness 10, 100, 52, 174, 10 um)
  const bands = [
    { x0: 0, x1: 10, color: 'rgba(245, 158, 11, 0.14)', label: 'Cu' },
    { x0: 10, x1: 110, color: 'rgba(56, 189, 248, 0.09)', label: 'Neg' },
    { x0: 110, x1: 162, color: 'rgba(148, 163, 184, 0.10)', label: 'Sep' },
    { x0: 162, x1: 336, color: 'rgba(244, 63, 94, 0.09)', label: 'Pos' },
    { x0: 336, x1: 346, color: 'rgba(226, 232, 240, 0.14)', label: 'Al' },
  ];

  for (const b of bands) {
    const x0 = xOf(b.x0);
    const x1 = xOf(b.x1);
    ctx.fillStyle = b.color;
    ctx.fillRect(x0, pad.t, x1 - x0, ph);
    ctx.fillStyle = 'rgba(148, 163, 184, 0.65)';
    ctx.font = '9px JetBrains Mono, monospace';
    ctx.textAlign = 'center';
    if (x1 - x0 > 14) {
      ctx.fillText(b.label, 0.5 * (x0 + x1), pad.t + 10);
    }
  }

  // Select array
  let vals = profile.phiE;
  let unit = 'V';
  let color = '#34d399';
  if (activeField === 'phiS') {
    vals = profile.phiS;
    unit = 'V';
    color = '#38bdf8';
  } else if (activeField === 'jLi') {
    vals = profile.jLi;
    unit = 'A/cm³';
    color = '#f43f5e';
  } else if (activeField === 'ce') {
    vals = profile.ce;
    unit = 'mol/m³';
    color = '#a78bfa';
  } else if (activeField === 'T') {
    vals = profile.T;
    unit = 'K';
    color = '#fb923c';
  } else if (activeField === 'stoich') {
    vals = profile.theta;
    unit = 'θ';
    color = '#facc15';
  }

  let minV = Infinity;
  let maxV = -Infinity;
  for (let i = 0; i < vals.length; i++) {
    if (!Number.isNaN(vals[i])) {
      if (vals[i] < minV) minV = vals[i];
      if (vals[i] > maxV) maxV = vals[i];
    }
  }
  if (!Number.isFinite(minV)) {
    minV = 0;
    maxV = 1;
  }
  if (Math.abs(maxV - minV) < 1e-6) {
    minV -= 0.01;
    maxV += 0.01;
  }
  const margin = 0.1 * (maxV - minV);
  minV -= margin;
  maxV += margin;

  const yOf = (v) => pad.t + ph - ((v - minV) / (maxV - minV)) * ph;

  // Grid & axes
  const vSpan = maxV - minV;
  const fmtTick = (val) => {
    if (vSpan < 0.02) return val.toFixed(4);
    if (vSpan < 2.0) return val.toFixed(2);
    if (Math.abs(val) >= 100) return val.toFixed(0);
    return val.toFixed(2);
  };
  ctx.strokeStyle = 'rgba(148, 163, 184, 0.15)';
  ctx.lineWidth = 1;
  for (let i = 0; i <= 4; i++) {
    const y = pad.t + (i / 4) * ph;
    ctx.beginPath();
    ctx.moveTo(pad.l, y);
    ctx.lineTo(pad.l + pw, y);
    ctx.stroke();
    ctx.fillStyle = '#94a3b8';
    ctx.font = '10px JetBrains Mono, monospace';
    ctx.textAlign = 'right';
    const val = maxV - (i / 4) * vSpan;
    ctx.fillText(fmtTick(val), pad.l - 5, y + 3);
  }

  ctx.textAlign = 'center';
  ctx.fillStyle = '#94a3b8';
  let lastTickX = -100;
  [0, 110, 162, 346].forEach((xv) => {
    const tx = xOf(xv);
    if (tx - lastTickX >= 42) {
      ctx.fillText(`${xv} µm`, tx, h - 8);
      lastTickX = tx;
    }
  });

  // Zero reference line for jLi
  if (minV < 0 && maxV > 0) {
    ctx.strokeStyle = 'rgba(226, 232, 240, 0.25)';
    ctx.setLineDash([3, 3]);
    ctx.beginPath();
    ctx.moveTo(pad.l, yOf(0));
    ctx.lineTo(pad.l + pw, yOf(0));
    ctx.stroke();
    ctx.setLineDash([]);
  }

  // Plot contiguous segments per domain (break across NaNs!)
  ctx.strokeStyle = color;
  ctx.lineWidth = 2.3;
  let drawing = false;
  ctx.beginPath();
  for (let i = 0; i < vals.length; i++) {
    if (Number.isNaN(vals[i])) {
      drawing = false;
      continue;
    }
    const x = xOf(xUm[i]);
    const y = yOf(vals[i]);
    if (!drawing) {
      ctx.moveTo(x, y);
      drawing = true;
    } else {
      ctx.lineTo(x, y);
    }
  }
  ctx.stroke();

  // Dots
  ctx.fillStyle = color;
  for (let i = 0; i < vals.length; i++) {
    if (Number.isNaN(vals[i])) continue;
    ctx.beginPath();
    ctx.arc(xOf(xUm[i]), yOf(vals[i]), 2.5, 0, Math.PI * 2);
    ctx.fill();
  }

  ctx.textAlign = 'right';
  ctx.fillStyle = '#cbd5e1';
  ctx.font = '10px Inter, sans-serif';
  ctx.fillText(`x-Profile (${unit})`, pad.l + pw - 4, pad.t + 10);
  void reg;
}

/**
 * 3. Solver Convergence Race Plot: log10(relative residual) vs iterations.
 */
export function drawConvergenceRaceChart(canvas, race) {
  const { ctx, w, h } = setupCanvas(canvas);
  const pad = { l: 54, r: 18, t: 24, b: 36 };
  const pw = w - pad.l - pad.r;
  const ph = h - pad.t - pad.b;

  if (!race) return;

  const series = [
    { name: 'Primal + Jacobi GMRES(60)', data: race.jacobi.history, color: '#94a3b8', dash: [4, 3] },
    { name: 'Primal + Point-ILU(0)/SSOR GMRES(60)', data: race.ilu.history, color: '#fb7185', dash: [] },
    { name: 'Fixed Sec. 4: Flux-Hyperbolic Block-MINRES', data: race.fluxHyperbolic.history, color: '#38bdf8', dash: [] },
    { name: 'Fixed Sec. 2.2: Coupled X-Line MG-PCG', data: race.xlineMG.history, color: '#34d399', dash: [] },
  ];

  const maxIter = Math.max(50, ...series.map((s) => s.data.length));
  const yMin = -11;
  const yMax = 0.5;
  const xOf = (it) => pad.l + (it / maxIter) * pw;
  const yOf = (logV) => pad.t + ph - ((Math.max(yMin, Math.min(yMax, logV)) - yMin) / (yMax - yMin)) * ph;

  // Grid
  ctx.strokeStyle = 'rgba(148, 163, 184, 0.14)';
  ctx.lineWidth = 1;
  for (let lv = 0; lv >= -10; lv -= 2) {
    const y = yOf(lv);
    ctx.beginPath();
    ctx.moveTo(pad.l, y);
    ctx.lineTo(pad.l + pw, y);
    ctx.stroke();
    ctx.fillStyle = '#94a3b8';
    ctx.font = '10px JetBrains Mono, monospace';
    ctx.textAlign = 'right';
    ctx.fillText(`1e${lv}`, pad.l - 6, y + 3);
  }

  for (let i = 0; i <= 5; i++) {
    const itVal = Math.round((i / 5) * maxIter);
    const x = xOf(itVal);
    ctx.beginPath();
    ctx.moveTo(x, pad.t);
    ctx.lineTo(x, pad.t + ph);
    ctx.stroke();
    ctx.fillStyle = '#94a3b8';
    ctx.textAlign = 'center';
    ctx.fillText(`${itVal}`, x, h - 14);
  }

  // Target tolerance line 1e-10
  ctx.strokeStyle = 'rgba(52, 211, 153, 0.35)';
  ctx.setLineDash([5, 4]);
  ctx.beginPath();
  ctx.moveTo(pad.l, yOf(-10));
  ctx.lineTo(pad.l + pw, yOf(-10));
  ctx.stroke();
  ctx.setLineDash([]);

  // Plot series
  for (const s of series) {
    ctx.strokeStyle = s.color;
    ctx.lineWidth = 2.3;
    ctx.setLineDash(s.dash);
    ctx.beginPath();
    s.data.forEach((val, idx) => {
      const logV = Math.log10(Math.max(1e-14, val));
      const x = xOf(idx);
      const y = yOf(logV);
      if (idx === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    });
    ctx.stroke();
  }
  ctx.setLineDash([]);

  // Legend bottom-left in clear space below breakthroughs
  const legW = 285;
  const legH = 66;
  const legX = pad.l + 10;
  const legY = pad.t + ph - legH - 8;
  ctx.fillStyle = 'rgba(9, 14, 26, 0.86)';
  ctx.strokeStyle = 'rgba(148, 163, 184, 0.22)';
  ctx.lineWidth = 1;
  ctx.fillRect(legX, legY, legW, legH);
  ctx.strokeRect(legX, legY, legW, legH);

  ctx.textAlign = 'left';
  ctx.font = '11px Inter, sans-serif';
  let ly = legY + 15;
  for (const s of series) {
    ctx.fillStyle = s.color;
    ctx.fillText(`━ ${s.name}`, legX + 8, ly);
    ly += 15;
  }

  ctx.fillStyle = '#94a3b8';
  ctx.font = '10px Inter, sans-serif';
  ctx.textAlign = 'center';
  ctx.fillText('Krylov / MINRES / PCG Iterations', pad.l + pw / 2, h - 2);
}

/**
 * 4. Log-Log Condition Number Scaling Chart (vs Contrast or vs Aspect Ratio).
 */
export function drawConditionScalingChart(canvas, rows, xLabel) {
  const { ctx, w, h } = setupCanvas(canvas);
  const pad = { l: 56, r: 16, t: 22, b: 36 };
  const pw = w - pad.l - pad.r;
  const ph = h - pad.t - pad.b;

  if (!rows || rows.length === 0) return;

  const logXs = rows.map((r) => Math.log10(r.x));
  const xMin = Math.min(...logXs);
  const xMax = Math.max(...logXs);
  const yMin = 0;
  const yMax = 16;

  const xOf = (lx) => pad.l + ((lx - xMin) / Math.max(1e-6, xMax - xMin)) * pw;
  const yOf = (ly) => pad.t + ph - ((Math.max(yMin, Math.min(yMax, ly)) - yMin) / (yMax - yMin)) * ph;

  ctx.strokeStyle = 'rgba(148, 163, 184, 0.14)';
  ctx.lineWidth = 1;
  for (let p = 0; p <= 16; p += 4) {
    const y = yOf(p);
    ctx.beginPath();
    ctx.moveTo(pad.l, y);
    ctx.lineTo(pad.l + pw, y);
    ctx.stroke();
    ctx.fillStyle = '#94a3b8';
    ctx.font = '10px JetBrains Mono, monospace';
    ctx.textAlign = 'right';
    ctx.fillText(`1e${p}`, pad.l - 6, y + 3);
  }

  ctx.textAlign = 'center';
  rows.forEach((r) => {
    const x = xOf(Math.log10(r.x));
    ctx.fillStyle = '#94a3b8';
    ctx.font = '10px JetBrains Mono, monospace';
    const lbl = r.x >= 1000 ? r.x.toExponential(0) : String(r.x);
    ctx.fillText(lbl, x, h - 14);
  });

  const curves = [
    { key: 'condNaiveHyp', label: 'Naive Sec. 4 (Tr = L²/σ)', color: '#f59e0b', dash: [4, 3] },
    { key: 'condElliptic', label: 'Primal Elliptic FVM', color: '#fb7185', dash: [] },
    { key: 'condScaledHyp', label: 'Fixed Sec. 4 Scaled Flux-Hyp', color: '#38bdf8', dash: [] },
    { key: 'condPrecondHyp', label: 'Block-Schur Precond (≤ φ² = 2.618)', color: '#34d399', dash: [] },
  ];

  for (const c of curves) {
    ctx.strokeStyle = c.color;
    ctx.lineWidth = 2.2;
    ctx.setLineDash(c.dash);
    ctx.beginPath();
    rows.forEach((r, idx) => {
      const x = xOf(Math.log10(r.x));
      const y = yOf(Math.log10(Math.max(1, r[c.key])));
      if (idx === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    });
    ctx.stroke();

    ctx.fillStyle = c.color;
    rows.forEach((r) => {
      const x = xOf(Math.log10(r.x));
      const y = yOf(Math.log10(Math.max(1, r[c.key])));
      ctx.beginPath();
      ctx.arc(x, y, 3, 0, Math.PI * 2);
      ctx.fill();
    });
  }
  ctx.setLineDash([]);

  const legW = 262;
  const legH = 58;
  const legX = pad.l + 10;
  const legY = pad.t + 6;
  ctx.fillStyle = 'rgba(9, 14, 26, 0.84)';
  ctx.strokeStyle = 'rgba(148, 163, 184, 0.2)';
  ctx.lineWidth = 1;
  ctx.fillRect(legX, legY, legW, legH);
  ctx.strokeRect(legX, legY, legW, legH);

  ctx.textAlign = 'left';
  ctx.font = '10px Inter, sans-serif';
  let ly = legY + 13;
  for (const c of curves) {
    ctx.fillStyle = c.color;
    ctx.fillText(`● ${c.label}`, legX + 8, ly);
    ly += 13;
  }

  ctx.fillStyle = '#94a3b8';
  ctx.textAlign = 'center';
  ctx.fillText(xLabel, pad.l + pw / 2, h - 2);
}

/**
 * 5. Report Fig. 29 Diagnostic: Left-Preconditioned Residual vs True Residual.
 */
export function drawFig29Diagnostic(canvas, fig29) {
  const { ctx, w, h } = setupCanvas(canvas);
  const pad = { l: 54, r: 18, t: 22, b: 32 };
  const pw = w - pad.l - pad.r;
  const ph = h - pad.t - pad.b;
  if (!fig29) return;

  const maxIt = Math.max(20, fig29.preHistory.length);
  const yMin = -8;
  const yMax = 0.5;
  const xOf = (it) => pad.l + (it / maxIt) * pw;
  const yOf = (lv) => pad.t + ph - ((Math.max(yMin, Math.min(yMax, lv)) - yMin) / (yMax - yMin)) * ph;

  ctx.strokeStyle = 'rgba(148, 163, 184, 0.14)';
  for (let lv = 0; lv >= -8; lv -= 2) {
    const y = yOf(lv);
    ctx.beginPath();
    ctx.moveTo(pad.l, y);
    ctx.lineTo(pad.l + pw, y);
    ctx.stroke();
    ctx.fillStyle = '#94a3b8';
    ctx.font = '10px JetBrains Mono, monospace';
    ctx.textAlign = 'right';
    ctx.fillText(`1e${lv}`, pad.l - 6, y + 3);
  }

  const drawLine = (data, color) => {
    ctx.strokeStyle = color;
    ctx.lineWidth = 2.2;
    ctx.beginPath();
    data.forEach((v, idx) => {
      const x = xOf(idx);
      const y = yOf(Math.log10(Math.max(1e-12, v)));
      if (idx === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    });
    ctx.stroke();
  };

  drawLine(fig29.preHistory, '#a78bfa');
  drawLine(fig29.trueHistory, '#fb7185');

  const legW = 268;
  const legH = 42;
  const legX = pad.l + 10;
  const legY = pad.t + ph - legH - 10;
  ctx.fillStyle = 'rgba(9, 14, 26, 0.86)';
  ctx.strokeStyle = 'rgba(148, 163, 184, 0.22)';
  ctx.lineWidth = 1;
  ctx.fillRect(legX, legY, legW, legH);
  ctx.strokeRect(legX, legY, legW, legH);

  ctx.textAlign = 'left';
  ctx.font = '11px Inter, sans-serif';
  ctx.fillStyle = '#a78bfa';
  ctx.fillText('━ Left-Preconditioned Residual ||M⁻¹r||', legX + 8, legY + 16);
  ctx.fillStyle = '#fb7185';
  ctx.fillText('━ True Residual ||b − Ax|| (Stalled!)', legX + 8, legY + 32);
}
