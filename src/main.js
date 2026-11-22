/**
 * Application Controller: wires the 3D WebGL Viewport, Solver Web Worker,
 * Telemetry Charts, Solver Benchmark Lab, and KaTeX Technical Monograph.
 */

import * as THREE from 'three';
import { GEOMETRY_PRESETS, MESH_TIERS } from './physics/batteryParams.js';
import { BatteryViewport3D, sampleColormap } from './viz/viewport3d.js';
import {
  drawDischargeChart,
  drawProfileChart,
  drawConvergenceRaceChart,
  drawConditionScalingChart,
  drawFig29Diagnostic,
} from './viz/charts.js';
import { renderMonographHTML } from './docs/monograph.js';

const worker = new Worker(new URL('./solver/solverWorker.js', import.meta.url), {
  type: 'module',
});

let viewport = null;
let latestSnapshot = null;
let latestBenchResult = null;
let isDischarging = false;
let katexRendered = false;

// DOM Elements
const statusPill = document.getElementById('solver-status-pill');
const presetSelect = document.getElementById('preset-select');
const benchPresetSelect = document.getElementById('bench-preset');
const presetDesc = document.getElementById('preset-desc');
const bcSelect = document.getElementById('bc-select');
const meshSelect = document.getElementById('mesh-select');
const clusteredCheck = document.getElementById('clustered-check');
const solverSelect = document.getElementById('solver-select');
const crateSlider = document.getElementById('crate-slider');
const crateVal = document.getElementById('crate-val');
const dtSelect = document.getElementById('dt-select');
const thermalSelect = document.getElementById('thermal-select');

const btnPlay = document.getElementById('btn-play');
const btnStep = document.getElementById('btn-step');
const btnReset = document.getElementById('btn-reset');
const btnCamReset = document.getElementById('btn-cam-reset');
const btnRunBench = document.getElementById('btn-run-bench');

const fieldSelect = document.getElementById('field-select');
const cmapSelect = document.getElementById('cmap-select');
const exagSlider = document.getElementById('exag-slider');
const exagVal = document.getElementById('exag-val');
const explodeCheck = document.getElementById('explode-check');
const vectorsCheck = document.getElementById('vectors-check');
const sliceXSlider = document.getElementById('slicex-slider');
const sliceYSlider = document.getElementById('slicey-slider');
const sliceZSlider = document.getElementById('slicez-slider');
const sliceXVal = document.getElementById('slicex-val');
const sliceYVal = document.getElementById('slicey-val');
const sliceZVal = document.getElementById('slicez-val');

const kpiVcell = document.getElementById('kpi-vcell');
const kpiTime = document.getElementById('kpi-time');
const kpiTemp = document.getElementById('kpi-temp');
const kpiIters = document.getElementById('kpi-iters');
const kpiCharge = document.getElementById('kpi-charge');

const hudTitle = document.getElementById('hud-field-title');
const hudBar = document.getElementById('hud-colorbar');
const hudMin = document.getElementById('hud-min');
const hudMid = document.getElementById('hud-mid');
const hudMax = document.getElementById('hud-max');
const probeBox = document.getElementById('probe-box');
const profileCoord = document.getElementById('profile-coord');

const chartDischarge = document.getElementById('chart-discharge');
const chartProfile = document.getElementById('chart-profile');
const regionStatsTbody = document.querySelector('#region-stats-table tbody');

const chartRace = document.getElementById('chart-race');
const chartFig29 = document.getElementById('chart-fig29');
const chartCondContrast = document.getElementById('chart-cond-contrast');
const chartCondAspect = document.getElementById('chart-cond-aspect');

function setStatus(text, kind = 'ready') {
  statusPill.textContent = text;
  statusPill.className = `status-pill status-${kind}`;
}

function populatePresets() {
  const keys = Object.keys(GEOMETRY_PRESETS);
  for (const key of keys) {
    const p = GEOMETRY_PRESETS[key];
    const opt1 = document.createElement('option');
    opt1.value = key;
    opt1.textContent = p.label;
    presetSelect.appendChild(opt1);

    const opt2 = document.createElement('option');
    opt2.value = key;
    opt2.textContent = p.label;
    benchPresetSelect.appendChild(opt2);
  }
  presetSelect.value = 'report_fig28';
  benchPresetSelect.value = 'report_fig28';
  syncPresetDefaults('report_fig28');
}

function syncPresetDefaults(presetId) {
  const p = GEOMETRY_PRESETS[presetId];
  if (!p) return;
  presetDesc.textContent = `${p.subtitle} — ${p.description}`;
  bcSelect.value = p.bc;
  crateSlider.value = String(p.defaultCRate);
  thermalSelect.value = p.defaultThermalBC;
  updateCrateLabel();
}

function updateCrateLabel() {
  const cRate = parseFloat(crateSlider.value);
  const mAcm2 = (cRate * 1.75).toFixed(2);
  crateVal.textContent = `${cRate.toFixed(1)} C (${mAcm2} mA/cm²)`;
}

function buildConfigFromUI(presetKey = presetSelect.value, useUIOptions = true) {
  const preset = GEOMETRY_PRESETS[presetKey] || GEOMETRY_PRESETS.report_fig28;
  const tier = MESH_TIERS[meshSelect.value] || MESH_TIERS.standard;
  return {
    Ly: preset.Ly,
    Lz: preset.Lz,
    nxCounts: tier.nxCounts,
    ny: tier.ny,
    nz: tier.nz,
    bc: useUIOptions ? bcSelect.value : preset.bc,
    clustered: useUIOptions ? clusteredCheck.checked : false,
    cRate: useUIOptions ? parseFloat(crateSlider.value) : preset.defaultCRate,
    solverMode: solverSelect.value,
    thermalBC: useUIOptions ? thermalSelect.value : preset.defaultThermalBC,
    hConv: 15.0,
  };
}

function initSimulation() {
  isDischarging = false;
  btnPlay.textContent = '▶ Discharge';
  setStatus('Solving 3D initial state…', 'running');
  const config = buildConfigFromUI(presetSelect.value, true);
  worker.postMessage({
    cmd: 'init',
    payload: {
      config,
      profileFrac: { yFrac: 0.5, zFrac: 0.5 },
    },
  });
}

function updateColorbarHUD(range) {
  hudTitle.textContent = `${range.label} (${range.unit})`;
  const stops = [];
  const tmp = new THREE.Color();
  for (let i = 0; i <= 10; i++) {
    const t = i / 10;
    sampleColormap(range.cmap, t, tmp);
    const r = Math.round(tmp.r * 255);
    const g = Math.round(tmp.g * 255);
    const b = Math.round(tmp.b * 255);
    stops.push(`rgb(${r},${g},${b}) ${i * 10}%`);
  }
  hudBar.style.background = `linear-gradient(90deg, ${stops.join(', ')})`;

  const span = Math.abs(range.max - range.min);
  const fmt = (v) => {
    if (span > 0 && span < 1.5) return v.toFixed(3);
    if (Math.abs(v) >= 100) return v.toFixed(1);
    return v.toFixed(4);
  };
  hudMin.textContent = fmt(range.min);
  hudMid.textContent = fmt(0.5 * (range.min + range.max));
  hudMax.textContent = fmt(range.max);
}

function updateTelemetryUI(snap) {
  latestSnapshot = snap;
  const d = snap.diagnostics;
  if (!d) return;

  kpiVcell.textContent = `${d.Vcell.toFixed(4)} V`;
  kpiTime.textContent = `${Math.round(snap.time)} s · ${(d.soc * 100).toFixed(1)}%`;
  kpiTemp.textContent = `${d.Tmax.toFixed(2)} K (${(d.Tmax - 273.15).toFixed(1)} °C)`;
  kpiIters.textContent = `${d.newtonIters} Nwt / ${d.lastLinIters} Lin (${d.solveMs.toFixed(1)} ms)`;
  kpiCharge.textContent = d.chargeErr.toExponential(2);

  // Region stats table (Report Table 18 format)
  regionStatsTbody.innerHTML = '';
  for (const row of d.regionStats) {
    const tr = document.createElement('tr');
    const phiSStr = row.hasPhiS ? `[${row.minPhiS.toFixed(4)}, ${row.maxPhiS.toFixed(4)}]` : 'n/a';
    const phiEStr = row.hasPhiE ? `[${row.minPhiE.toFixed(4)}, ${row.maxPhiE.toFixed(4)}]` : 'n/a';
    const jStr = row.hasJLi ? `[${row.minJLi.toFixed(3)}, ${row.maxJLi.toFixed(3)}]` : 'n/a';
    tr.innerHTML = `
      <td><strong>${row.region}</strong></td>
      <td>${phiSStr}</td>
      <td>${phiEStr}</td>
      <td>${jStr}</td>
      <td>${row.maxT.toFixed(2)}</td>
    `;
    regionStatsTbody.appendChild(tr);
  }

  drawDischargeChart(chartDischarge, snap.history);
  drawProfileChart(chartProfile, snap.profile, fieldSelect.value);
}

function renderProbeInfo(info) {
  const parts = [
    `<strong>${info.regionName}</strong> @ (x=${info.xUm.toFixed(1)} µm, y=${info.yCm.toFixed(2)} cm, z=${info.zCm.toFixed(2)} cm)`,
    `T=${info.T.toFixed(2)} K`,
  ];
  if (info.phiS !== null) parts.push(`φₛ=${info.phiS.toFixed(5)} V`);
  if (info.phiE !== null) parts.push(`φₑ=${info.phiE.toFixed(5)} V`);
  if (info.jLi !== null) parts.push(`jᴸⁱ=${info.jLi.toFixed(4)} A/cm³`);
  if (info.ce !== null) parts.push(`cₑ=${info.ce.toFixed(1)} mol/m³`);
  if (info.stoich !== null) parts.push(`θ=${info.stoich.toFixed(3)}`);
  parts.push(`|iₛ|=${info.iMagS.toFixed(3)} A/cm²`);
  probeBox.innerHTML = parts.join(' &nbsp;·&nbsp; ');

  // Also update 1D through-thickness profile to match probed (y, z)
  if (latestSnapshot) {
    const { ny, nz } = latestSnapshot.meshMeta;
    const yFrac = (info.j + 0.5) / ny;
    const zFrac = (info.k + 0.5) / nz;
    profileCoord.textContent = `y=${Math.round(yFrac * 100)}%, z=${Math.round(zFrac * 100)}%`;
    worker.postMessage({ cmd: 'update_profile_slice', payload: { yFrac, zFrac } });
  }
}

function renderBenchmarkCharts(res) {
  latestBenchResult = res;
  const { race, scaling } = res;

  const iluTag =
    race.ilu.finalRel <= 1e-8
      ? `${race.ilu.iters} iters`
      : `STALLED @ ${race.ilu.iters} (${race.ilu.finalRel.toExponential(1)})`;
  document.getElementById('bench-ilu-val').textContent = iluTag;
  document.getElementById('bench-hyp-val').textContent = `${race.fluxHyperbolic.iters} iters (${race.fluxHyperbolic.ms.toFixed(1)} ms)`;
  document.getElementById('bench-mg-val').textContent = `${race.xlineMG.iters} iters (${race.xlineMG.ms.toFixed(1)} ms)`;
  document.getElementById('bench-lap-val').textContent = `${race.scalarLaplace.iters} iters (rel ${race.scalarLaplace.finalRel.toExponential(1)})`;

  drawConvergenceRaceChart(chartRace, race);
  drawFig29Diagnostic(chartFig29, race.fig29);
  drawConditionScalingChart(
    chartCondContrast,
    scaling.contrastResults,
    'Collector-to-Cathode Conductivity Contrast (σ_CC / σ_p)'
  );
  drawConditionScalingChart(
    chartCondAspect,
    scaling.aspectResults,
    'Cell Height-to-Thickness Aspect Ratio (L_y / L_x)'
  );
}

function switchTab(tabId) {
  document.querySelectorAll('.nav-tab').forEach((btn) => {
    const active = btn.dataset.tab === tabId;
    btn.classList.toggle('active', active);
    btn.setAttribute('aria-selected', String(active));
  });
  document.querySelectorAll('.tab-panel').forEach((panel) => {
    panel.classList.toggle('active', panel.id === `tab-${tabId}`);
  });

  if (tabId === 'sim' && latestSnapshot) {
    drawDischargeChart(chartDischarge, latestSnapshot.history);
    drawProfileChart(chartProfile, latestSnapshot.profile, fieldSelect.value);
  } else if (tabId === 'bench') {
    if (latestBenchResult) {
      renderBenchmarkCharts(latestBenchResult);
    } else {
      triggerBenchmarkRun();
    }
  } else if (tabId === 'docs') {
    ensureKatexRendered();
  }
}

function triggerBenchmarkRun() {
  setStatus('Running live solver benchmark…', 'running');
  const key = benchPresetSelect.value || 'report_fig28';
  const preset = GEOMETRY_PRESETS[key];
  const tier = MESH_TIERS.fast;
  worker.postMessage({
    cmd: 'run_benchmark',
    payload: {
      config: {
        Ly: preset.Ly,
        Lz: preset.Lz,
        nxCounts: tier.nxCounts,
        ny: tier.ny,
        nz: tier.nz,
        bc: preset.bc,
        cRate: preset.defaultCRate,
        solverMode: 'hyperbolic_minres',
      },
    },
  });
}

function ensureKatexRendered() {
  if (katexRendered) return;
  const root = document.getElementById('monograph-content');
  if (!root || typeof window.renderMathInElement !== 'function') return;
  window.renderMathInElement(root, {
    delimiters: [
      { left: '$$', right: '$$', display: true },
      { left: '\\[', right: '\\]', display: true },
      { left: '$', right: '$', display: false },
      { left: '\\(', right: '\\)', display: false },
    ],
    throwOnError: false,
  });
  katexRendered = true;
}

// Worker Message Handler
worker.onmessage = (e) => {
  const { type, snapshot, profile, result, reason, message } = e.data;
  if (type === 'initialized' || type === 'snapshot') {
    viewport.updateSnapshot(snapshot);
    updateTelemetryUI(snapshot);
    if (type === 'initialized') {
      setStatus(`Ready · ${snapshot.meshMeta.N.toLocaleString()} cells`, 'ready');
    } else if (isDischarging) {
      setStatus(`Discharging · t=${Math.round(snapshot.time)}s`, 'running');
    } else {
      setStatus(`Ready · t=${Math.round(snapshot.time)}s`, 'ready');
    }
  } else if (type === 'profile') {
    if (latestSnapshot) {
      latestSnapshot.profile = profile;
      drawProfileChart(chartProfile, profile, fieldSelect.value);
    }
  } else if (type === 'completed') {
    isDischarging = false;
    btnPlay.textContent = '▶ Discharge';
    setStatus(reason || 'Discharge Complete', 'done');
  } else if (type === 'paused') {
    isDischarging = false;
    btnPlay.textContent = '▶ Discharge';
    setStatus('Paused', 'ready');
  } else if (type === 'benchmark_result') {
    renderBenchmarkCharts(result);
    setStatus('Benchmark Complete', 'ready');
  } else if (type === 'error') {
    isDischarging = false;
    btnPlay.textContent = '▶ Discharge';
    setStatus(`Error: ${message}`, 'done');
    console.error('Solver Worker Error:', message);
  }
};

// Bootstrap Application
function boot() {
  populatePresets();

  const vpEl = document.getElementById('viewport-3d');
  viewport = new BatteryViewport3D(vpEl, {
    onColorbarChange: updateColorbarHUD,
    onProbe: renderProbeInfo,
  });

  // Render Monograph HTML and wire jump links
  document.getElementById('monograph-root').innerHTML = renderMonographHTML();
  document.querySelectorAll('[data-jump-tab]').forEach((btn) => {
    btn.addEventListener('click', () => {
      const preset = btn.dataset.jumpPreset;
      if (preset && GEOMETRY_PRESETS[preset]) {
        presetSelect.value = preset;
        syncPresetDefaults(preset);
        initSimulation();
      }
      switchTab(btn.dataset.jumpTab);
    });
  });

  // Try KaTeX now and on load if deferred scripts finish shortly after
  ensureKatexRendered();
  window.addEventListener('load', ensureKatexRendered);

  // Navigation Tabs
  document.querySelectorAll('.nav-tab').forEach((btn) => {
    btn.addEventListener('click', () => switchTab(btn.dataset.tab));
  });

  // Geometry / Physics Controls
  presetSelect.addEventListener('change', () => {
    syncPresetDefaults(presetSelect.value);
    initSimulation();
  });
  [bcSelect, meshSelect, clusteredCheck, thermalSelect].forEach((el) => {
    el.addEventListener('change', () => initSimulation());
  });
  crateSlider.addEventListener('input', updateCrateLabel);
  crateSlider.addEventListener('change', () => initSimulation());

  solverSelect.addEventListener('change', () => {
    setStatus('Re-solving with selected method…', 'running');
    worker.postMessage({ cmd: 'update_solver_mode', payload: { solverMode: solverSelect.value } });
  });

  // Run Controls
  btnPlay.addEventListener('click', () => {
    if (isDischarging) {
      worker.postMessage({ cmd: 'pause' });
    } else {
      isDischarging = true;
      btnPlay.textContent = '⏸ Pause';
      worker.postMessage({ cmd: 'start', payload: { dt: parseFloat(dtSelect.value) } });
    }
  });

  btnStep.addEventListener('click', () => {
    if (isDischarging) return;
    worker.postMessage({ cmd: 'step', payload: { dt: parseFloat(dtSelect.value) } });
  });

  btnReset.addEventListener('click', () => {
    isDischarging = false;
    btnPlay.textContent = '▶ Discharge';
    worker.postMessage({ cmd: 'reset' });
  });

  btnCamReset.addEventListener('click', () => viewport.resetCamera());

  // 3D Visual Controls
  fieldSelect.addEventListener('change', () => {
    viewport.setViewOptions({ field: fieldSelect.value });
    if (latestSnapshot) drawProfileChart(chartProfile, latestSnapshot.profile, fieldSelect.value);
  });
  cmapSelect.addEventListener('change', () => viewport.setViewOptions({ colormap: cmapSelect.value }));
  exagSlider.addEventListener('input', () => {
    const v = parseInt(exagSlider.value, 10);
    exagVal.textContent = `${v}×`;
    viewport.setViewOptions({ exaggeration: v });
  });
  explodeCheck.addEventListener('change', () => viewport.setViewOptions({ exploded: explodeCheck.checked }));
  vectorsCheck.addEventListener('change', () => viewport.setViewOptions({ showVectors: vectorsCheck.checked }));

  const bindSlice = (slider, label, key) => {
    slider.addEventListener('input', () => {
      const val = parseFloat(slider.value);
      label.textContent = `${Math.round(val * 100)}%`;
      viewport.setViewOptions({ [key]: val });
    });
  };
  bindSlice(sliceXSlider, sliceXVal, 'sliceX');
  bindSlice(sliceYSlider, sliceYVal, 'sliceY');
  bindSlice(sliceZSlider, sliceZVal, 'sliceZ');

  // Benchmark Controls
  btnRunBench.addEventListener('click', triggerBenchmarkRun);
  benchPresetSelect.addEventListener('change', triggerBenchmarkRun);

  // Resize redraw for Canvas2D charts
  window.addEventListener('resize', () => {
    if (latestSnapshot) {
      drawDischargeChart(chartDischarge, latestSnapshot.history);
      drawProfileChart(chartProfile, latestSnapshot.profile, fieldSelect.value);
    }
    if (latestBenchResult) {
      renderBenchmarkCharts(latestBenchResult);
    }
  });

  // Launch initial 3D state on Report Fig. 28
  initSimulation();
}

boot();
