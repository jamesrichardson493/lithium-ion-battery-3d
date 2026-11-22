/**
 * Dedicated ES Module Web Worker running the 3D Coupled 6-Equation Solver
 * and Solver Benchmark Lab off the main UI / WebGL thread.
 */

import { Battery3DEngine } from './battery3dEngine.js';
import { computeConditionScalingBenchmarks } from './hyperbolicSolver.js';
import { runLiveSolverRace } from './naiveBaseline.js';

let engine = null;
let runTimer = null;
let isRunning = false;
let currentProfileFrac = { yFrac: 0.5, zFrac: 0.5 };

function stopLoop() {
  isRunning = false;
  if (runTimer !== null) {
    clearTimeout(runTimer);
    runTimer = null;
  }
}

function sendSnapshot(type = 'snapshot') {
  if (!engine) return;
  const snap = engine.getSnapshot(currentProfileFrac.yFrac, currentProfileFrac.zFrac);
  self.postMessage({ type, snapshot: snap });
}

self.onmessage = (e) => {
  const { cmd, payload } = e.data;

  try {
    if (cmd === 'init') {
      stopLoop();
      engine = new Battery3DEngine(payload.config);
      currentProfileFrac = payload.profileFrac || { yFrac: 0.5, zFrac: 0.5 };
      sendSnapshot('initialized');
    } else if (cmd === 'update_solver_mode') {
      if (engine) {
        engine.config.solverMode = payload.solverMode;
        engine.solvePotentials();
        sendSnapshot('snapshot');
      }
    } else if (cmd === 'update_profile_slice') {
      currentProfileFrac = payload;
      if (engine) {
        self.postMessage({
          type: 'profile',
          profile: engine.extractProfileX(currentProfileFrac.yFrac, currentProfileFrac.zFrac),
        });
      }
    } else if (cmd === 'step') {
      if (!engine || engine.terminated) return;
      engine.stepTime(payload.dt || 30.0);
      sendSnapshot('snapshot');
    } else if (cmd === 'start') {
      if (!engine) return;
      if (engine.terminated) {
        engine.reset();
      }
      isRunning = true;
      const dt = payload.dt || 30.0;
      const tick = () => {
        if (!isRunning || !engine) return;
        engine.stepTime(dt);
        sendSnapshot('snapshot');
        if (engine.terminated) {
          isRunning = false;
          self.postMessage({ type: 'completed', reason: engine.terminationReason });
          return;
        }
        runTimer = setTimeout(tick, 16);
      };
      tick();
    } else if (cmd === 'pause') {
      stopLoop();
      self.postMessage({ type: 'paused' });
    } else if (cmd === 'reset') {
      stopLoop();
      if (engine) {
        engine.reset();
        sendSnapshot('initialized');
      }
    } else if (cmd === 'run_benchmark') {
      // Run Benchmark Lab computations (fresh copy or current first-Newton benchmark)
      const benchEngine = payload?.config ? new Battery3DEngine(payload.config) : engine;
      if (!benchEngine) return;
      const snap0 = benchEngine.firstNewtonBenchmark;
      const race = runLiveSolverRace(benchEngine.mesh, snap0.kap, snap0.gbv, snap0.rhsPhi, snap0.R0);
      const scaling = computeConditionScalingBenchmarks();
      self.postMessage({
        type: 'benchmark_result',
        result: {
          race,
          scaling,
          config: benchEngine.config,
        },
      });
    }
  } catch (err) {
    stopLoop();
    self.postMessage({
      type: 'error',
      message: err?.message || String(err),
      stack: err?.stack,
    });
  }
};
