/**
 * Full 6-Equation 3D Time-Dependent Coupled Electrochemical-Thermal Battery Solver
 * (Li_x C_6 | Li_y Mn_2 O_4, 5 Regions, Doyle-Fuller-Newman / Wang-Gu-Liaw Model).
 *
 * Solves Report Eqs. (3)-(13) & Table 2 boundary conditions:
 *   - Solid species (Eq. 3): c_s & parabolic surface concentration c_se
 *   - Electrolyte species (Eqs. 4 & 9): implicit Backward Euler + ScalarXLineMG
 *   - Coupled Solid & Liquid Potentials + Butler-Volmer Reaction Current (Eqs. 5-7, 10, 12):
 *       Solved via exact analytical coupled Newton-Raphson with positive-domain
 *       gauge shift (eliminating 4.095 V cancellation in Al CC) using EITHER:
 *         (A) Fixed Sec. 4: Scaled Flux-Hyperbolic Block-Preconditioned MINRES
 *         (B) Fixed Sec. 2.2: Coupled X-Line Semi-Coarsened MG-PCG
 *   - 3D Energy / Temperature (Eqs. 8, 11, 13): implicit Backward Euler + ScalarXLineMG
 *     with exact zero-cancellation Joule heating F_f^2 R_f + reaction + entropic heat.
 */

import {
  F_CONST,
  R_GAS,
  T_AMBIENT,
  V_CUTOFF,
  C_RATE_1C,
  LAYER_SPECS,
  ELECTROLYTE_PARAMS,
  REGION_IDS,
  REGION_SHORT,
  uNeg,
  uPos,
  dUNegDT,
  dUPosDT,
  kappaBulk,
  kappaDiffusional,
  solidDiffusivity,
  reactionRate,
  electrolyteDiffusivity,
} from '../physics/batteryParams.js';
import { BatteryMesh3D } from './mesh3d.js';
import { CoupledXLineMG, ScalarXLineMG } from './multigrid.js';
import { FluxHyperbolicSolver } from './hyperbolicSolver.js';

export class Battery3DEngine {
  /**
   * @param {Object} config
   * @param {number} config.Ly - Cell height [m]
   * @param {number} config.Lz - Cell width [m]
   * @param {number[]} config.nxCounts - Cell counts per region [NegCC, Neg, Sep, Pos, PosCC]
   * @param {number} config.ny
   * @param {number} config.nz
   * @param {string} config.bc - 'BC1' | 'BC2' | 'BC3' | 'BC4'
   * @param {boolean} [config.clustered=false]
   * @param {number} [config.cRate=2.0]
   * @param {string} [config.solverMode='hyperbolic_minres']
   * @param {string} [config.thermalBC='convective']
   * @param {number} [config.hConv=15.0]
   */
  constructor(config) {
    this.config = {
      clustered: false,
      cRate: 2.0,
      solverMode: 'hyperbolic_minres',
      thermalBC: 'convective',
      hConv: 15.0,
      ...config,
    };

    this.mesh = new BatteryMesh3D({
      Ly: this.config.Ly,
      Lz: this.config.Lz,
      nxCounts: this.config.nxCounts,
      ny: this.config.ny,
      nz: this.config.nz,
      bc: this.config.bc,
      clustered: this.config.clustered,
    });

    const N = this.mesh.N;
    const Npot = this.mesh.Npot;

    // Primary 3D state fields (strict SI)
    this.phiS = new Float64Array(N);
    this.phiE = new Float64Array(N);
    this.jLi = new Float64Array(N);
    this.cs = new Float64Array(N);
    this.cse = new Float64Array(N);
    this.ce = new Float64Array(N);
    this.T = new Float64Array(N);
    this.etaField = new Float64Array(N);

    // 3D current density vectors [A/m^2] & heat sources [W/m^3]
    this.isVec = new Float64Array(3 * N);
    this.ieVec = new Float64Array(3 * N);
    this.qOhm = new Float64Array(N);
    this.qTotal = new Float64Array(N);

    // Internal solver arrays
    this.kap = new Float64Array(N);
    this.kapD = new Float64Array(N);
    this.gbv = new Float64Array(N);
    this.uPot = new Float64Array(Npot);
    this.resPot = new Float64Array(Npot);
    this.rhsPot = new Float64Array(Npot);
    this.duPot = new Float64Array(Npot);

    // Multigrid & Hyperbolic solvers
    this.coupledMG = new CoupledXLineMG(this.mesh);
    this.hypSolver = new FluxHyperbolicSolver(this.mesh, this.coupledMG);
    this.ceMG = new ScalarXLineMG(this.mesh.nxElec, this.mesh.ny, this.mesh.nz);
    this.tempMG = new ScalarXLineMG(this.mesh.nx, this.mesh.ny, this.mesh.nz);

    // Pre-allocated buffers for implicit c_e and T solves
    const Ne = this.mesh.Ne;
    this.ceX = new Float64Array(Ne);
    this.ceRhs = new Float64Array(Ne);
    this.ceDiag = new Float64Array(Ne);
    this.ceGx = new Float64Array(Ne);
    this.ceGy = new Float64Array(Math.max(1, this.mesh.ny - 1) * this.mesh.nz * this.mesh.nxElec);
    this.ceGz = new Float64Array(this.mesh.ny * Math.max(1, this.mesh.nz - 1) * this.mesh.nxElec);

    this.tX = new Float64Array(N);
    this.tRhs = new Float64Array(N);
    this.tDiag = new Float64Array(N);
    this.tGx = new Float64Array(N);
    this.tGy = new Float64Array(Math.max(1, this.mesh.ny - 1) * this.mesh.nz * this.mesh.nx);
    this.tGz = new Float64Array(this.mesh.ny * Math.max(1, this.mesh.nz - 1) * this.mesh.nx);

    // Positive-domain gauge shift reference (eliminates 4.0 V cancellation in Al CC)
    this.Vref = uPos(LAYER_SPECS[REGION_IDS.POS].theta0);

    // Reference resistance for Scaled Flux-Hyperbolic preconditioner
    this.R0 = this._computeReferenceResistance();

    // Simulation telemetry
    this.time = 0.0;
    this.stepCount = 0;
    this.terminated = false;
    this.terminationReason = '';
    this.history = [];
    this.lastDiagnostics = null;

    this.reset();
  }

  _computeReferenceResistance() {
    const { mesh } = this;
    const area = mesh.Ly * mesh.Lz;
    // Characteristic through-cell resistance [Ohm] dominated by electrolyte & electrodes
    const rEst =
      LAYER_SPECS[0].thickness / LAYER_SPECS[0].sigma +
      LAYER_SPECS[1].thickness / 0.1 +
      LAYER_SPECS[2].thickness / 0.1 +
      LAYER_SPECS[3].thickness / 0.1 +
      LAYER_SPECS[4].thickness / LAYER_SPECS[4].sigma;
    return rEst / area;
  }

  reset() {
    const { mesh, cs, cse, ce, T, phiS, phiE, jLi, etaField } = this;
    const { nx, ny, nz, reg } = mesh;
    const negSpec = LAYER_SPECS[REGION_IDS.NEG];
    const posSpec = LAYER_SPECS[REGION_IDS.POS];

    this.time = 0.0;
    this.stepCount = 0;
    this.terminated = false;
    this.terminationReason = '';
    this.history = [];
    this.Vref = uPos(posSpec.theta0);

    const Un0 = uNeg(negSpec.theta0);
    const Up0 = uPos(posSpec.theta0);

    for (let i = 0; i < nx; i++) {
      const r = reg[i];
      for (let j = 0; j < ny; j++) {
        for (let k = 0; k < nz; k++) {
          const c = mesh.cid(i, j, k);
          T[c] = T_AMBIENT;
          jLi[c] = 0.0;
          etaField[c] = 0.0;

          if (r === REGION_IDS.NEG) {
            cs[c] = negSpec.theta0 * negSpec.csMax;
            cse[c] = cs[c];
          } else if (r === REGION_IDS.POS) {
            cs[c] = posSpec.theta0 * posSpec.csMax;
            cse[c] = cs[c];
          } else {
            cs[c] = 0.0;
            cse[c] = 0.0;
          }

          if (r === REGION_IDS.NEG || r === REGION_IDS.SEP || r === REGION_IDS.POS) {
            ce[c] = ELECTROLYTE_PARAMS.ce0;
            phiE[c] = -Un0 - 0.03;
          } else {
            ce[c] = 0.0;
            phiE[c] = 0.0;
          }

          if (r === REGION_IDS.NEG_CC || r === REGION_IDS.NEG) {
            phiS[c] = 0.0;
          } else if (r === REGION_IDS.POS || r === REGION_IDS.POS_CC) {
            phiS[c] = Up0 - Un0 - 0.08;
          } else {
            phiS[c] = 0.0;
          }
        }
      }
    }

    // Solve initial consistent potential & reaction current state at t = 0
    const diag = this.solvePotentials();
    this._recordHistory(diag);
    return diag;
  }

  /**
   * Solve the coupled 3D nonlinear (phi_s, phi_e, j^Li) equations via exact Newton-Raphson
   * using positive-domain gauge shift and the selected breakthrough linear solver.
   */
  solvePotentials() {
    const {
      mesh,
      config,
      kap,
      kapD,
      gbv,
      uPot,
      resPot,
      rhsPot,
      duPot,
      coupledMG,
      hypSolver,
      Vref,
      R0,
    } = this;
    const { nx, ny, nz, dx, dy, dz, reg, sigma, epsE, sLoc, eLoc, negTabMask, posTabMask, negTabType } =
      mesh;

    const negSpec = LAYER_SPECS[REGION_IDS.NEG];
    const posSpec = LAYER_SPECS[REGION_IDS.POS];
    const areaX = dy * dz;
    const I_face = config.cRate * C_RATE_1C; // [A/m^2]
    const I_total = I_face * mesh.Ly * mesh.Lz; // [A]
    const I_pos_cell = I_total / mesh.posTabCells.length;
    const gxs0 = areaX / ((0.5 * dx[0]) / sigma[0]);

    // 1. Update effective electrolyte ionic conductivity kappa^eff and kappa_D^eff
    for (let i = 0; i < nx; i++) {
      const r = reg[i];
      if (r === REGION_IDS.NEG || r === REGION_IDS.SEP || r === REGION_IDS.POS) {
        const tort = Math.pow(epsE[i], ELECTROLYTE_PARAMS.brugg);
        for (let j = 0; j < ny; j++) {
          for (let k = 0; k < nz; k++) {
            const c = mesh.cid(i, j, k);
            const ke = kappaBulk(this.ce[c], this.T[c]) * tort;
            kap[c] = ke;
            kapD[c] = kappaDiffusional(ke, this.T[c]);
          }
        }
      } else {
        for (let j = 0; j < ny; j++) {
          for (let k = 0; k < nz; k++) {
            const c = mesh.cid(i, j, k);
            kap[c] = 0.0;
            kapD[c] = 0.0;
          }
        }
      }
    }

    // 2. Pack gauge-shifted potentials into uPot: phiS_tilde = phiS - Vref in Pos & PosCC
    for (let j = 0; j < ny; j++) {
      for (let k = 0; k < nz; k++) {
        for (let i = 0; i < nx; i++) {
          const c = mesh.cid(i, j, k);
          if (sLoc[i] >= 0) {
            const shift = reg[i] >= REGION_IDS.POS ? Vref : 0.0;
            uPot[mesh.potIdx(j, k, sLoc[i])] = this.phiS[c] - shift;
          }
          if (eLoc[i] >= 0) {
            uPot[mesh.potIdx(j, k, eLoc[i])] = this.phiE[c];
          }
        }
      }
    }

    // 3. Precompute static concentrated-solution diffusional current source:
    //    -div(kappa_D^eff grad ln c_e) added to electrolyte residual
    const diffSource = new Float64Array(mesh.Npot);
    for (let j = 0; j < ny; j++) {
      for (let k = 0; k < nz; k++) {
        for (let i = 0; i < nx; i++) {
          if (eLoc[i] < 0) continue;
          const c = mesh.cid(i, j, k);
          const lnC = Math.log(Math.max(50.0, this.ce[c]));
          const kd = kapD[c];
          const pA = mesh.potIdx(j, k, eLoc[i]);

          if (i < nx - 1 && eLoc[i + 1] >= 0) {
            const cR = mesh.cid(i + 1, j, k);
            const kdR = kapD[cR];
            const gDx = areaX / ((0.5 * dx[i]) / kd + (0.5 * dx[i + 1]) / kdR);
            const fl = gDx * (lnC - Math.log(Math.max(50.0, this.ce[cR])));
            diffSource[pA] += fl;
            diffSource[mesh.potIdx(j, k, eLoc[i + 1])] -= fl;
          }
          if (j < ny - 1) {
            const cT = mesh.cid(i, j + 1, k);
            const kdT = kapD[cT];
            const kdEff = (2.0 * kd * kdT) / (kd + kdT);
            const gDy = ((dx[i] * dz) / dy) * kdEff;
            const fl = gDy * (lnC - Math.log(Math.max(50.0, this.ce[cT])));
            diffSource[pA] += fl;
            diffSource[mesh.potIdx(j + 1, k, eLoc[i])] -= fl;
          }
          if (k < nz - 1) {
            const cF = mesh.cid(i, j, k + 1);
            const kdF = kapD[cF];
            const kdEff = (2.0 * kd * kdF) / (kd + kdF);
            const gDz = ((dx[i] * dy) / dz) * kdEff;
            const fl = gDz * (lnC - Math.log(Math.max(50.0, this.ce[cF])));
            diffSource[pA] += fl;
            diffSource[mesh.potIdx(j, k + 1, eLoc[i])] -= fl;
          }
        }
      }
    }

    if (config.solverMode === 'hyperbolic_minres') {
      hypSolver.updateFaces(kap, R0);
    }

    // 4. Newton-Raphson iteration
    const tSolveStart = performance.now();
    const newtonHistory = [];
    let totalLinIters = 0;
    let lastLinIters = 0;
    let lastLinRel = 0;
    let firstNewtonRhsSnapshot = null;
    let firstNewtonGbvSnapshot = null;

    for (let nit = 1; nit <= 18; nit++) {
      resPot.set(diffSource);
      gbv.fill(0);

      for (let j = 0; j < ny; j++) {
        for (let k = 0; k < nz; k++) {
          const tabIdx = j * nz + k;
          if (negTabType === 'left_face' && negTabMask[tabIdx]) {
            const p0 = mesh.potIdx(j, k, sLoc[0]);
            resPot[p0] += gxs0 * uPot[p0];
          } else if (negTabType === 'top_edge' && j === ny - 1 && negTabMask[tabIdx]) {
            for (let i = 0; i < mesh.nxCounts[0]; i++) {
              const p0 = mesh.potIdx(j, k, sLoc[i]);
              resPot[p0] += sigma[i] * ((dx[i] * dz) / (0.5 * dy)) * uPot[p0];
            }
          }

          if (posTabMask[tabIdx]) {
            const pLast = mesh.potIdx(j, k, sLoc[nx - 1]);
            resPot[pLast] += I_pos_cell;
          }

          for (let i = 0; i < nx - 1; i++) {
            if (sLoc[i] >= 0 && sLoc[i + 1] >= 0) {
              const gs = areaX / ((0.5 * dx[i]) / sigma[i] + (0.5 * dx[i + 1]) / sigma[i + 1]);
              const a = mesh.potIdx(j, k, sLoc[i]);
              const b = mesh.potIdx(j, k, sLoc[i + 1]);
              const fl = gs * (uPot[a] - uPot[b]);
              resPot[a] += fl;
              resPot[b] -= fl;
            }
            if (eLoc[i] >= 0 && eLoc[i + 1] >= 0) {
              const kL = kap[mesh.cid(i, j, k)];
              const kR = kap[mesh.cid(i + 1, j, k)];
              const ge = areaX / ((0.5 * dx[i]) / kL + (0.5 * dx[i + 1]) / kR);
              const a = mesh.potIdx(j, k, eLoc[i]);
              const b = mesh.potIdx(j, k, eLoc[i + 1]);
              const fl = ge * (uPot[a] - uPot[b]);
              resPot[a] += fl;
              resPot[b] -= fl;
            }
          }

          for (let i = 0; i < nx; i++) {
            const c = mesh.cid(i, j, k);
            const vol = mesh.cellVol[i];

            if (sLoc[i] >= 0) {
              const a = mesh.potIdx(j, k, sLoc[i]);
              if (j < ny - 1) {
                const g = (sigma[i] * dx[i] * dz) / dy;
                const b = mesh.potIdx(j + 1, k, sLoc[i]);
                const fl = g * (uPot[a] - uPot[b]);
                resPot[a] += fl;
                resPot[b] -= fl;
              }
              if (k < nz - 1) {
                const g = (sigma[i] * dx[i] * dy) / dz;
                const b = mesh.potIdx(j, k + 1, sLoc[i]);
                const fl = g * (uPot[a] - uPot[b]);
                resPot[a] += fl;
                resPot[b] -= fl;
              }
            }

            if (eLoc[i] >= 0) {
              const a = mesh.potIdx(j, k, eLoc[i]);
              if (j < ny - 1) {
                const kL = kap[c];
                const kT = kap[mesh.cid(i, j + 1, k)];
                const g = ((dx[i] * dz) / dy) * ((2.0 * kL * kT) / (kL + kT));
                const b = mesh.potIdx(j + 1, k, eLoc[i]);
                const fl = g * (uPot[a] - uPot[b]);
                resPot[a] += fl;
                resPot[b] -= fl;
              }
              if (k < nz - 1) {
                const kL = kap[c];
                const kF = kap[mesh.cid(i, j, k + 1)];
                const g = ((dx[i] * dy) / dz) * ((2.0 * kL * kF) / (kL + kF));
                const b = mesh.potIdx(j, k + 1, eLoc[i]);
                const fl = g * (uPot[a] - uPot[b]);
                resPot[a] += fl;
                resPot[b] -= fl;
              }
            }

            if (reg[i] === REGION_IDS.NEG || reg[i] === REGION_IDS.POS) {
              const psTilde = uPot[mesh.potIdx(j, k, sLoc[i])];
              const pe = uPot[mesh.potIdx(j, k, eLoc[i])];
              const Tloc = this.T[c];
              const frt = F_CONST / (R_GAS * Tloc);
              const spec = reg[i] === REGION_IDS.NEG ? negSpec : posSpec;
              const csMax = spec.csMax;
              const cSurf = Math.max(10.0, Math.min(csMax - 10.0, this.cse[c]));
              const ceLoc = Math.max(50.0, this.ce[c]);
              const kr = reactionRate(spec, Tloc);
              const i0 = F_CONST * kr * Math.sqrt(ceLoc * cSurf * (csMax - cSurf));
              const i0a = mesh.asArea[i] * i0;

              const UeqTilde =
                reg[i] === REGION_IDS.NEG
                  ? uNeg(cSurf / csMax)
                  : uPos(cSurf / csMax) - Vref;

              const eta = psTilde - pe - UeqTilde;
              const arg = Math.max(-30.0, Math.min(30.0, 0.5 * frt * eta));
              const jl = 2.0 * i0a * Math.sinh(arg);
              const dj = i0a * frt * Math.cosh(arg);

              resPot[mesh.potIdx(j, k, sLoc[i])] += vol * jl;
              resPot[mesh.potIdx(j, k, eLoc[i])] -= vol * jl;
              gbv[c] = vol * dj;
              this.jLi[c] = jl;
              this.etaField[c] = eta;
            }
          }
        }
      }

      let rNorm2 = 0.0;
      for (let q = 0; q < mesh.Npot; q++) {
        rNorm2 += resPot[q] * resPot[q];
        rhsPot[q] = -resPot[q];
      }
      rNorm2 = Math.sqrt(rNorm2);
      newtonHistory.push(rNorm2);

      if (nit === 1 && this.stepCount === 0) {
        firstNewtonRhsSnapshot = new Float64Array(rhsPot);
        firstNewtonGbvSnapshot = new Float64Array(gbv);
      }

      if (rNorm2 < 1e-12 * Math.max(1e-6, I_total)) {
        break;
      }

      coupledMG.updateOperator(kap, gbv);

      let linRes;
      if (config.solverMode === 'xline_mg_pcg') {
        linRes = coupledMG.solvePCG(rhsPot, duPot, 1e-10, 80);
      } else {
        linRes = hypSolver.solveNewtonStep(rhsPot, gbv, duPot, 1e-10, 80);
      }
      totalLinIters += linRes.iters;
      lastLinIters = linRes.iters;
      lastLinRel = linRes.relRes;

      let maxDu = 0.0;
      for (let q = 0; q < mesh.Npot; q++) {
        const av = Math.abs(duPot[q]);
        if (av > maxDu) maxDu = av;
      }

      const damp = maxDu > 0.25 ? 0.25 / maxDu : 1.0;
      for (let q = 0; q < mesh.Npot; q++) {
        uPot[q] += damp * duPot[q];
      }

      if (maxDu < 1e-11) {
        break;
      }
    }

    // 5. Unpack physical potentials & compute conservative 3D current vectors + Joule heat
    for (let j = 0; j < ny; j++) {
      for (let k = 0; k < nz; k++) {
        for (let i = 0; i < nx; i++) {
          const c = mesh.cid(i, j, k);
          if (sLoc[i] >= 0) {
            const shift = reg[i] >= REGION_IDS.POS ? Vref : 0.0;
            this.phiS[c] = uPot[mesh.potIdx(j, k, sLoc[i])] + shift;
          }
          if (eLoc[i] >= 0) {
            this.phiE[c] = uPot[mesh.potIdx(j, k, eLoc[i])];
          }
        }
      }
    }

    hypSolver.reconstructCurrentsAndHeating(uPot, kap, this.isVec, this.ieVec, this.qOhm);

    // 6. Compute terminal voltage, SOC, charge conservation, and region extrema (Table 18 format)
    let vPosSum = 0.0;
    for (const [j, k] of mesh.posTabCells) {
      vPosSum += this.phiS[mesh.cid(nx - 1, j, k)];
    }
    const vPosAvg = vPosSum / mesh.posTabCells.length;

    let vNegSum = 0.0;
    for (const [j, k] of mesh.negTabCells) {
      const iRef = negTabType === 'top_edge' ? Math.floor(mesh.nxCounts[0] / 2) : 0;
      vNegSum += this.phiS[mesh.cid(iRef, j, k)];
    }
    const vNegAvg = vNegSum / mesh.negTabCells.length;
    const Vcell = vPosAvg - vNegAvg;

    let I_neg = 0.0;
    let I_pos = 0.0;
    let negVol = 0.0;
    let csNegAvg = 0.0;
    for (let i = 0; i < nx; i++) {
      const vol = mesh.cellVol[i];
      if (reg[i] === REGION_IDS.NEG) {
        for (let j = 0; j < ny; j++) {
          for (let k = 0; k < nz; k++) {
            const c = mesh.cid(i, j, k);
            I_neg += this.jLi[c] * vol;
            csNegAvg += this.cs[c] * vol;
            negVol += vol;
          }
        }
      } else if (reg[i] === REGION_IDS.POS) {
        for (let j = 0; j < ny; j++) {
          for (let k = 0; k < nz; k++) {
            I_pos -= this.jLi[mesh.cid(i, j, k)] * vol;
          }
        }
      }
    }
    csNegAvg /= negVol;
    const soc = Math.max(0, Math.min(1, (csNegAvg / negSpec.csMax) / negSpec.theta0));
    const chargeErr =
      I_total > 0 ? Math.max(Math.abs(I_neg - I_total), Math.abs(I_pos - I_total)) / I_total : 0;

    const elapsedMs = performance.now() - tSolveStart;

    if (firstNewtonRhsSnapshot) {
      this.firstNewtonBenchmark = {
        kap: new Float64Array(kap),
        gbv: firstNewtonGbvSnapshot,
        rhsPhi: firstNewtonRhsSnapshot,
        R0,
      };
    }

    this.lastDiagnostics = {
      solverMode: config.solverMode,
      newtonIters: newtonHistory.length,
      totalLinIters,
      lastLinIters,
      lastLinRel,
      newtonHistory,
      solveMs: elapsedMs,
      Vcell,
      soc,
      I_total,
      I_neg,
      I_pos,
      chargeErr,
      regionStats: this.computeRegionStats(),
      ...this.computeThermalExtrema(),
    };
    return this.lastDiagnostics;
  }

  /**
   * Advance the full 6-equation 3D coupled system by dt seconds.
   */
  stepTime(dt) {
    if (this.terminated) {
      return this.lastDiagnostics;
    }
    const { mesh } = this;
    const { nx, ny, nz, reg } = mesh;
    const negSpec = LAYER_SPECS[REGION_IDS.NEG];
    const posSpec = LAYER_SPECS[REGION_IDS.POS];

    // 1. Update solid lithium concentration (Eq. 3)
    let hitStoichLimit = false;
    for (let i = 0; i < nx; i++) {
      const r = reg[i];
      if (r !== REGION_IDS.NEG && r !== REGION_IDS.POS) continue;
      const spec = r === REGION_IDS.NEG ? negSpec : posSpec;
      const lse = spec.Rs / 5.0; // Wang-Gu-Liaw parabolic diffusion length (Report Eq. 3)
      const invEpsF = 1.0 / (mesh.epsS[i] * F_CONST);
      const asF = mesh.asArea[i] * F_CONST;

      for (let j = 0; j < ny; j++) {
        for (let k = 0; k < nz; k++) {
          const c = mesh.cid(i, j, k);
          const jl = this.jLi[c];
          const Ds = solidDiffusivity(spec, this.T[c]);
          let csNew = this.cs[c] - dt * invEpsF * jl;
          let cseNew = csNew - (lse / (Ds * asF)) * jl;

          if (r === REGION_IDS.NEG) {
            if (cseNew <= 0.015 * spec.csMax) hitStoichLimit = true;
            csNew = Math.max(0.01 * spec.csMax, Math.min(0.99 * spec.csMax, csNew));
            cseNew = Math.max(0.008 * spec.csMax, Math.min(0.99 * spec.csMax, cseNew));
          } else {
            if (cseNew >= 0.985 * spec.csMax) hitStoichLimit = true;
            csNew = Math.max(0.05 * spec.csMax, Math.min(0.99 * spec.csMax, csNew));
            cseNew = Math.max(0.05 * spec.csMax, Math.min(0.992 * spec.csMax, cseNew));
          }
          this.cs[c] = csNew;
          this.cse[c] = cseNew;
        }
      }
    }

    // 2. Implicit Backward Euler 3D Electrolyte Species Transport (Eqs. 4 & 9)
    this._stepElectrolyteSpecies(dt);

    // 3. Implicit Backward Euler 3D Temperature (Eqs. 8, 11, 13)
    this._stepTemperature(dt);

    // 4. Re-solve coupled (phi_s, phi_e, j^Li) at updated (c_se, c_e, T)
    this.time += dt;
    this.stepCount += 1;
    const diag = this.solvePotentials();

    if (diag.Vcell <= V_CUTOFF || hitStoichLimit) {
      this.terminated = true;
      this.terminationReason =
        diag.Vcell <= V_CUTOFF ? `Reached cutoff voltage (${V_CUTOFF.toFixed(2)} V)` : 'Active material stoichiometry limit reached';
    }

    this._recordHistory(diag);
    return diag;
  }

  _stepElectrolyteSpecies(dt) {
    const { mesh, ceMG, ceX, ceRhs, ceDiag, ceGx, ceGy, ceGz } = this;
    const { nxElec, ny, nz, dx, dy, dz, eXIndices, reg, epsE } = mesh;
    const tTerm = (1.0 - ELECTROLYTE_PARAMS.tPlus) / F_CONST;
    const areaX = dy * dz;

    ceDiag.fill(0);
    ceGx.fill(0);
    ceGy.fill(0);
    ceGz.fill(0);

    // Local effective diffusivity per active electrolyte cell
    const De = new Float64Array(nxElec * ny * nz);
    for (let ie = 0; ie < nxElec; ie++) {
      const i = eXIndices[ie];
      const ee = epsE[i];
      for (let j = 0; j < ny; j++) {
        for (let k = 0; k < nz; k++) {
          const c = mesh.cid(i, j, k);
          De[mesh.elecIdx(ie, j, k)] = electrolyteDiffusivity(ee, this.T[c]);
        }
      }
    }

    for (let j = 0; j < ny; j++) {
      for (let k = 0; k < nz; k++) {
        const base = (j * nz + k) * nxElec;
        for (let ie = 0; ie < nxElec; ie++) {
          const i = eXIndices[ie];
          const c = mesh.cid(i, j, k);
          const vol = mesh.cellVol[i];
          const cap = (epsE[i] * vol) / dt;
          const idx = base + ie;

          ceX[idx] = this.ce[c];
          ceDiag[idx] = cap;
          const src = reg[i] === REGION_IDS.NEG || reg[i] === REGION_IDS.POS ? tTerm * vol * this.jLi[c] : 0.0;
          ceRhs[idx] = cap * this.ce[c] + src;

          if (ie > 0) {
            const ip = eXIndices[ie - 1];
            const dL = De[idx - 1];
            const dR = De[idx];
            ceGx[idx] = areaX / ((0.5 * dx[ip]) / dL + (0.5 * dx[i]) / dR);
          }
        }
      }
    }

    for (let ie = 0; ie < nxElec; ie++) {
      const i = eXIndices[ie];
      for (let j = 0; j < ny - 1; j++) {
        for (let k = 0; k < nz; k++) {
          const dL = De[mesh.elecIdx(ie, j, k)];
          const dR = De[mesh.elecIdx(ie, j + 1, k)];
          const deff = (2.0 * dL * dR) / (dL + dR);
          ceGy[(j * nz + k) * nxElec + ie] = (deff * dx[i] * dz) / dy;
        }
      }
      for (let j = 0; j < ny; j++) {
        for (let k = 0; k < nz - 1; k++) {
          const dL = De[mesh.elecIdx(ie, j, k)];
          const dR = De[mesh.elecIdx(ie, j, k + 1)];
          const deff = (2.0 * dL * dR) / (dL + dR);
          ceGz[(j * (nz - 1) + k) * nxElec + ie] = (deff * dx[i] * dy) / dz;
        }
      }
    }

    ceMG.buildFromConductances(ceDiag, ceGx, ceGy, ceGz);
    ceMG.solvePCG(ceRhs, ceX, 1e-10, 40);

    for (let ie = 0; ie < nxElec; ie++) {
      const i = eXIndices[ie];
      for (let j = 0; j < ny; j++) {
        for (let k = 0; k < nz; k++) {
          const v = Math.max(50.0, ceX[mesh.elecIdx(ie, j, k)]);
          this.ce[mesh.cid(i, j, k)] = v;
        }
      }
    }
  }

  _stepTemperature(dt) {
    const { mesh, config, tempMG, tX, tRhs, tDiag, tGx, tGy, tGz } = this;
    const { nx, ny, nz, dx, dy, dz, reg, rhoCp, lambda } = mesh;
    const areaX = dy * dz;
    const h = config.thermalBC === 'convective' ? config.hConv : 0.0;
    const isTopAmbient = config.thermalBC === 'top_ambient';

    tDiag.fill(0);
    tGx.fill(0);
    tGy.fill(0);
    tGz.fill(0);

    // Precompute volumetric heat generation qTotal = qOhm + qRxn + qRev [W/m^3]
    for (let i = 0; i < nx; i++) {
      const r = reg[i];
      const invCsMax =
        r === REGION_IDS.NEG
          ? 1.0 / LAYER_SPECS[REGION_IDS.NEG].csMax
          : r === REGION_IDS.POS
          ? 1.0 / LAYER_SPECS[REGION_IDS.POS].csMax
          : 0.0;
      for (let j = 0; j < ny; j++) {
        for (let k = 0; k < nz; k++) {
          const c = mesh.cid(i, j, k);
          let q = this.qOhm[c];
          if (r === REGION_IDS.NEG || r === REGION_IDS.POS) {
            const jl = this.jLi[c];
            const qRxn = jl * this.etaField[c]; // >= 0 irreversible
            const th = this.cse[c] * invCsMax;
            const dUdT = r === REGION_IDS.NEG ? dUNegDT(th) : dUPosDT(th);
            const qRev = jl * this.T[c] * dUdT;
            q += Math.max(0.0, qRxn) + qRev;
          }
          this.qTotal[c] = q;
        }
      }
    }

    for (let j = 0; j < ny; j++) {
      for (let k = 0; k < nz; k++) {
        const base = (j * nz + k) * nx;
        for (let i = 0; i < nx; i++) {
          const c = mesh.cid(i, j, k);
          const vol = mesh.cellVol[i];
          const cap = (rhoCp[i] * vol) / dt;
          const idx = base + i;

          tX[idx] = this.T[c];
          let diag = cap;
          let rhs = cap * this.T[c] + vol * this.qTotal[c];

          if (i > 0) {
            tGx[idx] = areaX / ((0.5 * dx[i - 1]) / lambda[i - 1] + (0.5 * dx[i]) / lambda[i]);
          }

          // External boundary heat transfer
          if (h > 0) {
            if (i === 0 || i === nx - 1) {
              const gExt = 1.0 / (1.0 / h + (0.5 * dx[i]) / lambda[i]) * areaX;
              diag += gExt;
              rhs += gExt * T_AMBIENT;
            }
            if (j === 0 || j === ny - 1) {
              const gExt = 1.0 / (1.0 / h + (0.5 * dy) / lambda[i]) * (dx[i] * dz);
              diag += gExt;
              rhs += gExt * T_AMBIENT;
            }
            if (k === 0 || k === nz - 1) {
              const gExt = 1.0 / (1.0 / h + (0.5 * dz) / lambda[i]) * (dx[i] * dy);
              diag += gExt;
              rhs += gExt * T_AMBIENT;
            }
          } else if (isTopAmbient && j === ny - 1) {
            // Report Fig. 9: top face clamped to T_AMBIENT, all other boundaries adiabatic
            const gTop = (lambda[i] * dx[i] * dz) / (0.5 * dy);
            diag += gTop;
            rhs += gTop * T_AMBIENT;
          }

          tDiag[idx] = diag;
          tRhs[idx] = rhs;
        }
      }
    }

    for (let i = 0; i < nx; i++) {
      const gyVal = (lambda[i] * dx[i] * dz) / dy;
      const gzVal = (lambda[i] * dx[i] * dy) / dz;
      for (let j = 0; j < ny - 1; j++) {
        for (let k = 0; k < nz; k++) tGy[(j * nz + k) * nx + i] = gyVal;
      }
      for (let j = 0; j < ny; j++) {
        for (let k = 0; k < nz - 1; k++) tGz[(j * (nz - 1) + k) * nx + i] = gzVal;
      }
    }

    tempMG.buildFromConductances(tDiag, tGx, tGy, tGz);
    tempMG.solvePCG(tRhs, tX, 1e-10, 40);

    for (let i = 0; i < nx; i++) {
      for (let j = 0; j < ny; j++) {
        for (let k = 0; k < nz; k++) {
          this.T[mesh.cid(i, j, k)] = tX[mesh.thermIdx(i, j, k)];
        }
      }
    }
  }

  computeThermalExtrema() {
    let tMin = Infinity;
    let tMax = -Infinity;
    for (let c = 0; c < this.T.length; c++) {
      if (this.T[c] < tMin) tMin = this.T[c];
      if (this.T[c] > tMax) tMax = this.T[c];
    }
    return { Tmin: tMin, Tmax: tMax };
  }

  /**
   * Region-by-region extrema matching Report Table 18 (page 30) with jLi in A/cm^3.
   */
  computeRegionStats() {
    const { mesh, phiS, phiE, jLi, T, ce } = this;
    const { nx, ny, nz, reg } = mesh;
    const stats = REGION_SHORT.map((name, r) => ({
      region: name,
      regionId: r,
      minPhiS: Infinity,
      maxPhiS: -Infinity,
      hasPhiS: r !== REGION_IDS.SEP,
      minPhiE: Infinity,
      maxPhiE: -Infinity,
      hasPhiE: r >= REGION_IDS.NEG && r <= REGION_IDS.POS,
      minJLi: Infinity,
      maxJLi: -Infinity,
      hasJLi: r === REGION_IDS.NEG || r === REGION_IDS.POS,
      minT: Infinity,
      maxT: -Infinity,
      minCe: Infinity,
      maxCe: -Infinity,
    }));

    for (let i = 0; i < nx; i++) {
      const st = stats[reg[i]];
      for (let j = 0; j < ny; j++) {
        for (let k = 0; k < nz; k++) {
          const c = mesh.cid(i, j, k);
          if (st.hasPhiS) {
            if (phiS[c] < st.minPhiS) st.minPhiS = phiS[c];
            if (phiS[c] > st.maxPhiS) st.maxPhiS = phiS[c];
          }
          if (st.hasPhiE) {
            if (phiE[c] < st.minPhiE) st.minPhiE = phiE[c];
            if (phiE[c] > st.maxPhiE) st.maxPhiE = phiE[c];
            if (ce[c] < st.minCe) st.minCe = ce[c];
            if (ce[c] > st.maxCe) st.maxCe = ce[c];
          }
          if (st.hasJLi) {
            const jCm3 = jLi[c] * 1e-6; // A/m^3 -> A/cm^3 (Report Table 18 units!)
            if (jCm3 < st.minJLi) st.minJLi = jCm3;
            if (jCm3 > st.maxJLi) st.maxJLi = jCm3;
          }
          if (T[c] < st.minT) st.minT = T[c];
          if (T[c] > st.maxT) st.maxT = T[c];
        }
      }
    }
    return stats;
  }

  /**
   * Extract 1D through-thickness profile along x at fractional (yFrac, zFrac)
   * for the line plot (matching Report Figs. 4, 5, 26, 27).
   */
  extractProfileX(yFrac = 0.5, zFrac = 0.5) {
    const { mesh, phiS, phiE, jLi, ce, T, cs } = this;
    const j = Math.max(0, Math.min(mesh.ny - 1, Math.floor(yFrac * mesh.ny)));
    const k = Math.max(0, Math.min(mesh.nz - 1, Math.floor(zFrac * mesh.nz)));
    const xsUm = new Float64Array(mesh.nx);
    const ps = new Float64Array(mesh.nx);
    const pe = new Float64Array(mesh.nx);
    const jl = new Float64Array(mesh.nx);
    const cArr = new Float64Array(mesh.nx);
    const tArr = new Float64Array(mesh.nx);
    const socArr = new Float64Array(mesh.nx);

    for (let i = 0; i < mesh.nx; i++) {
      const c = mesh.cid(i, j, k);
      const r = mesh.reg[i];
      xsUm[i] = mesh.xCenters[i] * 1e6;
      ps[i] = mesh.sLoc[i] >= 0 ? phiS[c] : NaN;
      pe[i] = mesh.eLoc[i] >= 0 ? phiE[c] : NaN;
      jl[i] = r === REGION_IDS.NEG || r === REGION_IDS.POS ? jLi[c] * 1e-6 : NaN;
      cArr[i] = mesh.eLoc[i] >= 0 ? ce[c] : NaN;
      tArr[i] = T[c];
      socArr[i] =
        r === REGION_IDS.NEG
          ? cs[c] / LAYER_SPECS[REGION_IDS.NEG].csMax
          : r === REGION_IDS.POS
          ? cs[c] / LAYER_SPECS[REGION_IDS.POS].csMax
          : NaN;
    }

    return {
      xUm: xsUm,
      reg: mesh.reg,
      phiS: ps,
      phiE: pe,
      jLi: jl,
      ce: cArr,
      T: tArr,
      theta: socArr,
    };
  }

  _recordHistory(diag) {
    this.history.push({
      t: this.time,
      V: diag.Vcell,
      Tmin: diag.Tmin,
      Tmax: diag.Tmax,
      soc: diag.soc,
      linIters: diag.lastLinIters,
      newtonIters: diag.newtonIters,
      chargeErr: diag.chargeErr,
    });
  }

  /**
   * Package lightweight typed-array snapshot for transferring across Web Worker to 3D UI.
   */
  getSnapshot(yFrac = 0.5, zFrac = 0.5) {
    const { mesh } = this;
    const negCsMax = LAYER_SPECS[REGION_IDS.NEG].csMax;
    const posCsMax = LAYER_SPECS[REGION_IDS.POS].csMax;
    const stoich = new Float64Array(mesh.N);
    for (let i = 0; i < mesh.nx; i++) {
      const r = mesh.reg[i];
      if (r === REGION_IDS.NEG || r === REGION_IDS.POS) {
        const invMax = 1.0 / (r === REGION_IDS.NEG ? negCsMax : posCsMax);
        for (let j = 0; j < mesh.ny; j++) {
          for (let k = 0; k < mesh.nz; k++) {
            const c = mesh.cid(i, j, k);
            stoich[c] = this.cs[c] * invMax;
          }
        }
      }
    }

    return {
      time: this.time,
      stepCount: this.stepCount,
      terminated: this.terminated,
      terminationReason: this.terminationReason,
      diagnostics: this.lastDiagnostics,
      history: this.history,
      profile: this.extractProfileX(yFrac, zFrac),
      meshMeta: {
        nx: mesh.nx,
        ny: mesh.ny,
        nz: mesh.nz,
        Lx: mesh.Lx,
        Ly: mesh.Ly,
        Lz: mesh.Lz,
        bc: mesh.bc,
        N: mesh.N,
        Npot: mesh.Npot,
        nFaces: mesh.nFaces,
        xEdges: mesh.xEdges,
        xCenters: mesh.xCenters,
        dx: mesh.dx,
        reg: mesh.reg,
        negTabMask: mesh.negTabMask,
        posTabMask: mesh.posTabMask,
        negTabType: mesh.negTabType,
        nxCounts: mesh.nxCounts,
      },
      fields: {
        phiS: this.phiS,
        phiE: this.phiE,
        jLi: this.jLi,
        ce: this.ce,
        T: this.T,
        stoich,
        isVec: this.isVec,
        ieVec: this.ieVec,
        qTotal: this.qTotal,
      },
    };
  }
}
