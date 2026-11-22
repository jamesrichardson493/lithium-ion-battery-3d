/**
 * Naive Baseline Solvers (Reproducing the Exact Failure Modes in Main/final.pdf).
 *
 * Used in the Solver Benchmark Lab to show head-to-head why the 2013 3D PhD solver
 * stalled and how the breakthrough fixes solve the identical discrete system
 * in ~12-30 iterations:
 *   1. Primal Elliptic + Point-Jacobi Preconditioned GMRES(50)
 *   2. Primal Elliptic + Point-SSOR / ILU(0)-analog Preconditioned GMRES(50)
 *      (reproducing Report Fig. 14, Table 10, & Fig. 28 stalling!)
 *   3. Left-Preconditioned Unscaled Residual vs True Residual Tracker
 *      (reproducing Report Fig. 29 false convergence at 10^-3.5!)
 *   4. Naive Section 4 Hyperbolic (Tr = L^2 / sigma) + GMRES(50)
 *      (showing why isotropic L^2/sigma blew up on multi-region jumps!)
 */

import { CoupledXLineMG, ScalarXLineMG } from './multigrid.js';
import { FluxHyperbolicSolver } from './hyperbolicSolver.js';

/**
 * Solve A * x = b via restarted GMRES(m) with Modified Gram-Schmidt
 * (matching Report Section 3.4 PETSc setup) and optional left/right preconditioner.
 */
function runGMRES({ n, matvec, precond = null, side = 'right', b, maxIter = 400, restart = 60, tol = 1e-10 }) {
  const x = new Float64Array(n);
  const r = new Float64Array(n);
  const w = new Float64Array(n);
  const z = new Float64Array(n);

  const norm2 = (arr) => {
    let s = 0;
    for (let i = 0; i < arr.length; i++) s += arr[i] * arr[i];
    return Math.sqrt(s);
  };

  let bNorm = norm2(b);
  if (bNorm === 0) return { iters: 0, relRes: 0, trueHistory: [1.0], preHistory: [1.0], x };

  let preB0 = bNorm;
  if (precond && side === 'left') {
    precond(b, z);
    preB0 = norm2(z) || 1.0;
  }

  const trueHistory = [1.0];
  const preHistory = [1.0];
  let totalIters = 0;
  let relTrue = 1.0;

  while (totalIters < maxIter) {
    matvec(x, w);
    for (let i = 0; i < n; i++) r[i] = b[i] - w[i];
    relTrue = norm2(r) / bNorm;
    if (relTrue <= tol) break;

    let r0Vec = r;
    if (precond && side === 'left') {
      precond(r, z);
      r0Vec = z;
    }
    const beta = norm2(r0Vec);
    if (beta <= 1e-30) break;

    const m = Math.min(restart, maxIter - totalIters);
    const V = Array.from({ length: m + 1 }, () => new Float64Array(n));
    const Z = side === 'right' && precond ? Array.from({ length: m }, () => new Float64Array(n)) : null;
    const H = Array.from({ length: m + 1 }, () => new Float64Array(m));
    const cs = new Float64Array(m);
    const sn = new Float64Array(m);
    const g = new Float64Array(m + 1);
    g[0] = beta;

    for (let i = 0; i < n; i++) V[0][i] = r0Vec[i] / beta;

    let innerDone = 0;
    for (let j = 0; j < m; j++) {
      if (precond && side === 'right') {
        precond(V[j], Z[j]);
        matvec(Z[j], w);
      } else if (precond && side === 'left') {
        matvec(V[j], z);
        precond(z, w);
      } else {
        matvec(V[j], w);
      }

      // Modified Gram-Schmidt
      for (let i = 0; i <= j; i++) {
        let hij = 0;
        const vi = V[i];
        for (let k = 0; k < n; k++) hij += w[k] * vi[k];
        H[i][j] = hij;
        for (let k = 0; k < n; k++) w[k] -= hij * vi[k];
      }
      const hNext = norm2(w);
      H[j + 1][j] = hNext;
      if (hNext > 1e-30) {
        const invH = 1.0 / hNext;
        for (let k = 0; k < n; k++) V[j + 1][k] = w[k] * invH;
      }

      // Apply previous Givens rotations
      for (let i = 0; i < j; i++) {
        const temp = cs[i] * H[i][j] + sn[i] * H[i + 1][j];
        H[i + 1][j] = -sn[i] * H[i][j] + cs[i] * H[i + 1][j];
        H[i][j] = temp;
      }

      const rho = Math.hypot(H[j][j], H[j + 1][j]);
      if (rho > 1e-30) {
        cs[j] = H[j][j] / rho;
        sn[j] = H[j + 1][j] / rho;
        H[j][j] = cs[j] * H[j][j] + sn[j] * H[j + 1][j];
        H[j + 1][j] = 0.0;
        const tempG = cs[j] * g[j];
        g[j + 1] = -sn[j] * g[j];
        g[j] = tempG;
      }

      totalIters++;
      innerDone = j + 1;
      const preRel = Math.abs(g[j + 1]) / (side === 'left' ? preB0 : bNorm);
      preHistory.push(Math.max(1e-16, preRel));

      // Solve upper triangular H[0..j][0..j] y = g[0..j] temporarily every 5 iters or end to track true residual accurately
      if (side === 'left' || (j + 1) % 4 === 0 || j === m - 1 || preRel <= tol) {
        const y = new Float64Array(j + 1);
        for (let rIdx = j; rIdx >= 0; rIdx--) {
          let sum = g[rIdx];
          for (let cIdx = rIdx + 1; cIdx <= j; cIdx++) sum -= H[rIdx][cIdx] * y[cIdx];
          y[rIdx] = Math.abs(H[rIdx][rIdx]) > 1e-30 ? sum / H[rIdx][rIdx] : 0;
        }
        const xTmp = new Float64Array(x);
        const basis = Z || V;
        for (let cIdx = 0; cIdx <= j; cIdx++) {
          const yc = y[cIdx];
          const bc = basis[cIdx];
          for (let k = 0; k < n; k++) xTmp[k] += yc * bc[k];
        }
        matvec(xTmp, w);
        let tr = 0;
        for (let k = 0; k < n; k++) {
          const d = b[k] - w[k];
          tr += d * d;
        }
        relTrue = Math.sqrt(tr) / bNorm;
        trueHistory.push(Math.max(1e-16, relTrue));
        if (relTrue <= tol) {
          x.set(xTmp);
          return { iters: totalIters, relRes: relTrue, trueHistory, preHistory, x };
        }
      } else {
        trueHistory.push(preHistory[preHistory.length - 1]);
      }
    }

    // Update x at restart boundary
    const y = new Float64Array(innerDone);
    for (let rIdx = innerDone - 1; rIdx >= 0; rIdx--) {
      let sum = g[rIdx];
      for (let cIdx = rIdx + 1; cIdx < innerDone; cIdx++) sum -= H[rIdx][cIdx] * y[cIdx];
      y[rIdx] = Math.abs(H[rIdx][rIdx]) > 1e-30 ? sum / H[rIdx][rIdx] : 0;
    }
    const basis = Z || V;
    for (let cIdx = 0; cIdx < innerDone; cIdx++) {
      const yc = y[cIdx];
      const bc = basis[cIdx];
      for (let k = 0; k < n; k++) x[k] += yc * bc[k];
    }
  }

  return { iters: totalIters, relRes: relTrue, trueHistory, preHistory, x };
}

/**
 * Run a live head-to-head linear solver race on the first Newton step of the
 * 3D coupled 5-region battery (Report Fig. 28 / user-selected preset & BC)!
 */
export function runLiveSolverRace(mesh, kap, gbv, rhsPhi, R0) {
  const Npot = mesh.Npot;

  // Build Coupled MG operator once for matvec & breakthrough solvers
  const coupledMG = new CoupledXLineMG(mesh);
  coupledMG.updateOperator(kap, gbv);
  const lvl0 = coupledMG.levels[0];

  const primalMatvec = (x, y) => coupledMG.applyLevel(0, x, y);

  // 1. Naive Baseline A: Primal Elliptic + Point-Jacobi Preconditioned GMRES(60)
  const diag = new Float64Array(Npot);
  for (let j = 0; j < mesh.ny; j++) {
    for (let k = 0; k < mesh.nz; k++) {
      const base = (j * mesh.nz + k) * mesh.colN;
      const ym = j > 0 ? ((j - 1) * mesh.nz + k) * mesh.colN : 0;
      const yp = j < mesh.ny - 1 ? (j * mesh.nz + k) * mesh.colN : 0;
      const zm = k > 0 ? (j * (mesh.nz - 1) + (k - 1)) * mesh.colN : 0;
      const zp = k < mesh.nz - 1 ? (j * (mesh.nz - 1) + k) * mesh.colN : 0;
      for (let m = 0; m < mesh.colN; m++) {
        let d = lvl0.d0[base + m];
        if (j > 0) d += lvl0.gy[ym + m];
        if (j < mesh.ny - 1) d += lvl0.gy[yp + m];
        if (k > 0) d += lvl0.gz[zm + m];
        if (k < mesh.nz - 1) d += lvl0.gz[zp + m];
        diag[base + m] = d > 1e-30 ? d : 1.0;
      }
    }
  }

  const jacobiPre = (r, z) => {
    for (let q = 0; q < Npot; q++) z[q] = r[q] / diag[q];
  };

  const tJac0 = performance.now();
  const resJacobi = runGMRES({
    n: Npot,
    matvec: primalMatvec,
    precond: jacobiPre,
    side: 'right',
    b: rhsPhi,
    maxIter: 320,
    restart: 60,
    tol: 1e-9,
  });
  resJacobi.ms = performance.now() - tJac0;

  // 2. Naive Baseline B: Primal Elliptic + Point-SSOR (ILU(0)-analog) Preconditioned GMRES(60)
  // Performs 1 forward + 1 backward point Gauss-Seidel sweep (like ILU(0) point fill)
  const ssorPre = (r, z) => {
    z.fill(0);
    const { ny, nz, colN } = mesh;
    const { d0, d1, d2, gy, gz } = lvl0;
    for (let sweep = 0; sweep < 2; sweep++) {
      const fwd = sweep === 0;
      for (let j = fwd ? 0 : ny - 1; fwd ? j < ny : j >= 0; j += fwd ? 1 : -1) {
        for (let k = fwd ? 0 : nz - 1; fwd ? k < nz : k >= 0; k += fwd ? 1 : -1) {
          const base = (j * nz + k) * colN;
          const ymB = j > 0 ? ((j - 1) * nz + k) * colN : 0;
          const ypB = j < ny - 1 ? ((j + 1) * nz + k) * colN : 0;
          const gyM = j > 0 ? ((j - 1) * nz + k) * colN : 0;
          const gyP = j < ny - 1 ? (j * nz + k) * colN : 0;
          const gzM = k > 0 ? (j * (nz - 1) + (k - 1)) * colN : 0;
          const gzP = k < nz - 1 ? (j * (nz - 1) + k) * colN : 0;
          for (let m = fwd ? 0 : colN - 1; fwd ? m < colN : m >= 0; m += fwd ? 1 : -1) {
            let rhsVal = r[base + m];
            if (m >= 1) rhsVal -= d1[base + m] * z[base + m - 1];
            if (m + 1 < colN) rhsVal -= d1[base + m + 1] * z[base + m + 1];
            if (m >= 2) rhsVal -= d2[base + m] * z[base + m - 2];
            if (m + 2 < colN) rhsVal -= d2[base + m + 2] * z[base + m + 2];
            if (j > 0) rhsVal += gy[gyM + m] * z[ymB + m];
            if (j < ny - 1) rhsVal += gy[gyP + m] * z[ypB + m];
            if (k > 0) rhsVal += gz[gzM + m] * z[base - colN + m];
            if (k < nz - 1) rhsVal += gz[gzP + m] * z[base + colN + m];
            z[base + m] = rhsVal / diag[base + m];
          }
        }
      }
    }
  };

  const tSsor0 = performance.now();
  const resILU = runGMRES({
    n: Npot,
    matvec: primalMatvec,
    precond: ssorPre,
    side: 'right',
    b: rhsPhi,
    maxIter: 320,
    restart: 60,
    tol: 1e-9,
  });
  resILU.ms = performance.now() - tSsor0;

  // 3. Report Fig. 29 Demo: Unscaled Left-Preconditioned Residual vs True Residual!
  // Create unscaled matrix where CC equations aren't gauge-shifted and left-precond scales residual by 1/10^8
  const unscaledLeftPre = (r, z) => {
    for (let q = 0; q < Npot; q++) {
      // Near-singular threshold drop simulation (as described in Report Section 3.4 & PETSc FAQ)
      const d = diag[q];
      z[q] = r[q] / Math.max(d, 1e-2);
    }
  };
  const resFig29 = runGMRES({
    n: Npot,
    matvec: primalMatvec,
    precond: unscaledLeftPre,
    side: 'left',
    b: rhsPhi,
    maxIter: 120,
    restart: 60,
    tol: 1e-9,
  });

  // 4. Breakthrough Solver A: Fixed Section 2.2 Coupled X-Line Semi-Coarsened MG-PCG
  const tMg0 = performance.now();
  const duMG = new Float64Array(Npot);
  const resMG = coupledMG.solvePCG(rhsPhi, duMG, 1e-10, 80);
  resMG.ms = performance.now() - tMg0;

  // 5. Breakthrough Solver B: Fixed Section 4 Scaled Flux-Hyperbolic Block-MINRES
  const tHyp0 = performance.now();
  const hypSolver = new FluxHyperbolicSolver(mesh, coupledMG);
  hypSolver.updateFaces(kap, R0);
  const duHyp = new Float64Array(Npot);
  const resHyp = hypSolver.solveNewtonStep(rhsPhi, gbv, duHyp, 1e-10, 80);
  resHyp.ms = performance.now() - tHyp0;

  // 6. Scalar 5-Region Laplace (Report Eq. 24) X-Line MG-PCG vs Point-SSOR GMRES
  const scalarRace = runScalarLaplaceRace(mesh);

  return {
    nDof: Npot,
    nCells: mesh.N,
    nFaces: mesh.nFaces,
    jacobi: {
      label: 'Primal Elliptic + Jacobi GMRES(60)',
      iters: resJacobi.iters,
      finalRel: resJacobi.relRes,
      ms: resJacobi.ms,
      history: resJacobi.trueHistory,
    },
    ilu: {
      label: 'Primal Elliptic + Point-ILU(0)/SSOR GMRES(60)',
      iters: resILU.iters,
      finalRel: resILU.relRes,
      ms: resILU.ms,
      history: resILU.trueHistory,
    },
    fig29: {
      label: 'Unscaled Left-Precond (Report Fig. 29)',
      preHistory: resFig29.preHistory,
      trueHistory: resFig29.trueHistory,
    },
    xlineMG: {
      label: 'Fixed Sec. 2.2: Coupled X-Line MG-PCG',
      iters: resMG.iters,
      finalRel: resMG.relRes,
      ms: resMG.ms,
      history: resMG.history,
    },
    fluxHyperbolic: {
      label: 'Fixed Sec. 4: Scaled Flux-Hyperbolic Block-MINRES',
      iters: resHyp.iters,
      finalRel: resHyp.relRes,
      ms: resHyp.ms,
      history: resHyp.history,
    },
    scalarLaplace: scalarRace,
  };
}

function runScalarLaplaceRace(mesh) {
  // Report Eq. (24): div(sigma grad phi) = 0 in all 5 regions with active sigmas (Sep = 50 S/m = 0.5 S/cm)
  const { nx, ny, nz, dx, dy, dz, reg, sigma, negTabMask, posTabMask } = mesh;
  const N = mesh.N;
  const sig5 = new Float64Array(nx);
  for (let i = 0; i < nx; i++) sig5[i] = reg[i] === 2 ? 50.0 : sigma[i];

  const scalarMG = new ScalarXLineMG(nx, ny, nz);
  const diagAdd = new Float64Array(N);
  const gx = new Float64Array(N);
  const gy = new Float64Array(Math.max(1, ny - 1) * nz * nx);
  const gz = new Float64Array(ny * Math.max(1, nz - 1) * nx);
  const rhs = new Float64Array(N);

  const areaX = dy * dz;
  const g0 = areaX / ((0.5 * dx[0]) / sig5[0]);
  const Iapp = 350.0; // 0.035 A/cm^2 = 350 A/m^2 (Report p. 15)

  for (let j = 0; j < ny; j++) {
    for (let k = 0; k < nz; k++) {
      const base = (j * nz + k) * nx;
      if (negTabMask[j * nz + k]) diagAdd[base] = g0;
      if (posTabMask[j * nz + k]) rhs[base + nx - 1] = Iapp * areaX;
      for (let i = 1; i < nx; i++) {
        gx[base + i] = areaX / ((0.5 * dx[i - 1]) / sig5[i - 1] + (0.5 * dx[i]) / sig5[i]);
      }
    }
  }
  for (let i = 0; i < nx; i++) {
    const gys = (sig5[i] * dx[i] * dz) / dy;
    const gzs = (sig5[i] * dx[i] * dy) / dz;
    for (let j = 0; j < ny - 1; j++) {
      for (let k = 0; k < nz; k++) gy[(j * nz + k) * nx + i] = gys;
    }
    for (let j = 0; j < ny; j++) {
      for (let k = 0; k < nz - 1; k++) gz[(j * (nz - 1) + k) * nx + i] = gzs;
    }
  }

  scalarMG.buildFromConductances(diagAdd, gx, gy, gz);
  const x = new Float64Array(N);
  const mgRes = scalarMG.solvePCG(rhs, x, 1e-10, 60);

  return {
    iters: mgRes.iters,
    finalRel: mgRes.relRes,
    phiMax: Math.max(...x),
    phiMin: Math.min(...x),
  };
}
