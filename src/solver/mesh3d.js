/**
 * Structured 5-Region 3D Finite-Volume Mesh & Active-Region DOF Mapping.
 *
 * Regions along x (thickness):
 *   0: Negative Current Collector (Cu)   - phi_s, T
 *   1: Negative Electrode (Li_x C_6)     - phi_s, phi_e, j^Li, c_s, c_e, T
 *   2: Separator (Polymer)               - phi_e, c_e, T
 *   3: Positive Electrode (Li_y Mn_2 O_4)- phi_s, phi_e, j^Li, c_s, c_e, T
 *   4: Positive Current Collector (Al)   - phi_s, T
 *
 * Crucially, NO dummy DOFs are allocated (no phi_s in separator, no phi_e in CCs),
 * and along every (j, k) column the active potential DOFs are stored contiguously
 * so that the coupled through-thickness X-line solve is strictly banded with
 * half-bandwidth 2 (pentadiagonal SPD).
 */

import { LAYER_SPECS, REGION_IDS } from '../physics/batteryParams.js';

export class BatteryMesh3D {
  /**
   * @param {Object} opts
   * @param {number} opts.Ly - Height along y [m]
   * @param {number} opts.Lz - Width along z [m]
   * @param {number[]} opts.nxCounts - [nNegCC, nNeg, nSep, nPos, nPosCC]
   * @param {number} opts.ny - Cells along height y
   * @param {number} opts.nz - Cells along width z
   * @param {string} opts.bc - 'BC1' | 'BC2' | 'BC3' | 'BC4'
   * @param {boolean} [opts.clustered=false] - Interface-clustered x spacing (Report Fig. 22)
   */
  constructor({ Ly, Lz, nxCounts, ny, nz, bc = 'BC1', clustered = false }) {
    this.Ly = Ly;
    this.Lz = Lz;
    this.nxCounts = [...nxCounts];
    this.ny = ny;
    this.nz = nz;
    this.bc = bc;
    this.clustered = clustered;
    this.dy = Ly / ny;
    this.dz = Lz / nz;

    const dxArr = [];
    const regArr = [];
    const sigArr = [];
    const epsEArr = [];
    const epsSArr = [];
    const rsArr = [];
    const asArr = [];
    const rhoCpArr = [];
    const lamArr = [];

    for (let r = 0; r < 5; r++) {
      const spec = LAYER_SPECS[r];
      const nc = nxCounts[r];
      const widths = computeRegionWidths(spec.thickness, nc, r, clustered);
      for (let m = 0; m < nc; m++) {
        dxArr.push(widths[m]);
        regArr.push(spec.id);
        sigArr.push(spec.sigma);
        epsEArr.push(spec.epsE);
        epsSArr.push(spec.epsS);
        rsArr.push(spec.Rs);
        asArr.push(spec.asArea);
        rhoCpArr.push(spec.rho * spec.cp);
        lamArr.push(spec.lambda);
      }
    }

    this.nx = dxArr.length;
    this.N = this.nx * ny * nz;
    this.dx = new Float64Array(dxArr);
    this.reg = new Int32Array(regArr);
    this.sigma = new Float64Array(sigArr);
    this.epsE = new Float64Array(epsEArr);
    this.epsS = new Float64Array(epsSArr);
    this.Rs = new Float64Array(rsArr);
    this.asArea = new Float64Array(asArr);
    this.rhoCp = new Float64Array(rhoCpArr);
    this.lambda = new Float64Array(lamArr);

    // Physical coordinates along thickness x
    this.xEdges = new Float64Array(this.nx + 1);
    this.xCenters = new Float64Array(this.nx);
    let xAcc = 0.0;
    for (let i = 0; i < this.nx; i++) {
      this.xEdges[i] = xAcc;
      this.xCenters[i] = xAcc + 0.5 * this.dx[i];
      xAcc += this.dx[i];
    }
    this.xEdges[this.nx] = xAcc;
    this.Lx = xAcc;

    // Region interface x coordinates
    this.regionBounds = [];
    let curX = 0.0;
    for (let r = 0; r < 5; r++) {
      const w = LAYER_SPECS[r].thickness;
      this.regionBounds.push({ start: curX, end: curX + w, width: w });
      curX += w;
    }

    // Active DOF mappings along x
    // 1. Coupled (phi_s, phi_e) column layout: bandwidth <= 2!
    this.sLoc = new Int32Array(this.nx).fill(-1);
    this.eLoc = new Int32Array(this.nx).fill(-1);
    let colN = 0;
    for (let i = 0; i < this.nx; i++) {
      if (this.reg[i] !== REGION_IDS.SEP) {
        this.sLoc[i] = colN++;
      }
      if (
        this.reg[i] === REGION_IDS.NEG ||
        this.reg[i] === REGION_IDS.SEP ||
        this.reg[i] === REGION_IDS.POS
      ) {
        this.eLoc[i] = colN++;
      }
    }
    this.colN = colN;
    this.Npot = ny * nz * colN;

    // 2. Electrolyte-only contiguous x indices (for scalar c_e implicit solve)
    const eIdxList = [];
    this.eOnlyLoc = new Int32Array(this.nx).fill(-1);
    for (let i = 0; i < this.nx; i++) {
      if (this.eLoc[i] >= 0) {
        this.eOnlyLoc[i] = eIdxList.length;
        eIdxList.push(i);
      }
    }
    this.eXIndices = new Int32Array(eIdxList);
    this.nxElec = eIdxList.length;
    this.Ne = this.nxElec * ny * nz;

    // Tab masks on negative (x = 0 or top edge of NegCC) & positive (x = Lx) collectors
    this.negTabMask = new Uint8Array(ny * nz); // face mask on x=0 (or y=Ly for BC3)
    this.posTabMask = new Uint8Array(ny * nz); // face mask on x=Lx
    this.negTabType = bc === 'BC3' ? 'top_edge' : 'left_face';
    this.negTabCells = [];
    this.posTabCells = [];
    this._buildTabMasks(bc);

    // Precompute cell volumes and active face list for Mixed/Hyperbolic & current reconstruction
    this.cellVol = new Float64Array(this.nx);
    for (let i = 0; i < this.nx; i++) {
      this.cellVol[i] = this.dx[i] * this.dy * this.dz;
    }

    this._buildFacesTopology();
  }

  _buildTabMasks(bc) {
    const { ny, nz } = this;
    const tabRows = Math.max(1, Math.round(ny / 10));
    const kHalf = Math.max(1, Math.floor(nz / 2));

    if (bc === 'BC1') {
      // Opposite corners: Neg tab bottom-front (y near 0, z < Lz/2), Pos tab top-back (y near Ly, z >= Lz/2)
      for (let j = 0; j < tabRows; j++) {
        for (let k = 0; k < kHalf; k++) {
          this.negTabMask[j * nz + k] = 1;
          this.negTabCells.push([j, k]);
        }
      }
      for (let j = ny - tabRows; j < ny; j++) {
        for (let k = kHalf; k < nz; k++) {
          this.posTabMask[j * nz + k] = 1;
          this.posTabCells.push([j, k]);
        }
      }
    } else if (bc === 'BC2') {
      // Same top edge: Neg tab (y near Ly, z < Lz/2), Pos tab (y near Ly, z >= Lz/2)
      for (let j = ny - tabRows; j < ny; j++) {
        for (let k = 0; k < kHalf; k++) {
          this.negTabMask[j * nz + k] = 1;
          this.negTabCells.push([j, k]);
        }
        for (let k = kHalf; k < nz; k++) {
          this.posTabMask[j * nz + k] = 1;
          this.posTabCells.push([j, k]);
        }
      }
    } else if (bc === 'BC3') {
      // BC3 (Report Fig. 13): Neg tab on top surface (y=Ly) of NegCC; Pos tab on top of right face
      // Map negTabMask onto (i_ncc, k) via negTabEdgeCells, and keep negTabCells as (j=ny-1, k)
      for (let k = 0; k < kHalf; k++) {
        this.negTabMask[(ny - 1) * nz + k] = 1;
        this.negTabCells.push([ny - 1, k]);
      }
      for (let j = ny - tabRows; j < ny; j++) {
        for (let k = kHalf; k < nz; k++) {
          this.posTabMask[j * nz + k] = 1;
          this.posTabCells.push([j, k]);
        }
      }
    } else {
      // BC4: Full faces (1D equivalence check, Report Fig. 10)
      for (let j = 0; j < ny; j++) {
        for (let k = 0; k < nz; k++) {
          this.negTabMask[j * nz + k] = 1;
          this.posTabMask[j * nz + k] = 1;
          this.negTabCells.push([j, k]);
          this.posTabCells.push([j, k]);
        }
      }
    }

    this.negTabArea =
      bc === 'BC3'
        ? this.nxCounts[0] > 0
          ? LAYER_SPECS[0].thickness * (kHalf * this.dz)
          : this.dy * this.dz
        : this.negTabCells.length * this.dy * this.dz;
    this.posTabArea = this.posTabCells.length * this.dy * this.dz;
  }

  _buildFacesTopology() {
    // Store static face arrays as flat typed arrays for maximum V8 throughput
    const { nx, ny, nz, sLoc, eLoc } = this;
    const cL = [];
    const cR = [];
    const dir = [];      // 0: x, 1: y, 2: z
    const isSolid = [];  // 1: solid phase, 0: electrolyte phase
    const cellL = [];    // primary 3D cell index (-1 if boundary)
    const cellR = [];    // primary 3D cell index
    const iArr = [];

    for (let j = 0; j < ny; j++) {
      for (let k = 0; k < nz; k++) {
        // Negative Dirichlet tab face
        if (this.negTabType === 'left_face' && this.negTabMask[j * nz + k]) {
          cL.push(-1);
          cR.push(this.potIdx(j, k, sLoc[0]));
          dir.push(0);
          isSolid.push(1);
          cellL.push(-1);
          cellR.push(this.cid(0, j, k));
          iArr.push(0);
        } else if (this.negTabType === 'top_edge' && j === ny - 1 && this.negTabMask[j * nz + k]) {
          // Grounded top y-face of NegCC cells
          for (let i = 0; i < this.nxCounts[0]; i++) {
            cL.push(-2); // special flag: y-top Dirichlet on cell i
            cR.push(this.potIdx(j, k, sLoc[i]));
            dir.push(1);
            isSolid.push(1);
            cellL.push(-2);
            cellR.push(this.cid(i, j, k));
            iArr.push(i);
          }
        }

        // Interior x-faces
        for (let i = 0; i < nx - 1; i++) {
          if (sLoc[i] >= 0 && sLoc[i + 1] >= 0) {
            cL.push(this.potIdx(j, k, sLoc[i]));
            cR.push(this.potIdx(j, k, sLoc[i + 1]));
            dir.push(0);
            isSolid.push(1);
            cellL.push(this.cid(i, j, k));
            cellR.push(this.cid(i + 1, j, k));
            iArr.push(i);
          }
          if (eLoc[i] >= 0 && eLoc[i + 1] >= 0) {
            cL.push(this.potIdx(j, k, eLoc[i]));
            cR.push(this.potIdx(j, k, eLoc[i + 1]));
            dir.push(0);
            isSolid.push(0);
            cellL.push(this.cid(i, j, k));
            cellR.push(this.cid(i + 1, j, k));
            iArr.push(i);
          }
        }

        // Interior y & z faces
        for (let i = 0; i < nx; i++) {
          if (j < ny - 1) {
            if (sLoc[i] >= 0) {
              cL.push(this.potIdx(j, k, sLoc[i]));
              cR.push(this.potIdx(j + 1, k, sLoc[i]));
              dir.push(1);
              isSolid.push(1);
              cellL.push(this.cid(i, j, k));
              cellR.push(this.cid(i, j + 1, k));
              iArr.push(i);
            }
            if (eLoc[i] >= 0) {
              cL.push(this.potIdx(j, k, eLoc[i]));
              cR.push(this.potIdx(j + 1, k, eLoc[i]));
              dir.push(1);
              isSolid.push(0);
              cellL.push(this.cid(i, j, k));
              cellR.push(this.cid(i, j + 1, k));
              iArr.push(i);
            }
          }
          if (k < nz - 1) {
            if (sLoc[i] >= 0) {
              cL.push(this.potIdx(j, k, sLoc[i]));
              cR.push(this.potIdx(j, k + 1, sLoc[i]));
              dir.push(2);
              isSolid.push(1);
              cellL.push(this.cid(i, j, k));
              cellR.push(this.cid(i, j, k + 1));
              iArr.push(i);
            }
            if (eLoc[i] >= 0) {
              cL.push(this.potIdx(j, k, eLoc[i]));
              cR.push(this.potIdx(j, k + 1, eLoc[i]));
              dir.push(2);
              isSolid.push(0);
              cellL.push(this.cid(i, j, k));
              cellR.push(this.cid(i, j, k + 1));
              iArr.push(i);
            }
          }
        }
      }
    }

    this.nFaces = cL.length;
    this.faceCL = new Int32Array(cL);
    this.faceCR = new Int32Array(cR);
    this.faceDir = new Uint8Array(dir);
    this.faceIsSolid = new Uint8Array(isSolid);
    this.faceCellL = new Int32Array(cellL);
    this.faceCellR = new Int32Array(cellR);
    this.faceI = new Int32Array(iArr);
  }

  cid(i, j, k) {
    return (i * this.ny + j) * this.nz + k;
  }

  potIdx(j, k, m) {
    return (j * this.nz + k) * this.colN + m;
  }

  elecIdx(ie, j, k) {
    return (j * this.nz + k) * this.nxElec + ie;
  }

  thermIdx(i, j, k) {
    return (j * this.nz + k) * this.nx + i;
  }
}

/**
 * Generates uniform or interface-clustered cell widths within a region.
 */
function computeRegionWidths(thickness, nc, regionId, clustered) {
  if (!clustered || nc <= 2) {
    return new Array(nc).fill(thickness / nc);
  }
  // Smooth geometric clustering toward interfaces (Report Fig. 22)
  const raw = new Array(nc);
  const ratio = 1.25;
  for (let m = 0; m < nc; m++) {
    if (regionId === 0) {
      // Cluster toward right interface (x = Lca)
      raw[m] = Math.pow(ratio, nc - 1 - m);
    } else if (regionId === 4) {
      // Cluster toward left interface (x = Lcc)
      raw[m] = Math.pow(ratio, m);
    } else {
      // Cluster toward both left & right interfaces
      const distEdge = Math.min(m, nc - 1 - m);
      raw[m] = Math.pow(ratio, distEdge);
    }
  }
  const sum = raw.reduce((a, b) => a + b, 0);
  return raw.map((v) => (v / sum) * thickness);
}
