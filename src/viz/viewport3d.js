/**
 * Interactive Three.js 3D Stratified Battery Viewport.
 *
 * Renders the 5-region 3D cell using a single-draw-call InstancedMesh with:
 *   - Per-instance perceptual scientific colormaps (Turbo, Viridis, Inferno, RdBu)
 *   - Adjustable thickness exaggeration (1x to 250x, default 100x matching Report Figs. 9/21/26)
 *   - Exploded 5-layer inspection mode
 *   - Interactive X / Y / Z cutaway slicers exposing interior gradients
 *   - Metallic 3D tab terminals reflecting BC1-BC4 tab placement
 *   - 3D Current Density Vector Field arrows reconstructed from conservative hyperbolic face fluxes
 *   - Raycaster click/hover probe inspecting all 6 coupled field values at any voxel
 */

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { REGION_IDS, REGION_NAMES, REGION_SHORT } from '../physics/batteryParams.js';

const COLORMAP_STOPS = {
  turbo: [
    [0.18995, 0.07176, 0.23217],
    [0.25107, 0.25237, 0.63374],
    [0.19376, 0.49108, 0.90362],
    [0.13777, 0.72011, 0.80702],
    [0.32049, 0.89811, 0.49635],
    [0.68089, 0.95857, 0.22866],
    [0.95029, 0.78102, 0.14502],
    [0.98525, 0.45350, 0.10722],
    [0.83180, 0.17266, 0.05660],
    [0.47960, 0.01583, 0.01055],
  ],
  viridis: [
    [0.26700, 0.00487, 0.32942],
    [0.28262, 0.14092, 0.45751],
    [0.25394, 0.26525, 0.52998],
    [0.20676, 0.37175, 0.55312],
    [0.16362, 0.47113, 0.55815],
    [0.12757, 0.56695, 0.55056],
    [0.13469, 0.65864, 0.51765],
    [0.26694, 0.74875, 0.44057],
    [0.47750, 0.82144, 0.31820],
    [0.74139, 0.87345, 0.14956],
    [0.99325, 0.90616, 0.14394],
  ],
  inferno: [
    [0.00146, 0.00047, 0.01387],
    [0.08741, 0.04456, 0.22481],
    [0.25823, 0.03857, 0.40648],
    [0.41633, 0.09020, 0.43294],
    [0.57830, 0.14804, 0.40441],
    [0.73568, 0.21591, 0.33024],
    [0.86501, 0.31682, 0.22605],
    [0.95451, 0.46874, 0.11906],
    [0.98762, 0.64532, 0.03989],
    [0.96489, 0.84385, 0.27339],
    [0.98836, 0.99836, 0.64492],
  ],
  rdbu: [
    [0.08, 0.25, 0.72],
    [0.18, 0.52, 0.92],
    [0.56, 0.78, 0.96],
    [0.93, 0.95, 0.97],
    [0.97, 0.73, 0.56],
    [0.92, 0.36, 0.22],
    [0.68, 0.05, 0.08],
  ],
};

export function sampleColormap(name, t, outColor) {
  const stops = COLORMAP_STOPS[name] || COLORMAP_STOPS.turbo;
  const clamped = Math.max(0, Math.min(1, t));
  const scaled = clamped * (stops.length - 1);
  const idx = Math.min(stops.length - 2, Math.floor(scaled));
  const frac = scaled - idx;
  const a = stops[idx];
  const b = stops[idx + 1];
  outColor.setRGB(
    a[0] + frac * (b[0] - a[0]),
    a[1] + frac * (b[1] - a[1]),
    a[2] + frac * (b[2] - a[2]),
    THREE.SRGBColorSpace
  );
}

export class BatteryViewport3D {
  /**
   * @param {HTMLElement} container
   * @param {Object} callbacks
   */
  constructor(container, callbacks = {}) {
    this.container = container;
    this.callbacks = callbacks;

    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color('#090e1a');

    const w = Math.max(container.clientWidth, 300);
    const h = Math.max(container.clientHeight, 300);

    this.camera = new THREE.PerspectiveCamera(38, w / h, 0.05, 200);
    this.camera.position.set(5.6, 3.6, 6.2);

    this.renderer = new THREE.WebGLRenderer({ antialias: true, powerPreference: 'high-performance' });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.renderer.setSize(w, h);
    this.renderer.toneMapping = THREE.ACESFilmicToneMapping;
    this.renderer.toneMappingExposure = 1.08;
    container.appendChild(this.renderer.domElement);

    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.07;
    this.controls.target.set(0, 0, 0);

    // Lighting
    this.scene.add(new THREE.AmbientLight('#c7d2fe', 0.85));
    const keyLight = new THREE.DirectionalLight('#ffffff', 1.35);
    keyLight.position.set(7, 10, 9);
    this.scene.add(keyLight);
    const rimLight = new THREE.DirectionalLight('#38bdf8', 0.55);
    rimLight.position.set(-8, -4, -6);
    this.scene.add(rimLight);

    // Subtle ground grid
    this.gridHelper = new THREE.GridHelper(12, 24, '#1e293b', '#131c2e');
    this.gridHelper.position.y = -3.4;
    this.scene.add(this.gridHelper);

    // Groups
    this.batteryGroup = new THREE.Group();
    this.scene.add(this.batteryGroup);
    this.decorGroup = new THREE.Group();
    this.batteryGroup.add(this.decorGroup);

    // Meshes
    this.cellGeo = new THREE.BoxGeometry(1, 1, 1);
    this.cellMat = new THREE.MeshStandardMaterial({
      roughness: 0.36,
      metalness: 0.06,
      color: '#ffffff',
    });
    this.ghostMat = new THREE.MeshStandardMaterial({
      color: '#64748b',
      transparent: true,
      opacity: 0.14,
      roughness: 0.25,
      depthWrite: false,
    });
    this.cellMesh = null;
    this.ghostMesh = null;
    this.vectorMesh = null;
    this.probeMarker = null;

    // View state
    this.viewOpts = {
      field: 'phiE',
      colormap: 'turbo',
      exaggeration: 100,
      exploded: false,
      showVectors: true,
      sliceX: 1.0,
      sliceY: 1.0,
      sliceZ: 1.0,
    };

    this.snapshot = null;
    this.meshSignature = '';
    this.fieldRange = { min: 0, max: 1, unit: 'V', label: 'Electrolyte Potential' };
    this.probedCell = null;

    // Raycasting
    this.raycaster = new THREE.Raycaster();
    this.pointer = new THREE.Vector2();
    this._bindEvents();

    this._resizeObs = new ResizeObserver(() => this._handleResize());
    this._resizeObs.observe(container);

    this._animate();
  }

  _bindEvents() {
    const dom = this.renderer.domElement;
    let downPos = { x: 0, y: 0 };
    dom.addEventListener('pointerdown', (e) => {
      downPos = { x: e.clientX, y: e.clientY };
    });
    dom.addEventListener('pointerup', (e) => {
      const dist = Math.hypot(e.clientX - downPos.x, e.clientY - downPos.y);
      if (dist < 5) {
        this._pickCell(e);
      }
    });
  }

  _pickCell(e) {
    if (!this.cellMesh || !this.snapshot) return;
    const rect = this.renderer.domElement.getBoundingClientRect();
    this.pointer.x = ((e.clientX - rect.left) / rect.width) * 2 - 1;
    this.pointer.y = -((e.clientY - rect.top) / rect.height) * 2 + 1;
    this.raycaster.setFromCamera(this.pointer, this.camera);
    const hits = this.raycaster.intersectObject(this.cellMesh, false);
    if (hits.length > 0 && hits[0].instanceId !== undefined) {
      const c = hits[0].instanceId;
      this.probedCell = c;
      this._emitProbe(c);
    }
  }

  _emitProbe(c) {
    if (!this.snapshot || c === null || c === undefined) return;
    const { meshMeta, fields } = this.snapshot;
    const { nx, ny, nz, xCenters, Ly, Lz, reg } = meshMeta;
    const i = Math.floor(c / (ny * nz));
    const rem = c % (ny * nz);
    const j = Math.floor(rem / nz);
    const k = rem % nz;
    if (i < 0 || i >= nx) return;

    const r = reg[i];
    const isx = fields.isVec[3 * c];
    const isy = fields.isVec[3 * c + 1];
    const isz = fields.isVec[3 * c + 2];
    const iex = fields.ieVec[3 * c];
    const iey = fields.ieVec[3 * c + 1];
    const iez = fields.ieVec[3 * c + 2];

    this.callbacks.onProbe?.({
      c,
      i,
      j,
      k,
      regionId: r,
      regionName: REGION_NAMES[r],
      xUm: xCenters[i] * 1e6,
      yCm: ((j + 0.5) / ny) * Ly * 100,
      zCm: ((k + 0.5) / nz) * Lz * 100,
      phiS: r !== REGION_IDS.SEP ? fields.phiS[c] : null,
      phiE: r >= REGION_IDS.NEG && r <= REGION_IDS.POS ? fields.phiE[c] : null,
      jLi: r === REGION_IDS.NEG || r === REGION_IDS.POS ? fields.jLi[c] * 1e-6 : null,
      ce: r >= REGION_IDS.NEG && r <= REGION_IDS.POS ? fields.ce[c] : null,
      T: fields.T[c],
      stoich: r === REGION_IDS.NEG || r === REGION_IDS.POS ? fields.stoich[c] : null,
      iMagS: Math.hypot(isx, isy, isz) * 1e-4,
      iMagE: Math.hypot(iex, iey, iez) * 1e-4,
    });
  }

  setViewOptions(partial) {
    const needGeometry =
      (partial.exaggeration !== undefined && partial.exaggeration !== this.viewOpts.exaggeration) ||
      (partial.exploded !== undefined && partial.exploded !== this.viewOpts.exploded) ||
      (partial.sliceX !== undefined && partial.sliceX !== this.viewOpts.sliceX) ||
      (partial.sliceY !== undefined && partial.sliceY !== this.viewOpts.sliceY) ||
      (partial.sliceZ !== undefined && partial.sliceZ !== this.viewOpts.sliceZ) ||
      (partial.field !== undefined && partial.field !== this.viewOpts.field) ||
      (partial.showVectors !== undefined && partial.showVectors !== this.viewOpts.showVectors);

    Object.assign(this.viewOpts, partial);
    if (!this.snapshot) return;

    if (needGeometry) {
      this._updateTransforms();
      this._updateVectors();
    }
    this._updateColors();
  }

  updateSnapshot(snapshot) {
    this.snapshot = snapshot;
    const { meshMeta } = snapshot;
    const sig = `${meshMeta.nx}-${meshMeta.ny}-${meshMeta.nz}-${meshMeta.Ly}-${meshMeta.Lz}-${meshMeta.bc}`;
    if (sig !== this.meshSignature) {
      this.meshSignature = sig;
      this.probedCell = null;
      this._rebuildMeshes();
    }
    this._updateTransforms();
    this._updateColors();
    this._updateVectors();
    if (this.probedCell !== null && this.probedCell < meshMeta.N) {
      this._emitProbe(this.probedCell);
    }
  }

  _rebuildMeshes() {
    if (this.cellMesh) {
      this.batteryGroup.remove(this.cellMesh);
      this.cellMesh.dispose();
    }
    if (this.ghostMesh) {
      this.batteryGroup.remove(this.ghostMesh);
      this.ghostMesh.dispose();
    }
    if (this.vectorMesh) {
      this.batteryGroup.remove(this.vectorMesh);
      this.vectorMesh.dispose();
    }

    const N = this.snapshot.meshMeta.N;
    this.cellMesh = new THREE.InstancedMesh(this.cellGeo, this.cellMat, N);
    this.cellMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    const initColors = new Float32Array(N * 3);
    this.cellMesh.instanceColor = new THREE.InstancedBufferAttribute(initColors, 3);
    this.batteryGroup.add(this.cellMesh);

    this.ghostMesh = new THREE.InstancedMesh(this.cellGeo, this.ghostMat, N);
    this.ghostMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    this.batteryGroup.add(this.ghostMesh);

    // Arrow geometry for 3D current density vectors
    const arrowGeo = createArrowGeometry();
    const arrowMat = new THREE.MeshStandardMaterial({
      color: '#ffffff',
      roughness: 0.25,
      metalness: 0.15,
      emissive: '#1e293b',
      emissiveIntensity: 0.35,
    });
    const maxArrows = Math.min(600, N);
    this.vectorMesh = new THREE.InstancedMesh(arrowGeo, arrowMat, maxArrows);
    this.vectorMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    this.vectorMesh.instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(maxArrows * 3), 3);
    this.batteryGroup.add(this.vectorMesh);
  }

  _computeLayoutScales() {
    const { meshMeta } = this.snapshot;
    const { Lx, Ly, Lz } = meshMeta;
    const maxPlane = Math.max(Ly, Lz);
    const scaleYZ = 5.4 / maxPlane;
    const worldLy = Ly * scaleYZ;
    const worldLz = Lz * scaleYZ;
    const worldLx = Lx * scaleYZ * this.viewOpts.exaggeration;
    const gapWorld = this.viewOpts.exploded ? Math.max(0.26, worldLx * 0.16) : 0.0;
    const totalWorldX = worldLx + 4 * gapWorld;
    return { scaleYZ, worldLx, worldLy, worldLz, gapWorld, totalWorldX };
  }

  _isFieldActiveInRegion(field, r) {
    if (field === 'phiS') return r !== REGION_IDS.SEP;
    if (field === 'phiE' || field === 'ce') return r >= REGION_IDS.NEG && r <= REGION_IDS.POS;
    if (field === 'jLi' || field === 'stoich') return r === REGION_IDS.NEG || r === REGION_IDS.POS;
    return true; // T, qTotal
  }

  _updateTransforms() {
    if (!this.snapshot || !this.cellMesh) return;
    const { meshMeta } = this.snapshot;
    const { nx, ny, nz, xCenters, dx, reg } = meshMeta;
    const { scaleYZ, worldLy, worldLz, gapWorld, totalWorldX } = this._computeLayoutScales();
    const { sliceX, sliceY, sliceZ, field } = this.viewOpts;

    const maxI = Math.max(0, Math.min(nx - 1, Math.ceil(sliceX * nx) - 1));
    const maxJ = Math.max(0, Math.min(ny - 1, Math.ceil(sliceY * ny) - 1));
    const maxK = Math.max(0, Math.min(nz - 1, Math.ceil(sliceZ * nz) - 1));

    const worldDy = worldLy / ny;
    const worldDz = worldLz / nz;
    const dummy = new THREE.Object3D();

    for (let i = 0; i < nx; i++) {
      const r = reg[i];
      const worldDx = dx[i] * scaleYZ * this.viewOpts.exaggeration;
      const wx = xCenters[i] * scaleYZ * this.viewOpts.exaggeration + r * gapWorld - 0.5 * totalWorldX;
      const active = this._isFieldActiveInRegion(field, r);

      for (let j = 0; j < ny; j++) {
        const wy = (j + 0.5) * worldDy - 0.5 * worldLy;
        for (let k = 0; k < nz; k++) {
          const wz = (k + 0.5) * worldDz - 0.5 * worldLz;
          const c = (i * ny + j) * nz + k;
          const inSlice = i <= maxI && j <= maxJ && k <= maxK;

          if (!inSlice) {
            dummy.scale.set(0, 0, 0);
            dummy.position.set(wx, wy, wz);
            dummy.updateMatrix();
            this.cellMesh.setMatrixAt(c, dummy.matrix);
            this.ghostMesh.setMatrixAt(c, dummy.matrix);
          } else if (active) {
            dummy.position.set(wx, wy, wz);
            dummy.scale.set(worldDx * 0.985, worldDy * 0.985, worldDz * 0.985);
            dummy.updateMatrix();
            this.cellMesh.setMatrixAt(c, dummy.matrix);
            dummy.scale.set(0, 0, 0);
            dummy.updateMatrix();
            this.ghostMesh.setMatrixAt(c, dummy.matrix);
          } else {
            dummy.position.set(wx, wy, wz);
            dummy.scale.set(worldDx * 0.97, worldDy * 0.97, worldDz * 0.97);
            dummy.updateMatrix();
            this.ghostMesh.setMatrixAt(c, dummy.matrix);
            dummy.scale.set(0, 0, 0);
            dummy.updateMatrix();
            this.cellMesh.setMatrixAt(c, dummy.matrix);
          }
        }
      }
    }

    this.cellMesh.instanceMatrix.needsUpdate = true;
    this.ghostMesh.instanceMatrix.needsUpdate = true;
    this.cellMesh.computeBoundingSphere();

    this._rebuildDecor(scaleYZ, worldLy, worldLz, gapWorld, totalWorldX);
  }

  _rebuildDecor(scaleYZ, worldLy, worldLz, gapWorld, totalWorldX) {
    while (this.decorGroup.children.length > 0) {
      const obj = this.decorGroup.children.pop();
      obj.geometry?.dispose?.();
      obj.material?.dispose?.();
    }

    const { meshMeta } = this.snapshot;
    const { nx, ny, nz, dx, xEdges, reg, negTabMask, posTabMask, negTabType } = meshMeta;

    // Region boundary outline wireframes
    const regionColors = ['#f59e0b', '#38bdf8', '#94a3b8', '#f43f5e', '#cbd5e1'];
    let xAcc = 0;
    for (let r = 0; r < 5; r++) {
      let rThick = 0;
      for (let i = 0; i < nx; i++) {
        if (reg[i] === r) rThick += dx[i];
      }
      const wRx = rThick * scaleYZ * this.viewOpts.exaggeration;
      const cx = (xAcc + 0.5 * rThick) * scaleYZ * this.viewOpts.exaggeration + r * gapWorld - 0.5 * totalWorldX;
      xAcc += rThick;

      const box = new THREE.BoxGeometry(wRx, worldLy, worldLz);
      const edges = new THREE.EdgesGeometry(box);
      const line = new THREE.LineSegments(
        edges,
        new THREE.LineBasicMaterial({ color: regionColors[r], transparent: true, opacity: 0.42 })
      );
      line.position.set(cx, 0, 0);
      this.decorGroup.add(line);
      box.dispose();
    }

    // Metallic 3D Tabs for Negative (Cu) and Positive (Al)
    const worldDy = worldLy / ny;
    const worldDz = worldLz / nz;
    const tabLen = 0.42;

    const negTabMat = new THREE.MeshStandardMaterial({
      color: '#d97706',
      metalness: 0.78,
      roughness: 0.25,
      emissive: '#78350f',
      emissiveIntensity: 0.25,
    });
    const posTabMat = new THREE.MeshStandardMaterial({
      color: '#e2e8f0',
      metalness: 0.85,
      roughness: 0.22,
      emissive: '#0284c7',
      emissiveIntensity: 0.2,
    });

    // Compute bounding boxes of Neg and Pos tabs
    let jMinN = ny, jMaxN = -1, kMinN = nz, kMaxN = -1;
    let jMinP = ny, jMaxP = -1, kMinP = nz, kMaxP = -1;
    for (let j = 0; j < ny; j++) {
      for (let k = 0; k < nz; k++) {
        if (negTabMask[j * nz + k]) {
          if (j < jMinN) jMinN = j;
          if (j > jMaxN) jMaxN = j;
          if (k < kMinN) kMinN = k;
          if (k > kMaxN) kMaxN = k;
        }
        if (posTabMask[j * nz + k]) {
          if (j < jMinP) jMinP = j;
          if (j > jMaxP) jMaxP = j;
          if (k < kMinP) kMinP = k;
          if (k > kMaxP) kMaxP = k;
        }
      }
    }

    if (jMaxN >= 0) {
      const cy = ((jMinN + jMaxN + 1) * 0.5) * worldDy - 0.5 * worldLy;
      const cz = ((kMinN + kMaxN + 1) * 0.5) * worldDz - 0.5 * worldLz;
      const hY = (jMaxN - jMinN + 1) * worldDy;
      const wZ = (kMaxN - kMinN + 1) * worldDz;

      if (negTabType === 'top_edge') {
        const wNcc = meshMeta.nxCounts[0] > 0 ? xEdges[meshMeta.nxCounts[0]] * scaleYZ * this.viewOpts.exaggeration : 0.08;
        const tabGeo = new THREE.BoxGeometry(wNcc, tabLen, wZ);
        const tab = new THREE.Mesh(tabGeo, negTabMat);
        tab.position.set(-0.5 * totalWorldX + 0.5 * wNcc, 0.5 * worldLy + 0.5 * tabLen, cz);
        this.decorGroup.add(tab);
      } else {
        const tabGeo = new THREE.BoxGeometry(tabLen, hY, wZ);
        const tab = new THREE.Mesh(tabGeo, negTabMat);
        tab.position.set(-0.5 * totalWorldX - 0.5 * tabLen, cy, cz);
        this.decorGroup.add(tab);
      }
    }

    if (jMaxP >= 0) {
      const cy = ((jMinP + jMaxP + 1) * 0.5) * worldDy - 0.5 * worldLy;
      const cz = ((kMinP + kMaxP + 1) * 0.5) * worldDz - 0.5 * worldLz;
      const hY = (jMaxP - jMinP + 1) * worldDy;
      const wZ = (kMaxP - kMinP + 1) * worldDz;
      const tabGeo = new THREE.BoxGeometry(tabLen, hY, wZ);
      const tab = new THREE.Mesh(tabGeo, posTabMat);
      tab.position.set(0.5 * totalWorldX + 0.5 * tabLen, cy, cz);
      this.decorGroup.add(tab);
    }

    this.gridHelper.position.y = -0.5 * worldLy - 0.4;
  }

  _updateColors() {
    if (!this.snapshot || !this.cellMesh) return;
    const { meshMeta, fields } = this.snapshot;
    const { nx, ny, nz, reg } = meshMeta;
    const field = this.viewOpts.field;

    const metaMap = {
      phiS: { arr: fields.phiS, unit: 'V', label: 'Solid Potential φ_s', defaultMap: 'turbo', sym: false },
      phiE: { arr: fields.phiE, unit: 'V', label: 'Electrolyte Potential φ_e', defaultMap: 'viridis', sym: false },
      jLi: { arr: fields.jLi, scale: 1e-6, unit: 'A/cm³', label: 'Reaction Current j^Li', defaultMap: 'rdbu', sym: true },
      ce: { arr: fields.ce, unit: 'mol/m³', label: 'Electrolyte Conc. c_e', defaultMap: 'turbo', sym: false },
      T: { arr: fields.T, unit: 'K', label: 'Temperature T', defaultMap: 'inferno', sym: false },
      stoich: { arr: fields.stoich, unit: 'θ', label: 'Solid Stoichiometry θ_s', defaultMap: 'viridis', sym: false },
    };

    const meta = metaMap[field] || metaMap.phiE;
    const arr = meta.arr;
    const scale = meta.scale || 1.0;

    let minVal = Infinity;
    let maxVal = -Infinity;
    for (let i = 0; i < nx; i++) {
      if (!this._isFieldActiveInRegion(field, reg[i])) continue;
      for (let j = 0; j < ny; j++) {
        for (let k = 0; k < nz; k++) {
          const v = arr[(i * ny + j) * nz + k] * scale;
          if (v < minVal) minVal = v;
          if (v > maxVal) maxVal = v;
        }
      }
    }

    if (!Number.isFinite(minVal)) {
      minVal = 0;
      maxVal = 1;
    }
    if (meta.sym) {
      const m = Math.max(Math.abs(minVal), Math.abs(maxVal), 1e-6);
      minVal = -m;
      maxVal = m;
    } else if (Math.abs(maxVal - minVal) < 1e-12) {
      minVal -= 1e-4;
      maxVal += 1e-4;
    }

    const span = maxVal - minVal;
    const cmap = this.viewOpts.colormap === 'auto' ? meta.defaultMap : this.viewOpts.colormap;
    const tmpColor = new THREE.Color();

    for (let i = 0; i < nx; i++) {
      const active = this._isFieldActiveInRegion(field, reg[i]);
      for (let j = 0; j < ny; j++) {
        for (let k = 0; k < nz; k++) {
          const c = (i * ny + j) * nz + k;
          if (!active) {
            tmpColor.setRGB(0.2, 0.25, 0.35);
          } else {
            const v = arr[c] * scale;
            sampleColormap(cmap, (v - minVal) / span, tmpColor);
          }
          this.cellMesh.setColorAt(c, tmpColor);
        }
      }
    }
    this.cellMesh.instanceColor.needsUpdate = true;

    this.fieldRange = {
      min: minVal,
      max: maxVal,
      unit: meta.unit,
      label: meta.label,
      cmap,
    };
    this.callbacks.onColorbarChange?.(this.fieldRange);
  }

  _updateVectors() {
    if (!this.snapshot || !this.vectorMesh) return;
    if (!this.viewOpts.showVectors) {
      this.vectorMesh.count = 0;
      return;
    }

    const { meshMeta, fields } = this.snapshot;
    const { nx, ny, nz, xCenters, reg } = meshMeta;
    const { scaleYZ, worldLy, worldLz, gapWorld, totalWorldX } = this._computeLayoutScales();
    const { sliceX, sliceY, sliceZ, field } = this.viewOpts;

    const maxI = Math.max(0, Math.min(nx - 1, Math.ceil(sliceX * nx) - 1));
    const maxJ = Math.max(0, Math.min(ny - 1, Math.ceil(sliceY * ny) - 1));
    const maxK = Math.max(0, Math.min(nz - 1, Math.ceil(sliceZ * nz) - 1));

    const worldDy = worldLy / ny;
    const worldDz = worldLz / nz;

    // Sample cleanly at up to ~6 x 7 x 5 grid points
    const stepI = Math.max(1, Math.floor(nx / 7));
    const stepJ = Math.max(1, Math.floor(ny / 6));
    const stepK = Math.max(1, Math.floor(nz / 5));

    const useElec = field === 'phiE' || field === 'ce';
    const vecArr = useElec ? fields.ieVec : fields.isVec;

    // Find max vector magnitude for clean scaling
    let maxMag = 0;
    const candidates = [];
    for (let i = 0; i <= maxI; i += stepI) {
      const r = reg[i];
      if (useElec && (r === REGION_IDS.NEG_CC || r === REGION_IDS.POS_CC)) continue;
      if (!useElec && r === REGION_IDS.SEP) continue;
      for (let j = Math.floor(stepJ / 2); j <= maxJ; j += stepJ) {
        for (let k = Math.floor(stepK / 2); k <= maxK; k += stepK) {
          const c = (i * ny + j) * nz + k;
          // Note: scale x-component visually with exaggeration so arrows point realistically in world coordinates
          const vx = vecArr[3 * c];
          const vy = vecArr[3 * c + 1];
          const vz = vecArr[3 * c + 2];
          const mag = Math.hypot(vx, vy, vz);
          if (mag > maxMag) maxMag = mag;
          if (mag > 1e-8) {
            candidates.push({ i, j, k, c, vx, vy, vz, mag });
          }
        }
      }
    }

    const dummy = new THREE.Object3D();
    const up = new THREE.Vector3(0, 1, 0);
    const dirVec = new THREE.Vector3();
    const color = new THREE.Color();
    const maxArrows = this.vectorMesh.instanceMatrix.count;
    let count = 0;

    for (const item of candidates) {
      if (count >= maxArrows) break;
      const { i, j, k, vx, vy, vz, mag } = item;
      const r = reg[i];
      const wx = xCenters[i] * scaleYZ * this.viewOpts.exaggeration + r * gapWorld - 0.5 * totalWorldX;
      const wy = (j + 0.5) * worldDy - 0.5 * worldLy;
      const wz = (k + 0.5) * worldDz - 0.5 * worldLz;

      // In world space, x is exaggerated, y & z are physical aspect ratio
      dirVec.set(vx * (this.viewOpts.exaggeration * 0.08), vy, vz);
      if (dirVec.lengthSq() < 1e-18) continue;
      dirVec.normalize();

      const normMag = Math.log1p((mag / Math.max(maxMag, 1e-12)) * 12.0) / Math.log1p(12.0);
      const len = 0.22 + 0.42 * normMag;

      dummy.position.set(wx, wy, wz);
      dummy.quaternion.setFromUnitVectors(up, dirVec);
      dummy.scale.set(len, len, len);
      dummy.updateMatrix();
      this.vectorMesh.setMatrixAt(count, dummy.matrix);

      color.setHSL(0.55 - 0.45 * normMag, 0.95, 0.62);
      this.vectorMesh.setColorAt(count, color);
      count++;
    }

    this.vectorMesh.count = count;
    this.vectorMesh.instanceMatrix.needsUpdate = true;
    if (this.vectorMesh.instanceColor) {
      this.vectorMesh.instanceColor.needsUpdate = true;
    }
  }

  resetCamera() {
    this.camera.position.set(5.6, 3.6, 6.2);
    this.controls.target.set(0, 0, 0);
    this.controls.update();
  }

  _handleResize() {
    const w = this.container.clientWidth;
    const h = this.container.clientHeight;
    if (w < 10 || h < 10) return;
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this.renderer.setSize(w, h);
  }

  _animate() {
    this._raf = requestAnimationFrame(() => this._animate());
    this.controls.update();
    this.renderer.render(this.scene, this.camera);
  }

  dispose() {
    cancelAnimationFrame(this._raf);
    this._resizeObs.disconnect();
    this.renderer.dispose();
  }
}

/**
 * Create a clean merged arrow BufferGeometry pointing along +Y (length 1.0).
 */
function createArrowGeometry() {
  const shaft = new THREE.CylinderGeometry(0.032, 0.032, 0.65, 8);
  shaft.translate(0, 0.05, 0);
  const head = new THREE.ConeGeometry(0.085, 0.35, 10);
  head.translate(0, 0.55, 0);

  const g1 = shaft.toNonIndexed();
  const g2 = head.toNonIndexed();
  const pos1 = g1.attributes.position.array;
  const pos2 = g2.attributes.position.array;
  const nrm1 = g1.attributes.normal.array;
  const nrm2 = g2.attributes.normal.array;

  const pos = new Float32Array(pos1.length + pos2.length);
  pos.set(pos1, 0);
  pos.set(pos2, pos1.length);

  const nrm = new Float32Array(nrm1.length + nrm2.length);
  nrm.set(nrm1, 0);
  nrm.set(nrm2, nrm1.length);

  const merged = new THREE.BufferGeometry();
  merged.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  merged.setAttribute('normal', new THREE.BufferAttribute(nrm, 3));
  shaft.dispose();
  head.dispose();
  g1.dispose();
  g2.dispose();
  return merged;
}
