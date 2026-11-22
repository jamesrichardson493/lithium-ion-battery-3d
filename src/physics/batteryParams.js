/**
 * Physical, Electrochemical, Kinetic, Transport, and Thermal Parameters
 * for the 5-Region Li_x C_6 | Li_y Mn_2 O_4 Lithium-Ion Cell.
 *
 * All internal values are strictly in SI units (m, s, A, V, mol, K, J, W)
 * derived directly from Doyle (1995/1996), Dualfoil, Gu & Wang (2000),
 * Srinivasan & Wang (2003), and Anderson & Kapadia (2013) Tables 5 & 14.
 */

export const F_CONST = 96487.0;   // Faraday constant [C/mol]
export const R_GAS = 8.3143;      // Universal gas constant [J/(mol*K)]
export const T_AMBIENT = 298.15;  // Reference ambient temperature [K] (25 C)
export const V_CUTOFF = 2.80;     // Discharge cutoff voltage [V]
export const C_RATE_1C = 17.5;    // 1C superficial current density [A/m^2] (= 0.00175 A/cm^2, Report p. 28)

export const REGION_IDS = {
  NEG_CC: 0,
  NEG: 1,
  SEP: 2,
  POS: 3,
  POS_CC: 4,
};

export const REGION_NAMES = [
  'Neg Current Collector (Cu)',
  'Negative Electrode (Li_x C_6)',
  'Separator (Polymer)',
  'Positive Electrode (Li_y Mn_2 O_4)',
  'Pos Current Collector (Al)',
];

export const REGION_SHORT = ['NegCC', 'Neg', 'Sep', 'Pos', 'PosCC'];

/**
 * Layer physical specifications (Report Table 5 & Table 14).
 * Total thickness = 346 um = 0.0346 cm.
 */
export const LAYER_SPECS = [
  {
    id: REGION_IDS.NEG_CC,
    name: 'Neg Current Collector (Cu)',
    shortName: 'NegCC',
    thickness: 10e-6,        // 0.001 cm
    sigma: 1.0e7,            // 1.0e5 S/cm (Table 5)
    epsE: 0.0,
    epsS: 0.0,
    Rs: 0.0,
    asArea: 0.0,
    rho: 8933.0,             // Copper density [kg/m^3]
    cp: 385.0,               // Copper specific heat [J/(kg*K)]
    lambda: 398.0,           // Copper thermal conductivity [W/(m*K)]
  },
  {
    id: REGION_IDS.NEG,
    name: 'Negative Electrode (Li_x C_6)',
    shortName: 'Neg',
    thickness: 100e-6,       // 0.01 cm
    sigma: 100.0,            // 1.0 S/cm (Table 5)
    epsE: 0.357,             // Doyle (1996) porosity
    epsS: 0.471,             // Active material volume fraction
    Rs: 12.5e-6,             // Particle radius [m]
    asArea: (3.0 * 0.471) / 12.5e-6, // Specific interfacial area [1/m]
    csMax: 26390.0,          // Max solid Li concentration [mol/m^3]
    theta0: 0.63,            // Initial stoichiometry
    DsRef: 3.9e-14,          // Solid diffusivity at 298.15 K [m^2/s]
    kRef: 2.0e-11,           // Reaction rate constant [m^2.5/(mol^0.5*s)]
    EaD: 30000.0,            // Solid diffusion activation energy [J/mol]
    EaK: 20000.0,            // Reaction activation energy [J/mol]
    rho: 2500.0,
    cp: 700.0,
    lambda: 5.0,
  },
  {
    id: REGION_IDS.SEP,
    name: 'Separator (Polymer)',
    shortName: 'Sep',
    thickness: 52e-6,        // 0.0052 cm
    sigma: 0.0,              // Electronic insulator
    epsE: 0.724,             // Doyle (1996) separator porosity
    epsS: 0.0,
    Rs: 0.0,
    asArea: 0.0,
    rho: 1200.0,
    cp: 700.0,
    lambda: 1.0,
  },
  {
    id: REGION_IDS.POS,
    name: 'Positive Electrode (Li_y Mn_2 O_4)',
    shortName: 'Pos',
    thickness: 174e-6,       // 0.0174 cm
    sigma: 3.8,              // 0.038 S/cm (Table 5)
    epsE: 0.444,             // Doyle (1996) porosity
    epsS: 0.297,             // Active material volume fraction
    Rs: 8.5e-6,              // Particle radius [m]
    asArea: (3.0 * 0.297) / 8.5e-6,
    csMax: 22860.0,          // Max solid Li concentration [mol/m^3]
    theta0: 0.17,            // Initial stoichiometry
    DsRef: 1.0e-13,          // Solid diffusivity at 298.15 K [m^2/s]
    kRef: 2.0e-11,           // Reaction rate constant [m^2.5/(mol^0.5*s)]
    EaD: 30000.0,
    EaK: 20000.0,
    rho: 2500.0,
    cp: 700.0,
    lambda: 5.0,
  },
  {
    id: REGION_IDS.POS_CC,
    name: 'Pos Current Collector (Al)',
    shortName: 'PosCC',
    thickness: 10e-6,        // 0.001 cm
    sigma: 3.8e7,            // 3.8e5 S/cm (Table 5)
    epsE: 0.0,
    epsS: 0.0,
    Rs: 0.0,
    asArea: 0.0,
    rho: 2702.0,             // Aluminum density [kg/m^3]
    cp: 903.0,               // Aluminum specific heat [J/(kg*K)]
    lambda: 238.0,           // Aluminum thermal conductivity [W/(m*K)]
  },
];

export const ELECTROLYTE_PARAMS = {
  ce0: 2000.0,               // Initial salt concentration [mol/m^3] (2.0 M)
  tPlus: 0.363,              // Li+ transference number (Doyle 1996)
  DeRef: 7.5e-11,            // Electrolyte salt diffusivity at 298.15 K [m^2/s]
  EaDe: 16000.0,             // Electrolyte diffusion activation energy [J/mol]
  EaKappa: 16500.0,          // Electrolyte conductivity activation energy [J/mol]
  brugg: 1.5,                // Bruggeman tortuosity exponent
  thermFactor: 1.0,          // (1 + d ln f_+- / d ln c_e)
};

/**
 * Geometry & Boundary Condition Presets directly from Main/final.pdf.
 */
export const GEOMETRY_PRESETS = {
  report_fig28: {
    id: 'report_fig28',
    label: 'Report Fig. 28 — 6 × 4 cm (Failed in 2013)',
    subtitle: '6 cm × 4 cm × 0.0346 cm · BC1 Opposite Corners (AR = 173)',
    Ly: 0.06,
    Lz: 0.04,
    bc: 'BC1',
    defaultCRate: 2.0,
    defaultThermalBC: 'convective',
    description:
      'The exact 3D geometry where ILU(1)+GMRES failed to converge across 6 MPI processors in Section 3.3 (Figure 28) of the original report.',
  },
  report_fig26: {
    id: 'report_fig26',
    label: 'Report Fig. 26 — 2 × 2 cm Pouch Cell',
    subtitle: '2 cm × 2 cm × 0.0346 cm · BC2 Top Tabs (AR = 58)',
    Ly: 0.02,
    Lz: 0.02,
    bc: 'BC2',
    defaultCRate: 1.0,
    defaultThermalBC: 'convective',
    description:
      'Matches the 2 cm × 2 cm × 0.0346 cm domain from Sections 3.1, 3.3 (Figs. 26–27) and 3.4 (ILUT false residual convergence, Fig. 29).',
  },
  report_tab18: {
    id: 'report_tab18',
    label: 'Report Table 18 — 1 × 1 cm Verification Cell',
    subtitle: '1 cm × 1 cm × 0.0346 cm · BC2 Top Tabs (AR = 29)',
    Ly: 0.01,
    Lz: 0.01,
    bc: 'BC2',
    defaultCRate: 1.0,
    defaultThermalBC: 'convective',
    description:
      'Matches the 1 cm × 1 cm 1C validation cell used to compare 2D and 3D extrema in Section 3.3.1 (Tables 14–18).',
  },
  report_fig9: {
    id: 'report_fig9',
    label: 'Report Fig. 9 — 50 × 20 cm EV Pouch Cell',
    subtitle: '50 cm × 20 cm × 0.0346 cm · BC1 (AR = 1,445!)',
    Ly: 0.50,
    Lz: 0.20,
    bc: 'BC1',
    defaultCRate: 2.0,
    defaultThermalBC: 'top_ambient',
    description:
      'Extreme EV-scale aspect ratio (1,445:1) from Section 2.3 (Fig. 9), where thermal & ohmic gradients along the 50 cm current collector dominate.',
  },
};

export const MESH_TIERS = {
  fast: {
    id: 'fast',
    label: 'Fast Interactive (3,360 cells)',
    nxCounts: [3, 6, 4, 8, 3],
    ny: 14,
    nz: 10,
  },
  standard: {
    id: 'standard',
    label: 'Standard (7,200 cells)',
    nxCounts: [4, 8, 5, 9, 4],
    ny: 20,
    nz: 12,
  },
  fine: {
    id: 'fine',
    label: 'Fine Table-14 Tier (13,824 cells)',
    nxCounts: [4, 10, 5, 13, 4],
    ny: 24,
    nz: 16,
  },
};

/**
 * Doyle / Dualfoil Open-Circuit Potential U_neg(theta) for Li_x C_6 [V vs Li/Li+].
 */
export function uNeg(theta) {
  const th = Math.max(0.005, Math.min(0.995, theta));
  return (
    0.7222 +
    0.1387 * th +
    0.029 * Math.sqrt(th) -
    0.0172 / th +
    0.0019 / Math.pow(th, 1.5) +
    0.2808 * Math.exp(0.9 - 15.0 * th) -
    0.7984 * Math.exp(0.4465 * th - 0.4108)
  );
}

/**
 * Entropic coefficient dU_neg/dT [V/K] (approximate reversible heat fit for Li_x C_6).
 */
export function dUNegDT(theta) {
  const th = Math.max(0.01, Math.min(0.99, theta));
  return -0.00015 * Math.sin(Math.PI * th) + 0.00005 * (0.5 - th);
}

/**
 * Doyle (1995/1996) Open-Circuit Potential U_pos(theta) for Li_y Mn_2 O_4 [V vs Li/Li+].
 */
export function uPos(theta) {
  const th = Math.max(0.05, Math.min(0.99, theta));
  return (
    4.19829 +
    0.0565661 * Math.tanh(-14.5546 * th + 8.60942) -
    0.0275479 / (Math.pow(0.998432 - th, 0.492465) - 1.90111) -
    0.157123 * Math.exp(-0.04738 * Math.pow(th, 8)) +
    0.810239 * Math.exp(-40.0 * (th - 0.133875))
  );
}

/**
 * Entropic coefficient dU_pos/dT [V/K] for Li_y Mn_2 O_4.
 */
export function dUPosDT(theta) {
  const th = Math.max(0.05, Math.min(0.99, theta));
  return -0.00025 + 0.00018 * (th - 0.5);
}

/**
 * Doyle (1996) 4th-order polynomial bulk ionic conductivity kappa(ce, T) [S/m].
 */
export function kappaBulk(ce, T) {
  const c = Math.max(50.0, Math.min(4000.0, ce));
  const c2 = c * c;
  const k25 =
    4.1253e-2 +
    5.007e-4 * c -
    4.7212e-7 * c2 +
    1.5094e-10 * c2 * c -
    1.6018e-14 * c2 * c2;
  const kSafe = Math.max(0.02, k25);
  return kSafe * Math.exp((ELECTROLYTE_PARAMS.EaKappa / R_GAS) * (1.0 / T_AMBIENT - 1.0 / T));
}

/**
 * Diffusional ionic conductivity kappa_D^eff [A/m] in Eq. (6) & Eq. (10).
 * Report Eq. (6): div(kappa_eff grad phi_e) + div(kappa_D_eff grad ln c_e) + j^Li = 0
 * Note: Under concentrated solution theory, kappa_D^eff = (2*R*T*kappa^eff / F)*(t_+^0 - 1)*(1 + d ln f/d ln c).
 */
export function kappaDiffusional(kappaEff, T) {
  return (
    ((2.0 * R_GAS * T * kappaEff) / F_CONST) *
    (ELECTROLYTE_PARAMS.tPlus - 1.0) *
    ELECTROLYTE_PARAMS.thermFactor
  );
}

/**
 * Temperature-dependent solid diffusivity D_s(T) [m^2/s].
 */
export function solidDiffusivity(spec, T) {
  return spec.DsRef * Math.exp((spec.EaD / R_GAS) * (1.0 / T_AMBIENT - 1.0 / T));
}

/**
 * Temperature-dependent reaction rate constant k_r(T).
 */
export function reactionRate(spec, T) {
  return spec.kRef * Math.exp((spec.EaK / R_GAS) * (1.0 / T_AMBIENT - 1.0 / T));
}

/**
 * Effective electrolyte diffusivity D_e^eff(epsE, T) [m^2/s].
 */
export function electrolyteDiffusivity(epsE, T) {
  const tort = Math.pow(epsE, ELECTROLYTE_PARAMS.brugg);
  return (
    ELECTROLYTE_PARAMS.DeRef *
    tort *
    Math.exp((ELECTROLYTE_PARAMS.EaDe / R_GAS) * (1.0 / T_AMBIENT - 1.0 / T))
  );
}
