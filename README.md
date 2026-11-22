# 3D 5-Region Lithium-Ion Battery Breakthrough Simulator & Monograph

Interactive 3D coupled electrochemical-thermal simulation (`Li_x C_6 | Li_y Mn_2 O_4`) solving the full 5-region, 6-equation Doyle-Fuller-Newman / Wang-Gu-Liaw battery model across realistic high-aspect-ratio pouch-cell geometries with copper and aluminum current collectors.

Developed to resolve the 3D matrix ill-conditioning bottleneck documented in Anderson & Kapadia (2013), *"Numerical Solution of Lithium Batteries"* (UTC Center of Excellence in Applied Computational Science and Engineering, File #14).

## The Mathematical & Numerical Breakthrough

1. **Fixing Section 4 ("Hyperbolicizing" Elliptic Equations)**:
   - Replacing isotropic $T_r = L^2/\sigma$ with conductivity-proportional anisotropic relaxation $T_{r,k} = \sigma / a_k^2$ (moving $\sigma$ from numerator to denominator $1/\sigma$).
   - Discretizing on a staggered dual (Raviart-Thomas / MAC) flux grid with integrated face currents $F_f = (\sigma\nabla\phi)\cdot\mathbf{A}_f$ (no continuous nodal tangential fluxes, no upwind artificial resistivity).
   - Applying symmetric diagonal resistance scaling $W_f = (1 + R_f/R_0)^{-1/2}$, collapsing matrix condition numbers by $>10^7\times$ (e.g., $7.9\times 10^{10} \to 1.8\times 10^3$ at aspect ratio $1,445$).
   - Solving with Block-Schur Preconditioned MINRES whose eigenvalues lie strictly inside $[-\varphi, -1] \cup [\varphi^{-1}, 1]$ ($\kappa \le \varphi^2 \approx 2.618$).

2. **Fixing Section 2.2 (Multigrid) & Section 3 (3D Newton-Krylov)**:
   - Active-region DOF compaction + positive-domain gauge shift ($\tilde{\phi}_s^+ = \phi_s^+ - V_{\text{ref}}$) eliminating $4.095\text{ V}$ cancellation inside the aluminum current collector.
   - Exact coupled Butler-Volmer Jacobian forming a strictly SPD 3D resistor-network graph Laplacian.
   - Coupled through-thickness ($x$-line) pentadiagonal $LDL^T$ smoother paired with conservative 2D $(y,z)$ Galerkin semi-coarsening multigrid.

## Running Locally (Zero Build Dependencies)

```bash
cd battery-3d-sim
python3 -m http.server 8080
# Open http://localhost:8080 in any modern browser
```

## Verification Suite

```bash
python3 tests/verify_solver.py
```

Verifies condition-number collapse across aspect ratios up to $1,445$, 2nd-order Method of Manufactured Solutions convergence across discontinuous interfaces, machine-precision equivalence between both 3D breakthrough solvers on Report Fig. 28 ($6\text{ cm} \times 4\text{ cm} \times 0.0346\text{ cm}$), and full 6-equation charge and lithium-salt conservation.
