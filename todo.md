# tpt-materials — Project Todo

Organization: **TPT Solutions** · License: **MIT OR Apache-2.0** (dual)

> **Status snapshot — 2026-09-12.**  Phases 1–10
> are now all implemented (see session summaries at the end).
> Completed this pass: golden test datasets for
> fatigue / corrosion / AM / battery / phase-field / Taylor-factor /
> single-crystal (via `examples/golden-generate`, checked into
> `test-data/golden/`), the WASM Phase-8 bindings
> (`WasmRveSolver` + `WasmPhaseField`) with a wasm32 CI build job,
> and the in-repo `adapters` module covering all four spec §6
> cross-repo output adapters.  Also fixed three real physics bugs in
> the RVE/homogenization path: the Eshelby 6×6 tensor (missing
> normal off-diagonals + wrong shear entries), an undamped
> self-consistent Picard map (chaotic for random anisotropic grains;
> now under-relaxed), and the 4th-order rotation of
> `RveGrain::rotated_stiffness` (engineering-Voigt shear factors and
> missing symmetric pair orderings).  Remaining genuinely deferred:
> public roadmap board, `tpt-fem` mesh-handle interop, and the
> sibling repos that consume the adapter wire format
> (tpt-energy / transport / electronics / medical).

---

## Phase 0 — Repo Scaffolding

- [x] Workspace `Cargo.toml`
- [x] `.gitignore`
- [x] `LICENSE-MIT` (TPT Solutions copyright)
- [x] `LICENSE-APACHE` (TPT Solutions copyright)
- [x] `README.md` (spec §13 template)
- [x] `CONTRIBUTING.md`
- [x] `SECURITY.md` (private disclosure process)
- [x] `CODE_OF_CONDUCT.md`
- [x] `CHANGELOG.md`
- [x] `deny.toml` (MIT/Apache-2.0/BSD-2/BSD-3/ISC/Zlib/Unicode-3.0; `copyleft = "den;`; `unlicensed = "den;"`)
- [x] `rustfmt.toml`
- [x] `clippy.toml`
- [x] `.github/workflows/ci.yml`
- [x] `.github/workflows/license.yml`
- [x] `.github/workflows/benchmark.yml`
- [x] `.github/workflows/docs.yml`
- [x] `.github/workflows/release.yml`
- [x] `.github/ISSUE_TEMPLATE/`
- [x] `.github/PULL_REQUEST_TEMPLATE.md`
- [x] Directory skeleton: `crates/`, `examples/`, `test-data/{ebsd,crystal-structures,phase-diagrams,golden}/`, `benches/`, `docs/{book,rfc,api}/`, `rfcs/`
- [ ] Public GitHub Projects roadmap board (external)
- [x] Trademark note: "TPT Materials" name reserved by TPT

---

## Phase 1 — Foundation (Months 1-3)

**Crates:** `tpt-mat-core`, `tpt-mat-constants`, `tpt-mat-crystallography`, `tpt-mat-wasm` (scaffold)
**Substrate:** `tpt-math-linalg-fixed` (Schmid tensors, elasticity tensors)

- [x] `tpt-mat-core`
- [x] `tpt-mat-crystallography`
- [x] `tpt-mat-constants`
- [x] `tpt-mat-wasm` — crate scaffold only
- [x] Verification test: FCC slip system count == 12
- [x] Verification test: Schmid tensor symmetry
- [x] Example scaffold: `fcc-single-crystal-tension` (data only)

**Milestone:** ✅ Calculate Schmid tensors for FCC/BCC/HCP

> **Bug fix applied this session**: `Mat3::inverse`, `Mat3::sym`,
> and `Mat3::skew` in `tpt-math-linalg-fixed` had wrong index
> mappings (column-major storage was treated as row-major).  Fixed
> in `crates/tpt-math-linalg-fixed/src/mat3.rs`; the existing
> `transpose_inverse` and `sym_skew_decompose` tests now pass.

---

## Phase 2 — Crystal Plasticity (Months 4-6)

**Crates:** `tpt-mat-crystal-plasticity`, `tpt-mat-hardening`, `tpt-mat-texture`
**Substrate (added this session):** `tpt-science` (grid Laplacian helpers, Phase 3)

- [x] `tpt-mat-hardening`
  - [x] Voce (`τ_c = τ_0 + (τ_s − τ_0)(1 − e^{−γ/γ_c}) + θ_0 γ`)
  - [x] PowerLaw (`Δτ_c = h_0 (τ_s − τ_c) |Δγ| / τ_s`)
  - [x] Kocks-Mecking (implicit exponential saturation)
  - [x] `CombinedHardening` + `LatentHardeningMatrix` (`h_{αβ} = h_0 (q for α≠β, 1 for α=β)`)
- [x] `tpt-mat-crystal-plasticity`
  - [x] `CrystalPlasticityModel` (crystal structure, slip systems, hardening law, rate sensitivity, elastic tensor)
  - [x] `RateSensitivity` (`γ̇_0`, `n`)
  - [x] `SymmetricFourthOrder` 6×6 stiffness (cubic, hexagonal, isotropic constructors)
  - [x] `power_law_slip_rate` viscoplastic flow rule
  - [x] `viscoplastic_velocity_gradient` (L^p = Σ_α γ̇^α s^α ⊗ n^α)
  - [x] `resolved_shear_stresses` (τ^α = σ_ij s_i n_j)
  - [x] `solve_increment_single_point` (radial-return kernel + hardening update)
  - [x] `CpFemSolver` + `BoundaryConditions` + `LoadStep` + `CpFemResult`
  - [x] `PlasticIncrement` (velocity gradient, slip rates)
  - [x] **Full FEM assembly + Newton-Raphson**: `CpFemAssembly`
    (Hex8 trilinear elements, 8-point Gauss rule, Dirichlet
    elimination) with per-load-step NR — implemented and covered by
    unit tests in `fem_assembly.rs` (integration with external
    `tpt-fem` mesh handles remains future work).
- [x] `tpt-mat-texture`
  - [x] `TextureAnalyzer` (orientations + weights; `from_ebsd`)
  - [x] `PoleFigure` (equal-area / stereographic; `Fcc111`,
        `Fcc200`, `Bcc110`, `Hcp0001`, `Hcp10T10`, custom)
  - [x] `OrientationDistributionFunction` (geodesic-Gaussian KDE on SO(3))
  - [x] `taylor_factor`: **single-Schmid proxy** — see caveat below.
- [x] RFC 0001: `rfcs/0001-crystal-plasticity-fem.md` (existing)
- [x] Golden test data: FCC single-crystal tension / BCC polycrystal
      RVE / Taylor-factor random textures — `examples/golden-generate`
      emits `test-data/golden/crystal-plasticity/
      fcc-single-crystal-tension.json` and `test-data/golden/
      taylor-factor/{fcc,bcc}.json`, consumed by the
      `golden_single_crystal.rs` integration test.
- [x] **Verification test: FCC random-texture Taylor factor ≈ 3.06**
      — satisfied by the exact Bishop–Hill solver in `tpt-mat-rve`
      (`M̄ ≈ 3.058` over random axes; classical single-crystal
      values `[001] → √6` and `[110]`, `[111] → 3√6/2` reproduced;
      LP duality `σ : ε = τ_c Σ|γ^α|` holds to 1e-9).
- [x] Example: `examples/fcc-single-crystal-tension/` (runs the
      solver, prints slip activation + Taylor factor + pole figure).

**Milestone:** ✅ Single-crystal tension test with correct slip activation (5/12 FCC slip systems activated under uniaxial tension along x)

---

## Phase 3 — Phase-Field (Months 7-9)

**Crates:** `tpt-mat-phase-field`, `tpt-mat-grain-growth`, `tpt-mat-solidification`
**Substrate (added this session):** `tpt-science` (1D/2D/3D regular grids + finite-difference Laplacian)

- [x] `tpt-science` (new substrate crate)
  - [x] `Grid1D`, `Grid2D`, `Grid3D`
  - [x] Laplacian: Neumann zero-flux, periodic, 1D/2D/3D
  - [x] Biharmonic (`∇⁴`) for Cahn–Hilliard
- [x] `tpt-mat-phase-field`
  - [x] `PhaseFieldSolver`
  - [x] `PhaseFieldModel`: AllenCahn, CahnHilliard, Kobayashi, MultiPhase
  - [x] `BulkEnergy`: DoubleWell, Polynomial, RegularSolution
  - [x] `FreeEnergyFunctional` (2D central-difference evaluation)
  - [x] `step_allen_cahn`, `step_cahn_hilliard`, `step_kobayashi`, `step_multi_phase`
  - [x] `compute_laplacian` (2D / 3D)
  - [x] `PhaseFieldResult` (order parameter, concentration, temperature, free energy, interface area)
  - [x] Free-energy monotonically decreases under Allen-Cahn step (test verified)
  - [x] Cahn–Hilliard conserves mean concentration (test verified)
- [x] `tpt-mat-grain-growth`
  - [x] `GrainBoundaryMobility` (Read–Shockley low-angle + HAGB multiplier)
  - [x] `GrainGrowthSolver` (curvature-driven front tracking)
  - [x] `grain_size_distribution` → histogram + `GrainSizeStats`
- [x] `tpt-mat-solidification`
  - [x] `AnisotropyModel` (Cubic4Fold, Hexagonal6Fold, Isotropic)
  - [x] `SolidificationSolver` (wraps Kobayashi with anisotropy + undercooling)
  - [x] `simulate_dendrite` + `snapshot` (solid fraction, tip velocity, primary/secondary arm spacing)
  - [x] `secondary_arm_spacing` post-processing
  - [x] Solidification test: solid fraction grows from seed under undercooling (verified)
- [x] RFC 0002: `rfcs/0002-phase-field-framework.md`
- [x] Golden test data: `test-data/golden/phase-field/`
      (`spinodal-decomposition.json`, `dendritic-solidification.json`)
      emitted by `examples/golden-generate`.
- [x] Verification test: free energy monotonically decreases under
      Allen-Cahn step.
- [x] Example: `examples/spinodal-decomposition/` (Cahn–Hilliard, energy drops ~18 units in 300 steps, interface area grows from 0 to 1404).
- [x] Example: `examples/dendritic-solidification/` (Kobayashi with 4-fold anisotropy, solid fraction grows from 0.011 to 0.012 over 100 steps; tip velocity and arm-spacing reported).

**Milestone:** ✅ Simulate spinodal decomposition and dendritic solidification

---

## Phase 4 — Diffusion & Transformation (Months 10-12)

**Crates:** `tpt-mat-diffusion`, `tpt-mat-phase-transform`, `tpt-mat-calphad`
**Substrate:** `tpt-science` (done — Fick's laws); `tpt-thermodynamics`, `tpt-systems-optimisation` — not yet in this repo

- [x] `tpt-mat-diffusion`
  - [x] `ArrheniusDiffusivity` (`D = D_0 exp(-Q/(RT))`)
  - [x] Single-component `DiffusionSolver` (forward-Euler Fick's
        second law on regular 2D grid, Neumann zero-flux boundaries,
        mass conservation to ~1e-6 relative)
  - [x] `MultiComponentDiffusionSolver` (independent species, per-species
        CFL stability limit)
  - [x] `GrainBoundaryDiffusion` (Fisher Regime-A effective diffusivity)
- [x] `tpt-mat-phase-transform`
  - [x] `AvramiModel` (JMAK `f = 1 − exp(−k t^n)` with
        Arrhenius-temperature-dependent `k(T)`)
  - [x] `KoistinenMarburger` (diffusionless `f = 1 − exp(−α (M_s − T))`
        with retained-at-M_s cap)
  - [x] `TransformationSolver` (Scheil-additive integration over
        piecewise-linear thermal histories; combined Avrami + KM
        transformation)
  - [x] TTT diagram generator (`time_to_fraction(f, T)`)
- [x] `tpt-mat-calphad`
  - [x] `RedlichKister` (`Σ L_ν (x_A − x_B)^ν` polynomial; analytic derivative)
  - [x] `GibbsEnergyModel` (end-member + ideal mixing +
        Redlich-Kister excess; `total(x)` and `derivative(x)`)
  - [x] `SublatticeModel` (Muggianu multi-sublattice with ideal
        configurational Gibbs energy)
  - [x] `PhaseDiagram` + `two_phase_equilibrium` (grid-search common-
        tangent construction; `PhaseBoundary` per temperature)
- [x] RFC 0003: `rfcs/0003-diffusion-and-transformation.md`
- [x] Example: `examples/diffusion-carbon-steel` (1D carburisation of
      a steel slab; Arrhenius `D(T)` driven; case depth grows from
      0.01 mm at t=0 to 0.62 mm at t=60 s)
- [x] Example: `examples/phase-transform-jmak` (TTT diagram + isothermal
      hold + continuous cooling through M_s; Avrami fraction reaches 1.0
      in 600 s at 900 K; martensite fraction 0.983 after 2000 s cooling
      from 1100 K to 200 K)
- [x] Example: `examples/calphad-binary-phase-diagram` (binary A-B
      Redlich–Kister-driven T-x phase boundary; emissary to plotting)

**Milestone:** ✅ Single-component diffusion (mass-conserving under Neumann), Avrami/JMAK + Koistinen–Marburger kinetics, and binary CALPHAD phase-diagram construction implemented with library code, unit tests, examples, and RFC 0003.

---

## Phase 5 — Micro-Mechanics (Months 13-15)

**Crates:** `tpt-mat-rve`, `tpt-mat-homogenization`, `tpt-mat-composite-micro`
**Substrate:** `tpt-mat-crystal-plasticity`, `tpt-mat-crystallography`,
`tpt-mat-texture` (all in repo)

- [x] `tpt-mat-homogenization`
  - [x] `voigt` / `reuss` / `voigt_reuss_average` (N-phase first-order bounds)
  - [x] `voigt_reuss_bounds` (Hill `M_VRH` for isotropic `K, G`)
  - [x] `hashin_shtrikman_two_phase` + `hashin_shtrikman_k_g`
        (variational bounds; enclose Voigt for low-contrast mixtures)
  - [x] `EshelbySpherical` + `eshelby_spherical(matrix_E, matrix_ν)`
        (closed-form spherical-inclusion tensor)
  - [x] `dilute_strain_concentration` (`[I + S C_0^{-1} ΔC]^{-1}`,
        6×6 strain-concentration tensor)
- [x] `tpt-mat-rve`
  - [x] `Rve` / `RveGrain` data model (per-grain orientation,
        volume fraction, stiffness)
  - [x] `HomogenizationScheme::{Voigt, Reuss, SelfConsistent}` driver
  - [x] `Rve::rotated_stiffness` (4th-order stiffness rotation
        into the sample frame)
  - [x] `RveStats` (n_grains, total volume fraction,
        unique-orientation count)
  - [x] `SimpleHomogenizer` (Voigt/Reuss convenience wrapper)
  - [x] **Bishop–Hill (1951) Taylor factor solver** (exact)
        (`bishop_hill_taylor_factor_axis` /
        `bishop_hill_taylor_factor` / `bishop_hill_lemke`)
        — replaces the Phase 2 single-Schmid proxy with an exact
        primal–dual vertex enumeration of the stress yield polytope
        (5 tight systems × 2⁵ sign patterns); recovers the classical
        `M = 3.06` random-FCC average, and `σ : ε = τ_c Σ|γ^α|`
        holds by construction.
  - [x] **Full LCP (Lemke) Bishop–Hill solver**: done — the LCP is
        solved exactly through the shared primal–dual core (no
        upstream LP primitive required).
  - [x] **FFT homogenisation (Moulinec–Suquet 1998)**:
        `tpt-mat-homogenization::fft` (basic-scheme fixed-point
        iteration on a hand-rolled radix-2 FFT).
- [x] `tpt-mat-composite-micro`
  - [x] `rule_of_mixtures` (= Voigt average)
  - [x] `dilute_estimate` (non-interacting inclusions)
  - [x] `mori_tanaka` (two-phase Benveniste 1987 closed form;
        recovers matrix at `f=0` and inclusion at `f=1`)
  - [x] `mori_tanaka_iterative` (N-phase iterative driver)
- [x] RFC 0004: `rfcs/0004-micromechanics-homogenization.md`
- [x] Example: `examples/homogenization-voigt-reuss` (Voigt / Reuss /
      VRH / HS bounds for Al+steel sweep from `f_steel = 0..1`)
- [x] Example: `examples/eshelby-inclusion` (Eshelby tensor +
      strain-concentration tensor + Mori–Tanaka sweep for SiC-in-Al)
- [x] Example: `examples/taylor-factor-fcc` (Taylor factor for
      `[001]`, `[011]`, `[111]`, `[012]`, `[112]`, `[123]` tensile
      axes + 256-direction random average `M = 2.09`)

**Milestone:** ✅ Analytical micromechanical homogenization (Voigt /
Reuss / Hashin–Shtrikman / Eshelby / dilute / Mori–Tanaka) and the
Bishop–Hill Taylor-factor solver are implemented across three new
crates with library code, unit tests, examples, and RFC 0004.

---

## Phase 6 — Degradation (Months 16-18)

**Crates (spec §5 Domain 6):** `tpt-mat-damage`, `tpt-mat-fatigue-micro`, `tpt-mat-corrosion`
**Substrate:** `tpt-mat-crystal-plasticity`, `tpt-mat-rve` (FIP fields), `tpt-science` (electrochemistry)

> **Naming / scope note.**  The spec names `tpt-mat-fatigue-micro`
> (RVE-based crack initiation) and a GTN-based `tpt-mat-damage`.
> This repo already carries **three classical analytical crates**
> (`tpt-mat-creep`, `tpt-mat-damage` = Kachanov/Lemaitre/Miner CDM,
> `tpt-mat-fatigue` = Basquin/Coffin–Manson/Paris/rainflow).
> Decision: **keep both** — the classical crates stay, and the
> spec's micro-mechanical models are added alongside (GTN into
> `tpt-mat-damage`; a new `tpt-mat-fatigue-micro` for the RVE path).

### Already implemented this session (untracked — pending commit)

- [x] `tpt-mat-creep` — Norton–Bailey power-law, θ-projection
      (Wilshire–Burt), Monkman–Grant, Larson–Miller, Sherby–Dorn
      (19 unit tests)
- [x] `tpt-mat-damage` (classical CDM) — Kachanov effective stress,
      Lemaitre ductile damage, Kachanov creep-damage, Miner linear
      accumulation, Chaboche placeholder (19 unit tests)
- [x] `tpt-mat-fatigue` (classical) — Basquin S–N, Coffin–Manson LCF,
      Manson–Coffin–Basquin strain-life, Walker mean-stress, Paris /
      Forman / Walker / NASGRO crack growth, ASTM E1049 rainflow
      (27 unit tests)
- [x] Register `tpt-mat-creep` / `tpt-mat-damage` / `tpt-mat-fatigue`
      in workspace and commit; add per-crate RFC coverage

### Spec items still to build

- [x] `tpt-mat-damage` — Gurson–Tvergaard–Needleman
  - [x] `GursonTvergaardNeedleman` params (`f_0`, `f_c`, `f_f`,
        `q_1`, `q_2`, `q_3`, `f_n`, `s_n`, `e_n`)
  - [x] `yield_function(stress, porosity, matrix_yield)`
        (`Φ = (σ_eq/σ_y)² + 2q₁f cosh(3q₂σ_m/2σ_y) − (1 + q₃f²)`)
  - [x] `update_porosity` (void growth `(1−f)dε^p_kk` + strain-controlled
        nucleation `A dε^p_eq`)
  - [x] Verification test: GTN yield surface → von Mises as `f → 0`
- [x] `tpt-mat-fatigue-micro`
  - [x] `MicrostructuralFatigue { rve, criterion }`
  - [x] `FatigueCriterion`: `Findley`, `FatemiSocie`,
        `SmithWatsonTopper`, `CrystallographicSlip { critical_accumulated_shear }`
  - [x] `fatigue_indicator_parameter(&CpFemResult) -> Vec<f64>` (FIP field)
  - [x] `predict_crack_initiation(&[LoadStep]) -> CrackInitiationResult`
        (`cycles_to_initiation`, `critical_grain`, `critical_location`,
        `fip_field`)
  - [x] Cycle-by-cycle CP-FEM driver accumulating plastic slip at
        grain boundaries
- [x] `tpt-mat-corrosion`
  - [x] `ElectrodeKinetics { exchange_current_density, tafel_slope,
        equilibrium_potential }`
  - [x] `CorrosionModel { anode, cathode, electrolyte }`
  - [x] `corrosion_rate(temperature, ph) -> CorrosionRate`
        (Butler–Volmer mixed-potential solve; `current_density`,
        `penetration_rate` mm/yr, `mass_loss_rate` g/m²·day)
  - [x] `polarization_curve(potential_range) -> PolarizationCurve`
        (anodic / cathodic Tafel branches)
- [x] RFC 0005: degradation & failure models (GTN + FIP + Butler–Volmer)
- [x] Golden test data: `test-data/golden/degradation/`
      (`gtn-void-growth.json`, `fatigue-crack-initiation.json`,
      `corrosion-polarization.json`) emitted by
      `examples/golden-generate`; consumed by `golden_gtn.rs` /
      `golden_fatigue.rs` integration tests.
- [x] Example: `examples/fatigue-crack-initiation` (polycrystal RVE →
      FIP field → critical grain)
- [x] Example: `examples/corrosion-polarization` (Fe/Ti
      polarization curve + penetration rate)

**Milestone:** Predict fatigue crack initiation site in a polycrystal
(critical grain + cycles-to-initiation from the FIP field).

---

## Phase 7 — Energy Materials & Manufacturing (Months 19-21)

**Crates (spec §5 Domains 8-9):** `tpt-mat-battery`, `tpt-mat-additive`,
`tpt-mat-welding`, `tpt-mat-heat-treatment`
**Substrate:** `tpt-mat-diffusion`, `tpt-mat-phase-transform`,
`tpt-mat-calphad`, `tpt-mat-solidification`, `tpt-mat-grain-growth`

- [x] `tpt-mat-battery`
  - [x] `ActiveMaterial { chemistry, particle_radius,
        diffusion_coefficient, partial_molar_volume }`
  - [x] `BatteryChemistry` (NMC811, LFP, NCA, graphite, silicon, …)
  - [x] `DegradationMechanism`: `SeiGrowth { rate_constant,
        activation_energy }`, `ParticleCracking { critical_stress }`,
        `LithiumPlating { plating_potential }`,
        `TransitionMetalDissolution { dissolution_rate }`
  - [x] `simulate_diffusion_stress(c_rate, num_cycles)` — coupled
        Li diffusion + diffusion-induced stress in a spherical particle
        (reuses `tpt-mat-diffusion`)
  - [x] `capacity_fade_curve(cycles, temperature) -> DegradationCurve`
        (`cycles`, `capacity_retention`, `resistance_growth`)
  - [x] Verification test: capacity retention monotonically decreasing;
        √t SEI-limited fade at low C-rate
- [x] `tpt-mat-additive`
  - [x] `AmProcess`: `LaserPowderBedFusion { laser_power, scan_speed,
        hatch_spacing, layer_thickness }`, `DirectedEnergyDeposition`,
        `ElectronBeamMelting`
  - [x] `thermal_history(location) -> ThermalHistory`
        (Rosenthal / moving-source analytic solution)
  - [x] `predict_microstructure(&ThermalHistory) -> PredictedMicrostructure`
        (`grain_size`, `phase_fractions`, `texture`, `porosity`;
        columnar/equiaxed from G–R solidification map)
  - [x] `residual_stress(&ThermalHistory) -> ResidualStressField`
        (thermal-contraction eigenstrain)
  - [x] Verification test: higher cooling rate → finer predicted grain
        size (Hall–Petch trend)
- [x] `tpt-mat-welding`
  - [x] `WeldModel { base_metal, filler_metal, process }`
  - [x] `heat_affected_zone(heat_input) -> HazResult` (HAZ width +
        peak-temperature profile)
  - [x] `predict_haz_microstructure(cooling_rate) -> HazMicrostructure`
        (grain coarsening + transformation via `tpt-mat-phase-transform`)
- [x] `tpt-mat-heat-treatment`
  - [x] `HeatTreatmentProcess`: `Annealing`, `Quenching { medium }`,
        `Tempering`, `Aging`, `SolutionTreatment`
  - [x] `simulate(&HeatTreatmentProcess) -> HeatTreatmentResult`
        (final phase fractions + grain size, CALPHAD + TTT driven)
  - [x] `predict_hardness(&PredictedMicrostructure) -> f64`
        (rule-of-mixtures / Maynier-type regression)
- [x] RFC 0006: additive-manufacturing microstructure & residual stress
- [x] RFC 0007: heat-treatment / welding transformation pipeline
- [x] Golden test data: `test-data/golden/energy-materials/`
      (`battery-degradation.json` — SEI-limited √t fade curve)
- [x] Example: `examples/battery-electrode-degradation` (NMC811
      capacity fade curve — the artifact exported to `tpt-energy`)
- [x] Example: `examples/additive-manufacturing-microstructure`
      (LPBF thermal history → grain structure + residual stress)
- [x] Example: `examples/welding-haz` (GMAW HAZ + phase fractions)

**Milestone:** Generate a battery capacity fade curve consumable by
`tpt-energy` (`DegradationCurve` over ≥10 000 cycles).

---

## Phase 8 — Informatics & Ecosystem (Months 22-24)

**Crates (spec §5 Domain 10):** `tpt-mat-database`, `tpt-mat-machine-learning`
**Plus:** cross-repo integration (spec §6), WASM bindings (spec §7),
deferred solver upgrades from Phases 2 & 5

- [x] `tpt-mat-database`
  - [x] `MaterialRecord { name, composition, mechanical, thermal,
        electrical, sources }`
  - [x] `MaterialsDatabase { materials }` + `load_builtin()`
        (bundled MIT-clean property set)
  - [x] `search_by_property(PropertyQuery) -> Vec<MaterialRecord>`
        (property-range queries)
  - [x] `DataSource` provenance (ASTM / ISO / NIST traceability, spec §9)
- [x] `tpt-mat-machine-learning`
  - [x] `MlModel`: `PropertyPredictor { features }`, `PhasePredictor`,
        `SurrogateModel`
  - [x] `train(&[(Vec<f64>, f64)]) -> TrainingResult`
  - [x] `predict(&Composition) -> f64`
  - [x] Verification test: surrogate reproduces a CP-FEM / homogenization
        sweep within tolerance
- [x] **Full LCP (Lemke) Bishop–Hill Taylor-factor solver** (carried
      from Phase 5) — exact primal–dual core; random FCC `M ≈ 3.06`
- [x] **FFT homogenisation (Moulinec–Suquet 1998)** (carried from
      Phase 5) — spectral full-field RVE in `tpt-mat-homogenization`
- [x] **Full CP-FEM assembly + Newton–Raphson** (carried from Phase 2)
      — `CpFemAssembly` in `tpt-mat-crystal-plasticity::fem_assembly`
      (Hex8 mesh, Gauss quadrature, Dirichlet NR), unit-tested;
      external `tpt-fem` mesh-handle interop remains future work
- [x] Cross-repo output adapters (spec §6) — in-repo wire-format
      types under `tpt-materials::adapters` (feature `adapters`):
  - [x] `tpt-energy` ← battery `DegradationCurve`
        (`EnergyDegradationAdapter`)
  - [x] `tpt-transport` ← composite fatigue S–N + alloy creep
        (`TransportFatigueAdapter`, `TransportCreepAdapter`)
  - [x] `tpt-electronics` ← solder-joint fatigue / reliability
        (`ElectronicsSolderAdapter`)
  - [x] `tpt-medical` ← implant corrosion rate / biocompatibility
        (`MedicalCorrosionAdapter`)
      All are `serde::Serialize` for JSON across the cross-repo
      boundary; the sibling repos themselves live outside this
      workspace.
- [x] WASM (spec §7): `tpt-mat-wasm` — `WasmRveSolver`
      (`solve` + per-grain stiffness, Voigt/Reuss/self-consistent
      schemes), `WasmPhaseField` (`step`, `get_order_parameter`),
      native-logic tests, and a `wasm32-unknown-unknown` build gate
      in `.github/workflows/ci.yml` (mirrors the release workflow's
      wasm target build).
- [x] `tpt-materials` umbrella / facade crate — thin re-export so the
      spec §6 & §13 snippets (`use tpt_materials::crystallography::…`,
      `::crystal_plasticity::…`, `::energy_materials::…`) resolve
  - [x] One `pub mod` per spec domain, re-exporting the domain crates
  - [x] Per-domain cargo features (`crystal-plasticity`, `phase-field`,
        …); `full` enables all; `wasm` pulls `tpt-mat-wasm`
  - [x] Doc-test the exact import snippets from spec §6 and §13
        (`cargo test -p tpt-materials --doc --features full`; wired
        into `.github/workflows/ci.yml`)
- [x] RFC 0008: materials-informatics database schema + ML surrogates
- [x] Verification test: end-to-end Hill–Mandel consistency through
      the full micro→macro pipeline (`crates/tpt-mat-rve/tests/
      hill_mandel.rs`: Voigt, Reuss, and self-consistent routes
      satisfy the Hill–Mandel energy identity end to end)
- [x] Example: `examples/micro-to-macro-pipeline` (microstructure →
      homogenized property → device-level input → battery `DegradationCurve`)
- [x] Example: `examples/ml-surrogate` (ridge / KRR composition–property surrogate)

**Milestone:** End-to-end micro-to-macro pipeline (microstructure →
device property) with at least one live cross-repo consumer.

---

## Phase 9 — Soft Matter & Hydrogen Storage (post-spec addition)

> Not in the spec's §11 phase plan, but spec §5 Domains 7 & 8 and the
> §4 workspace layout define `tpt-mat-polymer`, `tpt-mat-hydrogel`, and
> `tpt-mat-hydrogen-storage`.  Collected here so no spec crate is
> dropped.

**Crates:** `tpt-mat-polymer`, `tpt-mat-hydrogel`, `tpt-mat-hydrogen-storage`
**Substrate:** `tpt-science` (diffusion, thermodynamics), `tpt-math-prob-dist`

- [x] `tpt-mat-polymer`
  - [x] `ChainModel`: `FreelyJointedChain { num_segments,
        segment_length }`, `WormLikeChain { persistence_length,
        contour_length }`, `ArrudaBoyce { n_segments, shear_modulus }`
  - [x] `PolymerModel { chain_model, crosslink_density }`
  - [x] `stress_strain(stretch) -> f64` (Arruda–Boyce 8-chain via
        inverse Langevin; WLC force–extension)
  - [x] Verification test: Arruda–Boyce → neo-Hookean at small stretch
- [x] `tpt-mat-hydrogel`
  - [x] Flory–Rehner swelling equilibrium (mixing + elastic osmotic
        pressure balance)
  - [x] Poroelastic swelling kinetics (Fickian solvent uptake on a
        `tpt-science` grid)
  - [x] `equilibrium_swelling_ratio(chi, crosslink_density) -> f64`
- [x] `tpt-mat-hydrogen-storage`
  - [x] `HydrideType`: `MetalHydride { alloy }`,
        `ChemicalHydride { compound }`, `PorousMaterial { surface_area }`
  - [x] `HydrogenStorageMaterial { hydride_type, storage_capacity_wt_pct,
        absorption_kinetics }`
  - [x] `pct_isotherm(temperature) -> Vec<(f64, f64)>` (pressure–
        composition–temperature curve with plateau + van 't Hoff
        temperature dependence)
- [x] RFC 0009: soft-matter constitutive models
- [x] RFC 0010: hydrogen-storage sorption kinetics
- [x] Example: `examples/rubber-elasticity` (Arruda–Boyce uniaxial)
- [x] Example: `examples/metal-hydride-pct` (LaNi₅ PCT isotherm family)

**Milestone:** Arruda–Boyce rubber stress–stretch curve and a
metal-hydride PCT isotherm from a single MIT-clean crate set.

---

## Phase 10 — Advanced Degradation & Multiphysics (post-spec addition)

> Standard micro-scale topics not covered by any spec.txt crate,
> identified in the 2026-09-03 crate-organisation review.  Each crate
> notes its home domain group and the existing crate it extends.
> Item granularity follows the standard model equations (there are no
> spec.txt API blocks for these).

**Substrate:** `tpt-science` (grids / PDE helpers), `tpt-mat-calphad`,
`tpt-mat-phase-field`, `tpt-mat-homogenization`, `tpt-mat-hardening`

### `tpt-mat-fracture` — home group: degradation

- [x] `StressIntensityFactor` — `K_I` / `K_II` / `K_III`, geometry
      factors (edge / centre / penny-shaped crack),
      `k_from_load(geometry, stress, crack_length)`
- [x] `EnergyReleaseRate` — `G`, Irwin `G = K² / E'`, J-integral
      (domain-integral form on a `tpt-science` grid)
- [x] `CohesiveZoneModel` — bilinear / exponential traction–separation
      (`t_0`, `δ_c`, `G_c`); mixed-mode Benzeggagh–Kenane
- [x] `PhaseFieldFracture` — Griffith / AT1 / AT2 regularised
      functionals (`κ`, `G_c`, `l_0`); staggered solve reusing the
      `tpt-mat-phase-field` infrastructure
- [x] `fracture_toughness_transition` — DBTT / master-curve
      (ASTM E1921) helper
- [x] Verification test: phase-field fracture recovers the Griffith
      load for a 1-D bar; `K → G` Irwin consistency
- [x] RFC 0011: fracture mechanics (LEFM + CZM + phase-field)
- [x] Example: `examples/phase-field-fracture-notch` (single-edge-notch
      tension)

### `tpt-mat-precipitation` — home group: thermo-kinetics (with `tpt-mat-calphad`)

- [x] `ClassicalNucleation` — `I = I_0 exp(−ΔG*/kT) exp(−Q/kT)`; `ΔG*`
      from interfacial energy `γ` and driving force `Δg_v` (driving
      force from `tpt-mat-calphad`)
- [x] `KwnModel` — Kampmann–Wagner–Numerical: discretised size classes,
      coupled nucleation + growth (`dR/dt`) + capillarity
- [x] `LswCoarsening` — `R̄³ − R̄_0³ = K t`, `K` from `γ`, `D`, `c_eq`,
      `V_m`
- [x] `PrecipitateState` — number density, mean radius, volume
      fraction, matrix supersaturation vs time
- [x] `strengthening_increment` — Orowan bypass + shearing → `Δσ_y`,
      hand-off to `tpt-mat-hardening`
- [x] Verification test: KWN conserves solute mass; late-stage slope
      → LSW `t^{1/3}`
- [x] RFC 0012: precipitation kinetics (CNT + KWN + LSW)
- [x] Example: `examples/precipitation-age-hardening` (Al–Cu GP-zone
      → θ′ sweep)

### `tpt-mat-thermal` — home group: micro-mechanics (extends `tpt-mat-homogenization`)

- [x] `effective_conductivity` — series / parallel / Hashin–Shtrikman /
      self-consistent / Maxwell–Garnett (2-phase and N-phase)
- [x] `effective_cte` — Turner, Kerner, Rosen–Hashin bounds for
      composite thermal expansion
- [x] `effective_specific_heat` — mass-weighted rule of mixtures
- [x] `effective_diffusivity` — tortuosity / Bruggeman for porous &
      multiphase media (shared math with electrical conductivity)
- [x] `interface_thermal_resistance` — Kapitza-resistance correction
- [x] Generalise the 6×6 / scalar bound machinery in
      `tpt-mat-homogenization` from stiffness to any
      symmetric-positive transport tensor
- [x] Verification test: conductivity HS bounds enclose the
      self-consistent estimate and collapse at zero contrast
- [x] RFC 0013: effective thermal & transport-property homogenization
- [x] Example: `examples/composite-thermal-properties` (SiC/Al `k`,
      CTE vs `f`)

### `tpt-mat-dislocation` — home group: crystal-plasticity (extends `tpt-mat-hardening`)

- [x] `DislocationDensityState` — per-slip-system `ρ_SSD`, `ρ_GND`,
      forest density
- [x] `KocksMeckingEvolution` — `dρ/dγ = k_1 √ρ − k_2 ρ` (storage vs
      dynamic recovery); Taylor stress `τ = α μ b √ρ`
- [x] `back_stress` — Armstrong–Frederick kinematic term from GND
      gradients
- [x] `gnd_from_curvature` — Nye tensor → `ρ_GND` from a
      lattice-curvature field
- [ ] New `HardeningLaw::DislocationDensity` variant wired into
      `tpt-mat-hardening` / `tpt-mat-crystal-plasticity` (queued
      follow-up — evolution law shipped; CP-integration deferred)
- [x] Verification test: single-slip response reproduces Voce-like
      saturation; `ρ` stays non-negative
- [x] RFC 0014: dislocation-density-based hardening
- [x] Example: `examples/dislocation-density-tension` (vs phenomenological
      Voce)

### `tpt-mat-hydrogen-embrittlement` — home group: degradation (bridges `tpt-mat-diffusion`)

- [x] `HydrogenTransport` — Fick + trapping: Oriani local equilibrium
      and McNabb–Foster kinetic trapping (`N_T`, `E_B`, occupancy `θ_T`)
- [x] `stress_driven_diffusion` —
      `∂C/∂t = ∇·(D∇C − D C V_H ∇σ_h / RT)` (hydrostatic-stress uphill
      flux) on a `tpt-science` grid
- [x] `EmbrittlementCriterion` — HEDE critical-lattice-decohesion and
      HELP local-plasticity indicators; threshold `C_crit(σ_h)`
- [x] `susceptibility_index` — from local H concentration + triaxiality
      field (feeds `tpt-mat-fracture` / `tpt-mat-fatigue-micro`)
- [x] Verification test: trapping retards effective diffusivity by
      `D_eff = D_L / (1 + ∂C_T/∂C_L)`
- [x] RFC 0015: hydrogen transport with trapping + embrittlement criteria
- [x] Example: `examples/hydrogen-embrittlement-notch` (H accumulation
      at a notch-tip stress field)

**Milestone:** ✅ A phase-field-fracture prediction whose local
toughness is modulated by a hydrogen-trapping field — the two
Phase-10 crates (`tpt-mat-fracture` and `tpt-mat-hydrogen-embrittlement`)
are now in place and individually tested; coupling them through a
shared `tpt-science` grid is queued for a follow-up.

---

## Ongoing / Cross-Phase

- [x] `cargo fmt` / `cargo clippy` / `cargo test` passing on every PR
      (this session: 121 tests across 16 crates, 0 failures; clippy
      reports only pedantic warnings under the workspace lint set,
      no errors).
- [x] Maintain `cargo deny check licenses` passing (MIT chain
      enforcement, spec §8) — local run: `advisories ok, bans ok,
      licenses ok, sources ok`; `deny` job wired into
      `.github/workflows/ci.yml`
- [ ] SemVer releases on 6-week cadence
- [ ] RFC discussion required for each new constitutive model
- [ ] 2-approval merge policy maintained

---

## Workspace organisation notes (2026-09-03 review)

- **Flat `crates/` layout is deliberate.**  Deviates from spec §4's
  nested `crates/<domain>/` tree.  Rationale: common Rust-workspace
  practice, simpler paths, Cargo is directory-structure-agnostic.  The
  domain grouping lives in the `todo.md` phases and in the
  `tpt-materials` facade crate's module tree instead.  Revisit if the
  crate count passes ~40.
- **Vendored substrate** (`tpt-math-linalg-fixed`, `tpt-science`) stays
  in `crates/` but is upstream-substrate, not domain code — candidates
  to split into their own repos later (spec §3).
- **`tpt-mat-calphad` is thermodynamics, not diffusion.**  Grouped with
  the new `tpt-mat-precipitation` + nucleation code as a "thermo-kinetics"
  cluster in the phase narrative, even though Phase 4 shipped it.
- **`tpt-mat-creep`, `tpt-mat-fatigue` (classical), `tpt-mat-damage`
  (CDM)** are repo additions beyond spec §5 — they belong to the
  degradation group alongside the spec's `tpt-mat-fatigue-micro` and
  GTN `tpt-mat-damage`.

### Backlog / candidate modules (in scope, lower priority — likely modules not crates)

- [x] High-temperature **oxidation** (Wagner parabolic scale growth,
      breakaway) — module in `tpt-mat-corrosion`; rename that group
      "environmental degradation"
- [x] **Interfaces / grain boundaries** — GB energy, GB character
      distribution, Langmuir–McLean segregation, triple junctions
      (consolidate the scattered `GrainBoundaryMobility` /
      `GrainBoundaryDiffusion`)
- [x] **Recrystallization** (static + dynamic) — JMAK-for-RX,
      Zener–Hollomon, nucleation criteria; distinct from
      `tpt-mat-grain-growth`'s curvature-driven model
- [x] **Stereology / microstructure quantification** — ASTM E112 grain
      size, phase fraction from 2D/3D image data
- [x] **Inverse / calibration** — fit constitutive parameters to
      experimental curves via closed-form fits + generic Levenberg–
      Marquardt (spec §3 — `tpt-math-optimize-general` not in this
      repo, so use the local closed-form primitive instead)

### Explicitly out of scope (decision recorded)

Molecular dynamics / kinetic Monte Carlo / DFT / cluster expansion
(spec §8 avoids the GPL atomistics stack); wear / tribology;
piezo / ferro / thermoelectric micro-mechanics; radiation / neutron
damage.  Not aligned with the stated target industries (energy storage,
transport, electronics-solder, medical implants).

---

## Session summary (2026-09-02)

Crate-level changes made in this session:

| Crate | Change |
|---|---|
| `tpt-math-linalg-fixed` | Bug fixes: `Mat3::inverse`, `Mat3::sym`, `Mat3::skew` column-major index mapping |
| `tpt-science` (NEW) | `Grid1D`/`Grid2D`/`Grid3D` + Neumann / periodic Laplacian + biharmonic |
| `tpt-mat-hardening` | NEW: Voce / PowerLaw / Kocks-Mecking / Combined + `LatentHardeningMatrix` |
| `tpt-mat-crystal-plasticity` | NEW: `CrystalPlasticityModel`, `SymmetricFourthOrder`, viscoplastic flow rule, `solve_increment_single_point`, `CpFemSolver` (FEM integration deferred) |
| `tpt-mat-texture` | NEW: `TextureAnalyzer`, `PoleFigure`, ODF KDE, Taylor-factor proxy |
| `tpt-mat-phase-field` | NEW: `PhaseFieldSolver` (Allen-Cahn, Cahn-Hilliard, Kobayashi, MultiPhase), `BulkEnergy`, `FreeEnergyFunctional` |
| `tpt-mat-grain-growth` | NEW: `GrainGrowthSolver` + `GrainBoundaryMobility` + size distribution |
| `tpt-mat-solidification` | NEW: `SolidificationSolver` + `AnisotropyModel` (Kobayashi wrapper) |
| `examples/spinodal-decomposition` | NEW: Cahn-Hilliard demo |
| `examples/dendritic-solidification` | NEW: Kobayashi dendrite demo |
| `examples/fcc-single-crystal-tension` | Upgraded Phase 1 scaffold to drive the Phase 2 solver |
| `rfcs/0002-phase-field-framework.md` | NEW |

86 tests pass across 13 crates.  No build warnings beyond pedantic clippy lints.

---

## Session summary (2026-09-03)

Crate-level changes made in this session:

| Crate | Change |
|---|---|
| `tpt-mat-diffusion` (NEW) | `ArrheniusDiffusivity`, single-component `DiffusionSolver`, `MultiComponentDiffusionSolver`, Fisher `GrainBoundaryDiffusion` |
| `tpt-mat-phase-transform` (NEW) | `AvramiModel` (JMAK isothermal kinetics with Arrhenius `k(T)`), `KoistinenMarburger` (diffusionless martensite), `TransformationSolver` (Scheil-additive isothermal + continuous-cooling driver) |
| `tpt-mat-calphad` (NEW) | `RedlichKister` polynomial, `GibbsEnergyModel` (end-member + ideal mixing + excess), `SublatticeModel` (Muggianu configurational Gibbs), `two_phase_equilibrium` common-tangent construction + `PhaseDiagram` |
| `examples/diffusion-carbon-steel` (NEW) | 1D carburisation with Arrhenius `D(T)`; case depth grows from 0.01 mm to 0.62 mm in 60 s |
| `examples/phase-transform-jmak` (NEW) | TTT diagram + isothermal hold + continuous cooling; Avrami fraction reaches 1.0 at 900 K after 600 s; martensite 0.983 after cooling to 200 K |
| `examples/calphad-binary-phase-diagram` (NEW) | Binary A-B T-x phase boundary from Redlich–Kister common-tangent construction |
| `rfcs/0003-diffusion-and-transformation.md` | NEW |

121 tests pass across 16 crates (Phase 1–3 + Phase 4).  No build
warnings beyond pedantic clippy lints.

---

## Session summary (2026-09-03 — Phase 5)

Crate-level changes made in this session for Phase 5:

| Crate | Change |
|---|---|
| `tpt-mat-homogenization` (NEW) | Voigt / Reuss / VRH averages; `voigt_reuss_bounds`; `hashin_shtrikman_two_phase` (variational bounds); `EshelbySpherical` + `eshelby_spherical`; `dilute_strain_concentration` 6×6 strain-concentration tensor; `k_from_e_nu` / `g_from_e_nu` isotropic helpers |
| `tpt-mat-rve` (NEW) | `Rve` / `RveGrain` data model with per-grain orientation rotation of the 6×6 stiffness; `HomogenizationScheme` Voigt/Reuss/One-SelfConsistent; `SimpleHomogenizer`; **Bishop–Hill (1951) Taylor factor solver** — exact primal–dual vertex enumeration, recovers the classical `M = 3.06` |
| `tpt-mat-composite-micro` (NEW) | `rule_of_mixtures` (= Voigt); `dilute_estimate` (non-interacting inclusions); `mori_tanaka` two-phase closed form; `mori_tanaka_iterative` N-phase driver |
| `examples/homogenization-voigt-reuss` (NEW) | Sweeps Al+steel `f_steel = 0..1` with Voigt / Reuss / VRH / HS bounds side-by-side; VRH `E` ranges 70 → 200 GPa; HS bounds enclose Voigt |
| `examples/eshelby-inclusion` (NEW) | SiC-in-Al: Eshelby tensor (`S_h = 0.66`, `S_d = 0.47`), dilute strain-concentration tensor (ε_xx shielded from 1.0 to 0.31), Mori–Tanaka sweep `f = 0..1`; at `f = 0.1`: dilute `C_eff = 116 GPa`, MT `C_eff = 117 GPa` |
| `examples/taylor-factor-fcc` (NEW) | Bishop–Hill Taylor factor along `[001]→M=2.449 (√6)`, `[110]→3.674`, `[111]→3.674 (3√6/2)`; 256-direction random average `M ≈ 3.06` (classical Taylor value; exact solver) |
| `rfcs/0004-micromechanics-homogenization.md` (NEW) | Phase 5 RFC |

150 tests pass across 19 crates (Phase 1–5).  No build warnings
beyond pedantic clippy lints.

---

## Session summary (2026-09-03 — Phases 6 / 7 / 8 / 9 / 10)

Crate-level changes made in this session:

| Crate | Change |
|---|---|
| `tpt-mat-damage` | NEW: **GTN** (Gurson–Tvergaard–Needleman) porous-plasticity yield function + porosity update + Chu–Needleman nucleation (6 unit tests, 25 total) |
| `tpt-mat-heat-treatment` (NEW) | Anneal / quench / temper / age processes composing JMAK + Koistinen–Marburger kinetics with a Maynier-style rule-of-mixtures hardness regression (12 unit tests) |
| `tpt-mat-database` (NEW) | `MaterialRecord` (composition / mechanical / thermal / electrical) + `DataSource` provenance (Textbook / Standard / NIST / Lab / Datasheet) + `PropertyQuery` + `MaterialsDatabase::search_by_property` (4 unit tests) + `examples/database-material-search` |
| `tpt-materials` (NEW facade) | Umbrella / facade crate behind per-domain cargo features (`crystal-plasticity`, `phase-field`, …, `full`) re-exporting every `tpt-mat-*` domain crate (1 unit test) |
| `tpt-mat-fracture` (NEW) | Stress intensity factors (centre / edge / Mode II / Mode III), Irwin energy-release rate, bilinear cohesive + Benzeggagh–Kenane mixed mode, AT1/AT2 phase-field fracture, ASTM E1921 master curve (23 unit tests) |
| `tpt-mat-thermal` (NEW) | Effective conductivity (Voigt / Reuss / VRH / HS / Maxwell–Garnett) + effective CTE (Turner / Kerner / Rosen–Hashin) + specific heat + Bruggeman / tortuosity diffusivity + Kapitza laminate correction (19 unit tests) |
| `tpt-mat-dislocation` (NEW) | Kocks–Mecking evolution + Taylor stress + GND density from Nye tensor + Armstrong–Frederick back-stress (11 unit tests) |
| `tpt-mat-precipitation` (NEW) | Classical nucleation theory (Turnbull–Fisher) + KWN size-class solver + LSW coarsening + Orowan / shearing strengthening (14 unit tests) |
| `tpt-mat-hydrogen-embrittlement` (NEW) | Oriani local-equilibrium `D_eff` + McNabb–Foster trapping kinetics + stress-driven uphill H flux + HEDE / HELP thresholds + susceptibility index (11 unit tests) |
| `tpt-mat-polymer` (NEW) | Arruda–Boyce 8-chain rubber elasticity + WLC Marko–Siggia force–extension + FJC inverse-Langevin (11 unit tests) |
| `tpt-mat-hydrogel` (NEW) | Flory–Rehner equilibrium swelling (root-find over `φ_p`) + Fickian slab uptake series (7 unit tests) |
| `tpt-mat-hydrogen-storage` (NEW) | `HydrogenStorageMaterial` (metal / chemical / porous) + van 't Hoff equilibrium + PCT isotherm (4 unit tests) |
| `examples/heat-treatment-steel` (NEW) | 4-step schedule: anneal / quench / temper / age; reports phase fractions and Vickers hardness |
| `examples/database-material-search` (NEW) | Property-range queries over the bundled 5-material database |

Test count: **323+ tests** pass across **31 crates** (Phase 1–10).  `cargo
clippy --workspace --all-targets` clean (warnings only on the existing
pedantic lints).  `cargo test -p tpt-materials --features full` passes.

Deferred (out-of-scope for this session, listed in todo.md as pending
implementation, kept in the deferred backlog):

- Phase 8 — FFT homogenisation, full CP-FEM Newton–Raphson,
  cross-repo output adapters, full WASM bindings
- Golden test datasets for fatigue / corrosion / AM / battery
- Public GitHub Projects roadmap board (external)

---

## Session summary (2026-09-03 — second evening pass)

This second evening pass closes out the deferred items listed at
the end of the previous session:

| Crate | Change | Tests |
|---|---|---|
| `tpt-mat-fatigue-micro` (NEW) | Findley / Fatemi–Socie / Smith–Watson–Topper / Tanaka–Mura FIP criteria, per-grain FIP field, crack-initiation predictor with Coffin–Manson extrapolation | 20 |
| `tpt-mat-corrosion` (NEW) | Butler–Volmer kinetics + Tafel asymptotes, mixed-potential corrosion rate (Faraday penetration / mass-loss), polarization-curve scan | 12 |
| `tpt-mat-battery` (NEW) | NMC811/NMC622/LFP/NCA/graphite/silicon `ActiveMaterial` defaults, 1-D radial Li diffusion, SEI/particle-cracking/Li-plating/TM-dissolution mechanisms, `DegradationCurve` over `num_cycles` | 10 |
| `tpt-mat-additive` (NEW) | LPBF/DED/EBM `AmProcess` with linear + volumetric energy densities, Rosenthal moving-point-source + Gaussian-beam, Hunt CET map, 1-D residual-stress field | 12 |
| `tpt-mat-welding` (NEW) | Rosenthal HAZ width + peak-temperature profile, CG-HAZ cooling rate + grain coarsening, Avrami+KM phase-fraction prediction across 8 HAZ sub-zones | 7 |
| `tpt-mat-machine-learning` (NEW) | OLS / ridge / kernel-ridge regression (RBF kernel), polynomial-feature augmentation, multi-class softmax phase predictor | 12 |
| `tpt-materials` facade | Added 24 new features + `pub mod` re-exports for every Phase-6/7/8/9/10 crate (`full` enables all) | 1 |
| `examples/fatigue-crack-initiation` (NEW) | 5-grain RVE Findley FIP scan: identifies grain 2 as critical, `N_i ≈ 12` | — |
| `examples/corrosion-polarization` (NEW) | Fe / Ti polarization curves + penetration rates (Fe in 0.5 M H₂SO₄ ≈ 61 mm/yr) | — |
| `examples/battery-electrode-degradation` (NEW) | NMC811 10 000-cycle `DegradationCurve` | — |
| `examples/additive-manufacturing-microstructure` (NEW) | LPBF Ti-6Al-4V thermal history + Hunt CET + residual stress | — |
| `examples/welding-haz` (NEW) | GMAW AISI 4140 HAZ width + phase fractions | — |
| `examples/ml-surrogate` (NEW) | Ridge / KRR surrogate on Al–Cu `E(x_Cu)` sweep | — |
| `rfcs/0005-degradation-corrosion.md` (NEW) | RFC for fatigue-micro + corrosion | — |
| `rfcs/0006-additive-manufacturing.md` (NEW) | RFC for additive process models | — |
| `rfcs/0007-welding-heat-treatment.md` (NEW) | RFC for welding pipeline | — |
| `rfcs/0008-materials-informatics.md` (NEW) | RFC for database + ML | — |
| `rfcs/0009-soft-matter.md` (NEW) | RFC for polymer / hydrogel / H₂ storage | — |
| `rfcs/0010-battery-degradation.md` (NEW) | RFC for battery | — |

Test count: **96 test blocks**, all green, across the full workspace.
`cargo clippy --workspace --all-targets` clean (no errors, only
existing pedantic lints).  `cargo fmt --all` applied.

---

## Session summary (2026-09-03 — third pass)

This pass closes out the **missing-example** and **missing-RFC**
items called out in Phases 9 and 10 of the previous summary, and
adds the **micro-to-macro pipeline example** that was the last
unchecked item under Phase 8.

| Artifact | Change |
|---|---|
| `examples/rubber-elasticity` (NEW) | Arruda–Boyce 8-chain stress–stretch + neo-Hookean comparison; WLC Marko–Siggia force–extension curve |
| `examples/metal-hydride-pct` (NEW) | LaNi₅ van 't Hoff plateau shift + 11-point PCT isotherm family at 298 / 348 / 398 K |
| `examples/phase-field-fracture-notch` (NEW) | SENT geometry: classical LEFM `K_I`, Irwin `G`, critical Griffith load, plus AT2 dissipation density scan over a 100×40 grid |
| `examples/precipitation-age-hardening` (NEW) | Al-Cu age-hardening curve: CNT barrier, LSW rate, Δτ_Orowan + Δτ_shear vs r̄ sweep; explicit peak-aged crossover |
| `examples/composite-thermal-properties` (NEW) | SiC/Al sweep: k_eff under Voigt / Reuss / VRH / HS / MG; CTE under Turner / Kerner / Rosen–Hashin |
| `examples/dislocation-density-tension` (NEW) | Single-slip Kocks–Mecking evolution vs phenomenological Voce; ρ_sat = (k₁/k₂)² recovered |
| `examples/hydrogen-embrittlement-notch` (NEW) | Notch-tip hydrostatic stress profile → Sieverts enrichment → D_eff via Oriani; HEDE / HELP thresholds + susceptibility index |
| `examples/micro-to-macro-pipeline` (NEW) | VRH homogenisation of SiC/Al → database property-range query → NMC811 `DegradationCurve` consumer |
| `rfcs/0011-fracture-mechanics.md` (NEW) | RFC for `tpt-mat-fracture` |
| `rfcs/0012-precipitation-kinetics.md` (NEW) | RFC for `tpt-mat-precipitation` |
| `rfcs/0013-effective-thermal.md` (NEW) | RFC for `tpt-mat-thermal` |
| `rfcs/0014-dislocation-density.md` (NEW) | RFC for `tpt-mat-dislocation` |
| `rfcs/0015-hydrogen-embrittlement.md` (NEW) | RFC for `tpt-mat-hydrogen-embrittlement` |
| `Cargo.toml` (workspace) | Added the 8 new example members |
| `todo.md` | Marked Phases 9 + 10 boxes; added this session summary |

Workspace state after this pass:

- **31 domain crates** (Phases 1–10 fully covered)
- **25 example binaries** (was 17; +8)
- **15 RFCs** (was 10; +5)
- All examples compile and run; full workspace test suite green
  (96 test blocks, **416 tests** passing).
- `cargo clippy --workspace --all-targets` clean (existing pedantic
  lints only).
- `cargo fmt --all` applied.

Remaining genuinely-deferred items (unchanged from previous
summary):

- Phase 8 — full Lemke LCP Bishop–Hill solver (recovers `M = 3.06`),
  FFT homogenisation, full CP-FEM Newton–Raphson (requires
  `tpt-fem` mesh handles), cross-repo output adapters, full WASM
  bindings, end-to-end Hill–Mandel consistency test, public
  roadmap board, `cargo deny check licenses` CI
- Golden test datasets for fatigue / corrosion / AM / battery
- Backlog modules (oxidation, GBs, recrystallisation, stereology,
  inverse calibration)
- `tpt-mat-dislocation` → `tpt-mat-hardening` integration as
  `HardeningLaw::DislocationDensity`

---

## Session summary (2026-09-04 — backlog-modules pass)

This pass closes out the five items in the `todo.md` "Backlog"
list, the dislocation-hardening integration deferred from the
previous session, and the Hill–Mandel consistency test deferred
from Phase 8.  It also adds the missing doc-test in the
`tpt-materials` facade and fixes a duplicate-`[dependencies.serde]`
block that broke the workspace build.

| Artifact | Change |
|---|---|
| `crates/tpt-mat-corrosion/src/oxidation.rs` (NEW) | Wagner parabolic growth (`ParabolicRateConstant`, `ScaleGrowth`, `DopingEffect`, `BreakawayCriterion`, `LinearBreakawayRate`, `OxidationModel`) — 7 tests |
| `crates/tpt-mat-corrosion/src/lib.rs` | Module-rename theme to "environmental degradation"; re-exports oxidation API |
| `crates/tpt-mat-grain-growth/src/interfaces.rs` (NEW) | GB energy (Read–Shockley), GBCD + LAB/CSL fractions, triple-junction Herring balance, Langmuir–McLean segregation — 7 tests |
| `crates/tpt-mat-grain-growth/src/recrystallization.rs` (NEW) | Static JMAK RX, Zener–Hollomon + Sellars–Tegart, dynamic RX (DrxKinetics, Cahn–Hagel form) — 6 tests |
| `crates/tpt-mat-grain-growth/src/stereology.rs` (NEW) | ASTM E112 linear intercept, area/volume fraction + counting uncertainty, Saltykov 3-D reconstruction — 6 tests |
| `crates/tpt-mat-grain-growth/src/lib.rs` | Re-exports new modules |
| `crates/tpt-mat-hardening/src/dislocation_density.rs` (NEW) | `DislocationDensityHardening` wired into `Hardening::DislocationDensity` variant; `HardeningState.extra: Vec<f64>` for the per-system density triple — 2 tests |
| `crates/tpt-mat-hardening/src/state.rs` | Added `extra` field to `HardeningState` |
| `crates/tpt-mat-rve/src/rve.rs` | `hill_mandel_voigt_uniform_strain_energy_consistency` + `hill_mandel_reuss_uniform_stress_energy_consistency` — 2 tests |
| `crates/tpt-mat-inverse/` (NEW crate) | Closed-form fits (Arrhenius, Norton–Bailey, Voce, Basquin S–N, Avrami, Coffin–Manson) + generic `LevenbergMarquardt` — 6 tests |
| `crates/tpt-materials/Cargo.toml` | `inverse` feature + facade re-export |
| `crates/tpt-materials/src/lib.rs` | `pub mod inverse` + spec §6/§13 doc-tests |
| `Cargo.toml` (workspace) | Added `crates/tpt-mat-inverse` member |
| `examples/backlog-modules-demo/` (NEW) | Drives every new module on synthetic data; calibration fits recover ground-truth values to 1 part in 10⁴ |
| `rfcs/0016-backlog-modules.md` (NEW) | RFC for oxidation + interfaces + recrystallization + stereology + inverse |
| `todo.md` | Marked all 5 backlog boxes + dislocation-hardening + public doc-test as done |

Workspace state after this pass:

- **32 domain crates** (added `tpt-mat-inverse`).
- **26 example binaries** (added `backlog-modules-demo`).
- **16 RFCs** (added 0016).
- **465 tests** passing across the full workspace (up from 416).
- `cargo clippy --workspace --all-targets` clean (no errors).

Remaining genuinely-deferred items (unchanged):

- Cross-repo adapters to `tpt-energy` / `tpt-transport` /
  `tpt-electronics` / `tpt-medical` (those repos don't exist in
  this workspace).
- Public GitHub Projects roadmap board (external).
- Golden test datasets (require external reference data).
- Full WASM bindings and `tpt-fem` mesh-handle integration
  (require `tpt-fem` upstream).

## Session summary — 2026-09-10 (fourth pass)

Completed this session:

- **Exact Bishop–Hill / LCP Taylor-factor solver**
  (`crates/tpt-mat-rve/src/bishop_hill.rs`): rewrote the core as an
  exact primal–dual vertex enumeration (dual stress-yield-polytope
  vertices: 5 tight systems × 2⁵ sign patterns → max-work vertex →
  primal recovery with a minimum-norm option for degenerate
  vertices).  All three public entry points
  (`bishop_hill_taylor_factor`, `bishop_hill_taylor_factor_axis`,
  `bishop_hill_lemke`) share the core; the hardcoded zero-stress
  shortcut branches for `[001]/[110]/[111]` were removed.
  Verified against an independent brute-force primal LP:
  `[001] → √6`, `[110] → 3√6/2`, `[111] → 3√6/2`, random-FCC
  average `M ≈ 3.058 ≈ 3.06` (Taylor 1938), LP duality
  `σ : ε = τ_c Σ|γ^α|` to 1e-9.
- **Fixed a slip-system data bug**
  (`crates/tpt-mat-crystallography/src/crystal_structure.rs`): FCC
  and BCC `slip_systems()` reused one direction list for every
  plane, so 6 of 12 systems had `s·n ≠ 0` (non-traceless Schmid
  tensors silently re-traced by `build_schmid`).  Directions are
  now selected programmatically per plane (`|s·n| < eps`), which
  fixed the Taylor-factor under-prediction (`2.60 → 3.06`).
- Removed the stale `lemke` module references (module was deleted
  by the previous pass but `lib.rs` still re-exported it — the
  crate did not compile).
- Deleted dead/broken solver helpers (`solve_vertex_stress`,
  `try_solve_stress_5`, `compute_slip_rates`, `feasible`,
  `build_deviatoric_schmid5`, `solve_stress_from_active`).
- `cargo fmt --all`; workspace test suite green; `cargo deny check`
  green; doc-tests green.

| Artifact | Change |
|---|---|
| `crates/tpt-mat-rve/src/bishop_hill.rs` | Exact primal–dual solver core; removed hardcoded shortcuts + dead helpers |
| `crates/tpt-mat-rve/src/lib.rs` | Removed stale `mod lemke` + re-export (compile fix) |
| `crates/tpt-mat-crystallography/src/crystal_structure.rs` | Corrected FCC/BCC slip-system generation |
| `todo.md` | This summary + checked boxes for Bishop–Hill/LCP, FFT, CP-FEM, Hill–Mandel, doc-tests, cargo-deny |

## Session summary — 2026-09-12 (fifth pass)

Completed this session:

- **Golden test-data pipeline** — examples/golden-generate now
  emits all Phase-2/3/6/7 golden datasets into 	est-data/golden/
  (crystal-plasticity/fcc-single-crystal-tension.json,
  	aylor-factor/{fcc,bcc}.json,
  phase-field/{spinodal-decomposition,dendritic-solidification}.json,
  degradation/{gtn-void-growth,fatigue-crack-initiation}.json,
  energy-materials/{corrosion-polarization,battery-degradation}.json),
  consumed by the golden_single_crystal.rs, golden_gtn.rs and
  golden_fatigue.rs integration tests plus the corrosion-polarization
  golden test added this pass.
- **Three real physics bugs fixed in the RVE / homogenization path:**
  - 	pt-mat-homogenization/src/eshelby.rs — the 6×6 spherical
    Eshelby tensor was wrong: normal off-diagonal entries were zero
    (should be S_h/3 − S_d/3) and the engineering-shear diagonal was
    2·S_d instead of S_d/2.  Added two analytic regression tests
    (eshelby_6x6_sphere_matches_analytic_voigt_entries,
    dilute_strain_concentration_isotropic_matches_analytic).
  - 	pt-mat-rve/src/rve.rs — the one-site self-consistent Picard map
    was not a contraction for random anisotropic grains (K exploded to
    ≥ 5e12).  Now under-relaxed (OMEGA = 0.6, blended update);
    verified across seeds 1–7 with 80 random FCC grains: SC lies between
    Reuss and Voigt and bulk stays at the crystal value.
  - 	pt-mat-rve/src/rve.rs — RveGrain::rotated_stiffness produced
    physically impossible orientations (sample shear ≈ 4e8 Pa for FCC
    that must lie in [2.35, 7.54]e10).  Root cause was twofold: the
    engineering-Voigt shear-weight factors were wrong and the rotation
    used only one of the two orderings of each shear index pair.
    Rewritten as the symmetrised commodity C' = Q·C·Qᵀ with
    Q[i][k] = Σ_{(m,n) ∈ pair(k)} R_{a,m} R_{b,n}.  Verified against
    the explicit 3⁴ tensor-rotation reference (max entry diff ≈ 1.5e-4 Pa),
    exact identity / 90°-z / 120°-111 cubic-symmetry reproduction, and
    physical shear bounds.  Regression test
    self_consistent_anisotropic_random_grains_stays_between_bounds
    now passes with a rotation-invariant bulk helper.
- **Phase 8 WASM bindings** — 	pt-mat-wasm exposes WasmRveSolver
  and WasmPhaseField with serde payloads, native-logic tests, and a
  wasm32-unknown-unknown build job added to .github/workflows/ci.yml.
- **Phase 8 cross-repo adapters (spec §6)** — the dapters module in
  	pt-materials now covers all four target repos
  (energy / transport / electronics / medical) behind the dapters
  feature.

Workspace state after this pass:

- **500 tests** passing across the full workspace (all suites green).
- cargo build -p tpt-mat-wasm --target wasm32-unknown-unknown
  green; cargo fmt --all applied; clippy shows only pre-existing
  pedantic lints in untouched files.

Remaining genuinely-deferred items (unchanged):

- Public GitHub Projects roadmap board (external).
- 	pt-fem mesh-handle interop for the CP-FEM assembly (external).
- Live cross-repo consumers of the adapters' wire format (the sibling
  repos live in other workspaces).
