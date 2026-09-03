# tpt-materials — Project Todo

Organization: **TPT Solutions** · License: **MIT OR Apache-2.0** (dual)

> **Status snapshot — 2026-09-03.**  Phases 1–3 from the original
> scaffold are complete; Phases 1 (foundation), 2 (crystal
> plasticity) and 3 (phase-field) were implemented in the
> 2026-09-02 session, and Phase 4 (diffusion & transformation) was
> implemented earlier in this session with library code, unit tests,
> examples, and RFC 0003.  **Phase 5 (micro-mechanics) is
> implemented in this session** with three new crates
> (`tpt-mat-homogenization`, `tpt-mat-rve`, `tpt-mat-composite-micro`),
> unit tests, examples, and RFC 0004.  **Phases 6–8 are still
> pending implementation but are now fully task-broken-down from
> `spec.txt`** (this session), and a **Phase 9** was added to cover
> the soft-matter / hydrogen-storage crates the spec defines but
> never schedules.  Phase 6 also folds in three classical
> analytical crates already drafted this session (`tpt-mat-creep`,
> `tpt-mat-damage`, `tpt-mat-fatigue`) that are not yet committed.
> A **crate-organisation review** (2026-09-03) then added **Phase 10**
> (fracture, precipitation, effective thermal properties,
> dislocation-density plasticity, hydrogen embrittlement — standard
> topics no spec crate covered), a `tpt-materials` facade-crate task
> under Phase 8, and a "Workspace organisation notes" section
> (flat-layout rationale + backlog + out-of-scope decisions).

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
  - [ ] **Full FEM assembly + Newton-Raphson**: deferred.  The
    constitutive model kernel is complete; FEM assembly is wired
    through optional `tpt-fem-solve` / `tpt-fem-assembly`
    dependencies (off by default) — they are not yet exercised by
    an integration test because the `tpt-fem` mesh handle types are
    not in this workspace.
- [x] `tpt-mat-texture`
  - [x] `TextureAnalyzer` (orientations + weights; `from_ebsd`)
  - [x] `PoleFigure` (equal-area / stereographic; `Fcc111`,
        `Fcc200`, `Bcc110`, `Hcp0001`, `Hcp10T10`, custom)
  - [x] `OrientationDistributionFunction` (geodesic-Gaussian KDE on SO(3))
  - [x] `taylor_factor`: **single-Schmid proxy** — see caveat below.
- [x] RFC 0001: `rfcs/0001-crystal-plasticity-fem.md` (existing)
- [ ] Golden test data: FCC single-crystal tension / BCC polycrystal
      RVE / Taylor-factor random textures — *not produced*; they
      require FEM-integration test results that do not exist yet.
- [ ] **Verification test: FCC random-texture Taylor factor ≈ 3.06**
      is **not** satisfied by the current proxy (M ≈ 2.0–2.3).
      Documented as a limitation: the proper Taylor factor requires
      a Bishop–Hill LCP solver to find the 5 slip systems that
      accommodate the macroscopic strain; the single-Schmid
      `1 / max_schmid` proxy over-counts because it assumes
      single-system activation.  The full solver is queued for
      Phase 5 alongside the RVE/homogenization stack.
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
- [ ] Golden test data: not produced (would require real
      spinodal-dendrite numerical reference data, which is
      outside scope).
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
  - [x] **Bishop–Hill (1951) Taylor factor solver**
        (`bishop_hill_taylor_factor_axis` /
        `bishop_hill_taylor_factor`)
        — replaces the Phase 2 single-Schmid proxy with a true
        // primal solver; documented limitation: uses L2
        // pseudo-inverse, recovers `M ≈ 2.0–2.5` for random FCC
        // (vs Taylor's classical `3.06`); L1 solver deferred.
  - [ ] **Full LCP (Lemke) Bishop–Hill solver**: queued for the
        Phase 8 ecosystem / informatics work to recover `M = 3.06`.
  - [ ] **FFT homogenisation (Moulinec–Suquet 1998)**: queued
        for Phase 8.
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
- [ ] Register `tpt-mat-creep` / `tpt-mat-damage` / `tpt-mat-fatigue`
      in workspace and commit; add per-crate RFC coverage

### Spec items still to build

- [ ] `tpt-mat-damage` — Gurson–Tvergaard–Needleman
  - [ ] `GursonTvergaardNeedleman` params (`f_0`, `f_c`, `f_f`,
        `q_1`, `q_2`, `q_3`, `f_n`, `s_n`, `e_n`)
  - [ ] `yield_function(stress, porosity, matrix_yield)`
        (`Φ = (σ_eq/σ_y)² + 2q₁f cosh(3q₂σ_m/2σ_y) − (1 + q₃f²)`)
  - [ ] `update_porosity` (void growth `(1−f)dε^p_kk` + strain-controlled
        nucleation `A dε^p_eq`)
  - [ ] Verification test: GTN yield surface → von Mises as `f → 0`
- [ ] `tpt-mat-fatigue-micro`
  - [ ] `MicrostructuralFatigue { rve, criterion }`
  - [ ] `FatigueCriterion`: `Findley`, `FatemiSocie`,
        `SmithWatsonTopper`, `CrystallographicSlip { critical_accumulated_shear }`
  - [ ] `fatigue_indicator_parameter(&CpFemResult) -> Vec<f64>` (FIP field)
  - [ ] `predict_crack_initiation(&[LoadStep]) -> CrackInitiationResult`
        (`cycles_to_initiation`, `critical_grain`, `critical_location`,
        `fip_field`)
  - [ ] Cycle-by-cycle CP-FEM driver accumulating plastic slip at
        grain boundaries
- [ ] `tpt-mat-corrosion`
  - [ ] `ElectrodeKinetics { exchange_current_density, tafel_slope,
        equilibrium_potential }`
  - [ ] `CorrosionModel { anode, cathode, electrolyte }`
  - [ ] `corrosion_rate(temperature, ph) -> CorrosionRate`
        (Butler–Volmer mixed-potential solve; `current_density`,
        `penetration_rate` mm/yr, `mass_loss_rate` g/m²·day)
  - [ ] `polarization_curve(potential_range) -> PolarizationCurve`
        (anodic / cathodic Tafel branches)
- [ ] RFC 0005: degradation & failure models (GTN + FIP + Butler–Volmer)
- [ ] Golden test data: `test-data/golden/degradation/`
      (`gtn-void-growth.json`, `fatigue-crack-initiation.json`,
      `corrosion-polarization.json`)
- [ ] Example: `examples/fatigue-crack-initiation` (polycrystal RVE →
      FIP field → critical grain)
- [ ] Example: `examples/corrosion-polarization` (Ti-6Al-4V
      polarization curve + penetration rate)

**Milestone:** Predict fatigue crack initiation site in a polycrystal
(critical grain + cycles-to-initiation from the FIP field).

---

## Phase 7 — Energy Materials & Manufacturing (Months 19-21)

**Crates (spec §5 Domains 8-9):** `tpt-mat-battery`, `tpt-mat-additive`,
`tpt-mat-welding`, `tpt-mat-heat-treatment`
**Substrate:** `tpt-mat-diffusion`, `tpt-mat-phase-transform`,
`tpt-mat-calphad`, `tpt-mat-solidification`, `tpt-mat-grain-growth`

- [ ] `tpt-mat-battery`
  - [ ] `ActiveMaterial { chemistry, particle_radius,
        diffusion_coefficient, partial_molar_volume }`
  - [ ] `BatteryChemistry` (NMC811, LFP, NCA, graphite, silicon, …)
  - [ ] `DegradationMechanism`: `SeiGrowth { rate_constant,
        activation_energy }`, `ParticleCracking { critical_stress }`,
        `LithiumPlating { plating_potential }`,
        `TransitionMetalDissolution { dissolution_rate }`
  - [ ] `simulate_diffusion_stress(c_rate, num_cycles)` — coupled
        Li diffusion + diffusion-induced stress in a spherical particle
        (reuses `tpt-mat-diffusion`)
  - [ ] `capacity_fade_curve(cycles, temperature) -> DegradationCurve`
        (`cycles`, `capacity_retention`, `resistance_growth`)
  - [ ] Verification test: capacity retention monotonically decreasing;
        √t SEI-limited fade at low C-rate
- [ ] `tpt-mat-additive`
  - [ ] `AmProcess`: `LaserPowderBedFusion { laser_power, scan_speed,
        hatch_spacing, layer_thickness }`, `DirectedEnergyDeposition`,
        `ElectronBeamMelting`
  - [ ] `thermal_history(location) -> ThermalHistory`
        (Rosenthal / moving-source analytic solution)
  - [ ] `predict_microstructure(&ThermalHistory) -> PredictedMicrostructure`
        (`grain_size`, `phase_fractions`, `texture`, `porosity`;
        columnar/equiaxed from G–R solidification map)
  - [ ] `residual_stress(&ThermalHistory) -> ResidualStressField`
        (thermal-contraction eigenstrain)
  - [ ] Verification test: higher cooling rate → finer predicted grain
        size (Hall–Petch trend)
- [ ] `tpt-mat-welding`
  - [ ] `WeldModel { base_metal, filler_metal, process }`
  - [ ] `heat_affected_zone(heat_input) -> HazResult` (HAZ width +
        peak-temperature profile)
  - [ ] `predict_haz_microstructure(cooling_rate) -> HazMicrostructure`
        (grain coarsening + transformation via `tpt-mat-phase-transform`)
- [ ] `tpt-mat-heat-treatment`
  - [ ] `HeatTreatmentProcess`: `Annealing`, `Quenching { medium }`,
        `Tempering`, `Aging`, `SolutionTreatment`
  - [ ] `simulate(&HeatTreatmentProcess) -> HeatTreatmentResult`
        (final phase fractions + grain size, CALPHAD + TTT driven)
  - [ ] `predict_hardness(&PredictedMicrostructure) -> f64`
        (rule-of-mixtures / Maynier-type regression)
- [ ] RFC 0006: additive-manufacturing microstructure & residual stress
- [ ] RFC 0007: heat-treatment / welding transformation pipeline
- [ ] Golden test data: `test-data/golden/energy-materials/`
      (`battery-sei-growth.json`, `electrode-capacity-fade.json`)
- [ ] Example: `examples/battery-electrode-degradation` (NMC811
      capacity fade curve — the artifact exported to `tpt-energy`)
- [ ] Example: `examples/additive-manufacturing-microstructure`
      (LPBF thermal history → grain structure + residual stress)

**Milestone:** Generate a battery capacity fade curve consumable by
`tpt-energy` (`DegradationCurve` over ≥10 000 cycles).

---

## Phase 8 — Informatics & Ecosystem (Months 22-24)

**Crates (spec §5 Domain 10):** `tpt-mat-database`, `tpt-mat-machine-learning`
**Plus:** cross-repo integration (spec §6), WASM bindings (spec §7),
deferred solver upgrades from Phases 2 & 5

- [ ] `tpt-mat-database`
  - [ ] `MaterialRecord { name, composition, mechanical, thermal,
        electrical, sources }`
  - [ ] `MaterialsDatabase { materials }` + `load_builtin()`
        (bundled MIT-clean property set)
  - [ ] `search_by_property(PropertyQuery) -> Vec<MaterialRecord>`
        (property-range queries)
  - [ ] `DataSource` provenance (ASTM / ISO / NIST traceability, spec §9)
- [ ] `tpt-mat-machine-learning`
  - [ ] `MlModel`: `PropertyPredictor { features }`, `PhasePredictor`,
        `SurrogateModel`
  - [ ] `train(&[(Vec<f64>, f64)]) -> TrainingResult`
  - [ ] `predict(&Composition) -> f64`
  - [ ] Verification test: surrogate reproduces a CP-FEM / homogenization
        sweep within tolerance
- [ ] **Full LCP (Lemke) Bishop–Hill Taylor-factor solver** (carried
      from Phase 5) — recover `M ≈ 3.06` for random FCC
- [ ] **FFT homogenisation (Moulinec–Suquet 1998)** (carried from
      Phase 5) — spectral full-field RVE in `tpt-mat-homogenization`
- [ ] **Full CP-FEM assembly + Newton–Raphson** (carried from Phase 2)
      once `tpt-fem` mesh handles are available
- [ ] Cross-repo output adapters (spec §6):
  - [ ] `tpt-energy` ← battery `DegradationCurve`
  - [ ] `tpt-transport` ← composite fatigue S–N + alloy creep (`tpt-mat-creep`)
  - [ ] `tpt-electronics` ← solder-joint fatigue / reliability
  - [ ] `tpt-medical` ← implant corrosion rate / biocompatibility
- [ ] WASM (spec §7): `tpt-mat-wasm` — `WasmRveSolver`,
      `WasmPhaseField` (`step`, `get_order_parameter`), build in
      `benchmark`/`docs` CI
- [ ] `tpt-materials` umbrella / facade crate — thin re-export so the
      spec §6 & §13 snippets (`use tpt_materials::crystallography::…`,
      `::crystal_plasticity::…`, `::energy_materials::…`) resolve
  - [ ] One `pub mod` per spec domain, re-exporting the domain crates
  - [ ] Per-domain cargo features (`crystal-plasticity`, `phase-field`,
        …); `full` enables all; `wasm` pulls `tpt-mat-wasm`
  - [ ] Doc-test the exact import snippets from spec §6 and §13
- [ ] RFC 0008: materials-informatics database schema + ML surrogates
- [ ] Verification test: end-to-end Hill–Mandel consistency through
      the full micro→macro pipeline
- [ ] Example: `examples/micro-to-macro-pipeline` (microstructure →
      homogenized property → device-level input)

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

- [ ] `tpt-mat-polymer`
  - [ ] `ChainModel`: `FreelyJointedChain { num_segments,
        segment_length }`, `WormLikeChain { persistence_length,
        contour_length }`, `ArrudaBoyce { n_segments, shear_modulus }`
  - [ ] `PolymerModel { chain_model, crosslink_density }`
  - [ ] `stress_strain(stretch) -> f64` (Arruda–Boyce 8-chain via
        inverse Langevin; WLC force–extension)
  - [ ] Verification test: Arruda–Boyce → neo-Hookean at small stretch
- [ ] `tpt-mat-hydrogel`
  - [ ] Flory–Rehner swelling equilibrium (mixing + elastic osmotic
        pressure balance)
  - [ ] Poroelastic swelling kinetics (Fickian solvent uptake on a
        `tpt-science` grid)
  - [ ] `equilibrium_swelling_ratio(chi, crosslink_density) -> f64`
- [ ] `tpt-mat-hydrogen-storage`
  - [ ] `HydrideType`: `MetalHydride { alloy }`,
        `ChemicalHydride { compound }`, `PorousMaterial { surface_area }`
  - [ ] `HydrogenStorageMaterial { hydride_type, storage_capacity_wt_pct,
        absorption_kinetics }`
  - [ ] `pct_isotherm(temperature) -> Vec<(f64, f64)>` (pressure–
        composition–temperature curve with plateau + van 't Hoff
        temperature dependence)
- [ ] RFC 0009: soft-matter constitutive models
- [ ] RFC 0010: hydrogen-storage sorption kinetics
- [ ] Example: `examples/rubber-elasticity` (Arruda–Boyce uniaxial)
- [ ] Example: `examples/metal-hydride-pct` (LaNi₅ PCT isotherm family)

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

- [ ] `StressIntensityFactor` — `K_I` / `K_II` / `K_III`, geometry
      factors (edge / centre / penny-shaped crack),
      `k_from_load(geometry, stress, crack_length)`
- [ ] `EnergyReleaseRate` — `G`, Irwin `G = K² / E'`, J-integral
      (domain-integral form on a `tpt-science` grid)
- [ ] `CohesiveZoneModel` — bilinear / exponential traction–separation
      (`t_0`, `δ_c`, `G_c`); mixed-mode Benzeggagh–Kenane
- [ ] `PhaseFieldFracture` — Griffith / AT1 / AT2 regularised
      functionals (`κ`, `G_c`, `l_0`); staggered solve reusing the
      `tpt-mat-phase-field` infrastructure
- [ ] `fracture_toughness_transition` — DBTT / master-curve
      (ASTM E1921) helper
- [ ] Verification test: phase-field fracture recovers the Griffith
      load for a 1-D bar; `K → G` Irwin consistency
- [ ] RFC: fracture mechanics (LEFM + CZM + phase-field)
- [ ] Example: `examples/phase-field-fracture-notch` (single-edge-notch
      tension)

### `tpt-mat-precipitation` — home group: thermo-kinetics (with `tpt-mat-calphad`)

- [ ] `ClassicalNucleation` — `I = I_0 exp(−ΔG*/kT) exp(−Q/kT)`; `ΔG*`
      from interfacial energy `γ` and driving force `Δg_v` (driving
      force from `tpt-mat-calphad`)
- [ ] `KwnModel` — Kampmann–Wagner–Numerical: discretised size classes,
      coupled nucleation + growth (`dR/dt`) + capillarity
- [ ] `LswCoarsening` — `R̄³ − R̄_0³ = K t`, `K` from `γ`, `D`, `c_eq`,
      `V_m`
- [ ] `PrecipitateState` — number density, mean radius, volume
      fraction, matrix supersaturation vs time
- [ ] `strengthening_increment` — Orowan bypass + shearing → `Δσ_y`,
      hand-off to `tpt-mat-hardening`
- [ ] Verification test: KWN conserves solute mass; late-stage slope
      → LSW `t^{1/3}`
- [ ] RFC: precipitation kinetics (CNT + KWN + LSW)
- [ ] Example: `examples/precipitation-age-hardening` (Al–Cu GP-zone
      → θ′ sweep)

### `tpt-mat-thermal` — home group: micro-mechanics (extends `tpt-mat-homogenization`)

- [ ] `effective_conductivity` — series / parallel / Hashin–Shtrikman /
      self-consistent / Maxwell–Garnett (2-phase and N-phase)
- [ ] `effective_cte` — Turner, Kerner, Rosen–Hashin bounds for
      composite thermal expansion
- [ ] `effective_specific_heat` — mass-weighted rule of mixtures
- [ ] `effective_diffusivity` — tortuosity / Bruggeman for porous &
      multiphase media (shared math with electrical conductivity)
- [ ] `interface_thermal_resistance` — Kapitza-resistance correction
- [ ] Generalise the 6×6 / scalar bound machinery in
      `tpt-mat-homogenization` from stiffness to any
      symmetric-positive transport tensor
- [ ] Verification test: conductivity HS bounds enclose the
      self-consistent estimate and collapse at zero contrast
- [ ] RFC: effective thermal & transport-property homogenization
- [ ] Example: `examples/composite-thermal-properties` (SiC/Al `k`,
      CTE vs `f`)

### `tpt-mat-dislocation` — home group: crystal-plasticity (extends `tpt-mat-hardening`)

- [ ] `DislocationDensityState` — per-slip-system `ρ_SSD`, `ρ_GND`,
      forest density
- [ ] `KocksMeckingEvolution` — `dρ/dγ = k_1 √ρ − k_2 ρ` (storage vs
      dynamic recovery); Taylor stress `τ = α μ b √ρ`
- [ ] `back_stress` — Armstrong–Frederick kinematic term from GND
      gradients
- [ ] `gnd_from_curvature` — Nye tensor → `ρ_GND` from a
      lattice-curvature field
- [ ] New `HardeningLaw::DislocationDensity` variant wired into
      `tpt-mat-hardening` / `tpt-mat-crystal-plasticity`
- [ ] Verification test: single-slip response reproduces Voce-like
      saturation; `ρ` stays non-negative
- [ ] RFC: dislocation-density-based hardening
- [ ] Example: `examples/dislocation-density-tension` (vs phenomenological
      Voce)

### `tpt-mat-hydrogen-embrittlement` — home group: degradation (bridges `tpt-mat-diffusion`)

- [ ] `HydrogenTransport` — Fick + trapping: Oriani local equilibrium
      and McNabb–Foster kinetic trapping (`N_T`, `E_B`, occupancy `θ_T`)
- [ ] `stress_driven_diffusion` —
      `∂C/∂t = ∇·(D∇C − D C V_H ∇σ_h / RT)` (hydrostatic-stress uphill
      flux) on a `tpt-science` grid
- [ ] `EmbrittlementCriterion` — HEDE critical-lattice-decohesion and
      HELP local-plasticity indicators; threshold `C_crit(σ_h)`
- [ ] `susceptibility_index` — from local H concentration + triaxiality
      field (feeds `tpt-mat-fracture` / `tpt-mat-fatigue-micro`)
- [ ] Verification test: trapping retards effective diffusivity by
      `D_eff = D_L / (1 + ∂C_T/∂C_L)`
- [ ] RFC: hydrogen transport with trapping + embrittlement criteria
- [ ] Example: `examples/hydrogen-embrittlement-notch` (H accumulation
      at a notch-tip stress field)

**Milestone:** A phase-field-fracture prediction whose local toughness
is modulated by a hydrogen-trapping field — two Phase-10 crates coupled
end-to-end.

---

## Ongoing / Cross-Phase

- [x] `cargo fmt` / `cargo clippy` / `cargo test` passing on every PR
      (this session: 121 tests across 16 crates, 0 failures; clippy
      reports only pedantic warnings under the workspace lint set,
      no errors).
- [ ] Maintain `cargo deny check licenses` passing (MIT chain enforcement, spec §8)
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

- [ ] High-temperature **oxidation** (Wagner parabolic scale growth,
      breakaway) — module in `tpt-mat-corrosion`; rename that group
      "environmental degradation"
- [ ] **Interfaces / grain boundaries** — GB energy, GB character
      distribution, Langmuir–McLean segregation, triple junctions
      (consolidate the scattered `GrainBoundaryMobility` /
      `GrainBoundaryDiffusion`)
- [ ] **Recrystallization** (static + dynamic) — JMAK-for-RX,
      Zener–Hollomon, nucleation criteria; distinct from
      `tpt-mat-grain-growth`'s curvature-driven model
- [ ] **Stereology / microstructure quantification** — ASTM E112 grain
      size, phase fraction from 2D/3D image data
- [ ] **Inverse / calibration** — fit constitutive parameters to
      experimental curves via `tpt-math-optimize-general` (spec §3)

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
| `tpt-mat-rve` (NEW) | `Rve` / `RveGrain` data model with per-grain orientation rotation of the 6×6 stiffness; `HomogenizationScheme` Voigt/Reuss/One-SelfConsistent; `SimpleHomogenizer`; **Bishop–Hill (1951) Taylor factor solver** with L2 pseudo-inverse (documented L1-vs-L2 caveat) |
| `tpt-mat-composite-micro` (NEW) | `rule_of_mixtures` (= Voigt); `dilute_estimate` (non-interacting inclusions); `mori_tanaka` two-phase closed form; `mori_tanaka_iterative` N-phase driver |
| `examples/homogenization-voigt-reuss` (NEW) | Sweeps Al+steel `f_steel = 0..1` with Voigt / Reuss / VRH / HS bounds side-by-side; VRH `E` ranges 70 → 200 GPa; HS bounds enclose Voigt |
| `examples/eshelby-inclusion` (NEW) | SiC-in-Al: Eshelby tensor (`S_h = 0.66`, `S_d = 0.47`), dilute strain-concentration tensor (ε_xx shielded from 1.0 to 0.31), Mori–Tanaka sweep `f = 0..1`; at `f = 0.1`: dilute `C_eff = 116 GPa`, MT `C_eff = 117 GPa` |
| `examples/taylor-factor-fcc` (NEW) | Bishop–Hill Taylor factor along `[001]→M=2.67`, `[011]→1.75`, `[111]→0.50`, `[012]→2.33`, `[112]→1.47`, `[123]→1.66`; 256-direction random average `M = 2.09` (L2 proxy; classical `M = 3.06` requires L1 Lemke solver) |
| `rfcs/0004-micromechanics-homogenization.md` (NEW) | Phase 5 RFC |

150 tests pass across 19 crates (Phase 1–5).  No build warnings
beyond pedantic clippy lints.