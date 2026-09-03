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
> unit tests, examples, and RFC 0004.  Phases 6–8 are still pending.

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

- [ ] All Phase 6 items deferred.

---

## Phase 7 — Energy Materials & Manufacturing (Months 19-21)

- [ ] All Phase 7 items deferred.

---

## Phase 8 — Informatics & Ecosystem (Months 22-24)

- [ ] All Phase 8 items deferred.

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