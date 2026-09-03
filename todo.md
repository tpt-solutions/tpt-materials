# tpt-materials — Project Todo

Organization: **TPT Solutions** · License: **MIT OR Apache-2.0** (dual)

> **Status snapshot — 2026-09-02.**  Phase 1 (foundation) is complete
> from the original scaffold; Phase 2 (crystal plasticity) and Phase 3
> (phase-field) have been **implemented in this session** with library
> code, unit tests, examples, and RFC stubs.  Phases 4–8 are still
> pending.  See the per-phase sections for an honest breakdown of what
> is implemented, partially implemented, or stubbed.

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

- [ ] All Phase 4 items deferred.

---

## Phase 5 — Micro-Mechanics (Months 13-15)

**Crates:** `tpt-mat-rve`, `tpt-mat-homogenization`, `tpt-mat-composite-micro`
**Substrate:** `tpt-fem` (cross-repo; published on crates.io)

- [ ] All Phase 5 items deferred.
- [ ] **TODO**: a proper Taylor-factor Bishop–Hill LCP solver
      (lands here with the homogenization stack).

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
      (this session: 86 tests across 13 crates, 0 failures; clippy
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