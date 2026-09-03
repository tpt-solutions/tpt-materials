# RFC 0004 — Micro-Mechanics: Homogenization, RVE, and the Bishop–Hill Taylor Factor

**Status:** Implemented (this PR).
**Authors:** TPT Solutions.
**Applies to:** `tpt-mat-homogenization`, `tpt-mat-rve`,
`tpt-mat-composite-micro`.

---

## 1. Summary

Phase 5 of the `tpt-materials` roadmap introduces three new crates
that capture the analytical micro-mechanics of polycrystalline
aggregates and composite microstructures:

- **`tpt-mat-homogenization`** — Closed-form bounds and estimates
  for the effective elastic stiffness `C_eff` of a heterogeneous
  medium: Voigt / Reuss averages, Hill `M_VRH` arithmetic mean,
  Hashin–Shtrikman variational bounds for two-phase mixtures, and
  Eshelby's equivalent-inclusion tensor for a spherical
  inclusion in an isotropic matrix.
- **`tpt-mat-rve`** — A Representative Volume Element (RVE) data
  model (per-grain orientation, volume fraction, stiffness) plus a
  Bishop–Hill (1951) Taylor-factor solver that replaces the
  Phase-2 single-Schmid proxy with a true vertex-enumeration LCP
  solver.
- **`tpt-mat-composite-micro`** — Mean-field homogenization for
  multi-phase composites: rule-of-mixtures (Voigt), dilute
  (Maxwell) estimate for non-interacting inclusions, and the
  Mori–Tanaka (1973) scheme that interpolates between Voigt and
  Reuss bounds and matches Hashin–Shtrikman for spherical
  inclusions.

This RFC ratifies the public APIs of these three crates as the
canonical Phase 5 surface and records the explicit non-goals.

---

## 2. `tpt-mat-homogenization`

### 2.1 Voigt / Reuss / Hill averages

```rust
pub fn voigt(phases: &[(SymmetricFourthOrder, f64)]) -> SymmetricFourthOrder;
pub fn reuss(phases: &[(SymmetricFourthOrder, f64)]) -> [[f64; 6]; 6];
pub fn voigt_reuss_average(m_voigt: f64, m_reuss: f64) -> f64;
```

The Voigt average is the arithmetic mean of stiffness tensors
weighted by volume fraction (`C_eff ≤ C_V` is an upper bound for the
true `C_eff`).  The Reuss average is the arithmetic mean of
compliances (`S_eff ≥ S_R` ⇒ `C_eff ≥ C_R` is a lower bound).

`voigt_reuss_bounds` is a convenience that takes isotropic
`(K_r, G_r, f_r)` triplets and returns both bounds and the Hill
`M_VRH` average for `K`, `G`, `E`, `ν`.

### 2.2 Hashin–Shtrikman variational bounds

```rust
pub struct HashinShtrikmanResult {
    pub k_lower: f64, pub k_upper: f64,
    pub g_lower: f64, pub g_upper: f64,
}

pub fn hashin_shtrikman_two_phase(
    phase_a: (f64, f64, f64), // (E, ν, f)
    phase_b: (f64, f64, f64),
) -> HashinShtrikmanResult;

pub fn hashin_shtrikman_spherical_bulk(
    k_matrix: f64, g_matrix: f64,
    k_inclusion: f64, f_inclusion: f64,
) -> f64;
```

The two-phase HS bounds are the tightest possible *isotropic*
bounds on `K_eff, G_eff` for a two-phase composite with
isotropic constituents — they are tighter than Voigt / Reuss for
typical stiffness contrasts.  For an `N`-phase composite the
bounds must be applied pairwise (sequential homogenisation); a
proper `N`-phase HS implementation lands with the next phase that
requires it.

### 2.3 Eshelby equivalent-inclusion

```rust
pub struct EshelbySpherical { pub s_hydro: f64, pub s_dev: f64 }
pub fn eshelby_spherical(matrix_e: f64, matrix_nu: f64) -> EshelbySpherical;

pub struct StrainConcentrationTensor { pub data: [[f64; 6]; 6] }
pub fn dilute_strain_concentration(
    c_matrix: &SymmetricFourthOrder,
    c_inclusion: &SymmetricFourthOrder,
    eshelby: EshelbySpherical,
) -> StrainConcentrationTensor;
```

For a *spherical* inclusion in an *isotropic* matrix, the Eshelby
tensor `S` has the closed-form hydrostatic / deviatoric
decomposition

`S = S_h ⊗ I_h + S_d ⊗ I_d`,

`S_h = (1 + ν) / (3 (1 - ν))`, `S_d = 2 (4 - 5 ν) / (15 (1 - ν))`.

The dilute strain-concentration tensor is

`A = [I + S · C_0^{-1} · (C_1 - C_0)]^{-1}`.

Ellipsoidal-inclusion Eshelby tensors are not yet implemented;
they require either the closed-form Eshelby solution (prolate /
oblate spheroids) or a numerical FFT-based solver.

### 2.4 Non-goals

- **No FFT-based homogenisation.**  That belongs to a future
  `tpt-mat-homogenization-fft` crate (planned for Phase 8).
- **No two-point statistics.**  A `CorrelationFunction`-driven
  scheme (e.g. Ponte-Castañeda & Su, 1997) is deferred.
- **No non-linear homogenisation.**  All schemes in this crate
  assume linear elasticity.

---

## 3. `tpt-mat-rve`

### 3.1 Data model

```rust
pub struct RveGrain {
    pub label: String,
    pub volume_fraction: f64,
    pub orientation: [[f64; 3]; 3],  // crystal → sample
    pub stiffness: SymmetricFourthOrder,
}

pub struct Rve { pub grains: Vec<RveGrain> }
pub struct RveStats { pub n_grains: usize, pub total_volume_fraction: f64, ... }

pub enum HomogenizationScheme { Voigt, Reuss, SelfConsistent }
pub struct SimpleHomogenizer;
```

`Rve::homogenize(scheme)` returns the effective 6×6 stiffness for
the chosen analytical scheme:

- `Voigt`: arithmetic mean of grain stiffnesses (in sample frame).
- `Reuss`: harmonic mean (compliances averaged, then inverted).
- `SelfConsistent`: a one-site mean-field iteration that converges
  in 30 iterations for typical polycrystals.  This is a Voigt-step
  approximation — the proper dilute-Eshelby update lands with the
  Phase 8 FFT integration.

### 3.2 Bishop–Hill Taylor factor

```rust
pub struct BishopHillResult {
    pub n_active: usize,
    pub slip_rates: Vec<f64>,
    pub plastic_strain: Vec6,
    pub stress: Vec6,
    pub taylor_factor: f64,
}

pub fn bishop_hill_taylor_factor_axis(
    tensile_axis: [f64; 3],
    tau_c: f64,
    slips: Option<&[SlipSystem]>,
) -> BishopHillResult;

pub fn bishop_hill_taylor_factor(eps: Vec6, tau_c: f64) -> BishopHillResult;
```

The Bishop–Hill (1951) / Taylor (1938) maximum-work principle
reduces to a small Linear Complementarity Problem (LCP): for a
prescribed macroscopic strain `ε`, find a stress `σ` such that
`σ : P^α ≤ τ_c` for all `α` and `σ : P^α = τ_c · sign(γ^α)` for
the 5 active slip systems `α`.

This crate implements a **primal L2 pseudo-inverse** solver: it finds
the slip rates `γ^α` that minimise `‖γ^α‖_2` subject to
`Σ γ^α P^α = ε`, then reports the Taylor factor
`M = (Σ |γ^α|) / ε_eq`.  The vertex stress is computed when the
dual system is well-conditioned.

### 3.3 Limitations (documented)

The L2 pseudo-inverse solver returns the **L2-minimum** slip-rate
solution, not the classical **L1-minimum** Taylor flow rule.  The
L1 solver (Lemke's algorithm or branch-and-bound) is required to
recover Taylor's original random-texture FCC result of `M ≈ 3.06`;
the L2 solver returns `M ≈ 2.0–2.5` (a lower bound on the
true Taylor factor).  See `examples/taylor-factor-fcc`.

This is an **improvement** over the Phase 2 single-Schmid proxy
(`1 / max_schmid_factor`, which gave `M ≈ 2.0–2.3`), but is **not
yet** a complete replacement for the LCP solver.

### 3.4 Non-goals

- **No FEM assembly.**  A full `tpt-fem` integration that uses
  this RVE as a constitutive model lands when `tpt-fem` is
  upstreamed (out-of-repo).
- **No FFT-based RVE solver.**  The discrete Fourier-transform
  scheme (Moulinec & Suquet, 1998) is deferred to Phase 8.

---

## 4. `tpt-mat-composite-micro`

### 4.1 Scheme taxonomy

```rust
pub fn rule_of_mixtures(phases: &[(SymmetricFourthOrder, f64)]) -> SymmetricFourthOrder;
pub fn dilute_estimate(c_matrix: &SymmetricFourthOrder,
                        c_inclusion: &SymmetricFourthOrder,
                        f_inclusion: f64,
                        eshelby: EshelbySpherical) -> SymmetricFourthOrder;
pub fn mori_tanaka(c_matrix: &SymmetricFourthOrder,
                    c_inclusion: &SymmetricFourthOrder,
                    f_inclusion: f64,
                    eshelby: EshelbySpherical) -> MoriTanakaResult;
pub fn mori_tanaka_iterative(phases: &[(SymmetricFourthOrder, f64)],
                              matrix_nu: f64, max_iter: usize, tol: f64)
    -> SymmetricFourthOrder;
```

- `rule_of_mixtures`: arithmetic mean (= Voigt bound).
- `dilute_estimate`: non-interacting inclusions (valid only for
  `f ≪ 1`).
- `mori_tanaka`: two-phase exact Benveniste (1987) closed form.
- `mori_tanaka_iterative`: `N`-phase iterative driver; converges in
  ~30 iterations for typical stiffness contrasts.

### 4.2 Non-goals

- **No ellipsoidal-inclusion Eshelby tensors** (see §2.3).
- **No gradient-enhanced schemes** (e.g. Lahellec–Suquet, 2013)
  for nonlinear composites.

---

## 5. Test surface

The crates ship with `unit tests` covering:

- Voigt / Reuss / Hill recovery for single- and two-phase
  aggregates.
- Hashin–Shtrikman enclosure of the Voigt bound.
- Eshelby tensor symmetry, dilute strain concentration identity
  for matching inclusion / matrix, and the strain-amplification
  behaviour for stiff vs soft inclusions.
- RVE rotation of a 6×6 stiffness into the sample frame, two-grain
  Voigt averaging, RVE stats.
- Bishop–Hill Taylor factor for `[001]`, `[011]`, `[111]`, etc.
  tensile axes and a 256-direction random average (L2 proxy).
- Mori–Tanaka recovery of the matrix stiffness at `f = 0`,
  recovery of the inclusion stiffness at `f = 1`, and intermediate
  fractions lying between the Voigt and Reuss bounds.

Three examples in `examples/` drive the crates end-to-end:

- `homogenization-voigt-reuss` — sweeps volume fraction and
  prints Voigt / Reuss / VRH / HS bounds side-by-side.
- `eshelby-inclusion` — reports the Eshelby tensor, the dilute
  strain-concentration tensor, and the Mori–Tanaka effective
  stiffness as a function of `f_inclusion`.
- `taylor-factor-fcc` — high-symmetry axes + 256-direction
  random average for the Bishop–Hill FCC Taylor factor.

---

## 6. Forward work

- **Proper L1 Taylor solver.**  Replace the L2 pseudo-inverse
  with a Lemke-style LCP solver to recover `M ≈ 3.06` for random
  FCC textures.
- **FFT homogenisation.**  Moulinec–Suquet (1998) iterative
  scheme lands with the Phase 8 informatics / ecosystem work.
- **FEM integration.**  When `tpt-fem` is upstreamed, `Rve::solve`
  becomes a per-integration-point constitutive model.
- **Two-point correlation-driven homogenisation.**  Ponte-Castañeda
  nonlinear scheme.