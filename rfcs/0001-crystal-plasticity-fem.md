# RFC 0001 — Crystal-Plasticity FEM

| | |
|---|---|
| **Status** | Accepted |
| **Author** | TPT Solutions |
| **Created** | 2026-09-02 |
| **Phase** | 2 |

## Summary

Adopt a standard rate-dependent single-crystal plasticity model
(Pécha–Hutchinson–Bronkhorst / Asaro–Needleman) as the constitutive
backbone of `tpt-mat-crystal-plasticity`, and wire it into the
`tpt-fem` Newton-Raphson loop as `CpFemSolver`.

## Motivation

`todo.md` Phase 2 specifies a single-crystal tension example with
**correct slip activation**. Today the data model (`MaterialMicrostructure`,
`CrystalOrientation`, slip systems) is in place from Phase 1, but the
constitutive update that takes a strain increment into slip rates,
hardening, and lattice rotation does not exist. Phase 5 homogenization
(RVE) and Phase 6 fatigue both depend on this layer.

## Constitutive model

Multiplicative decomposition of the deformation gradient

```
F = F^e · F^p
```

with the plastic flow rule

```
L^p = Ḟ^p · (F^p)^(-1) = Σ_α  γ̇_α  (s_α ⊗ n_α)
```

where `α` indexes slip systems, `γ̇_α` is the slip rate, `s_α` and `n_α`
are the slip direction and plane normal in the crystal frame.

Rate-dependent (viscoplastic) flow rule:

```
γ̇_α = γ̇_0 · |τ_α / g_α|^(1/m) · sign(τ_α)
```

with `γ̇_0` the reference strain rate, `m` the rate sensitivity,
`g_α` the current slip resistance (CRSS), and `τ_α` the resolved
shear stress on slip system `α`.

Hardening (Voce, default):

```
ġ_α = Σ_β  h_αβ · |γ̇_β|
h_αβ = h_0 · q_αβ · (g_sat − g_β)            (diagonal part)
```

with the latent-hardening matrix `q_αβ` (`q = 1.0` self, `q = 1.4`
latent, default).

## Lattice rotation update

The elastic spin `W^e` is computed from the polar decomposition of
`F^e`, and the crystal orientation `R` is updated as

```
R_{n+1} = ΔR · R_n
```

where `ΔR = exp(Δt · (W^e − W^p))` and `W^p` is the plastic spin.

## Solver API

```rust
let model = CrystalPlasticityModel::fcc(
    ElasticTensor::cubic_anisotropic(...),
    HardeningLaw::Voce(VoceParameters { h_0: 500.0, g_sat: 250.0, g_0: 80.0 }),
    RateSensitivity { gamma_0: 1e-3, m: 0.02 },
);

let solver = CpFemSolver::new(mesh, model, boundary_conditions);
let result: CpFemResult = solver.solve_increment(strain_increment)?;
```

`CpFemSolver::solve_increment` uses
[`tpt_fem_solve::newton`](https://docs.rs/tpt-fem-solve) with a
physics-crate-supplied residual and Jacobian that integrates the
constitutive update at every Gauss point between calls.

### Per-GP integration scheme

Implicit backward Euler at each GP, one Newton iteration per time
step (sufficient for `m ≤ 0.05`).

## Phase 2 scope

This RFC covers:

- The `CrystalPlasticityModel` and `HardeningLaw` types in
  `tpt-mat-crystal-plasticity` and `tpt-mat-hardening`.
- The single-GP integration update (rate-dependent, Voce).
- `CpFemSolver` for **single-element** meshes: one hex8, one GP, a
  uniaxial tension BC. This is enough to validate the constitutive
  update end-to-end against the closed-form Schmid prediction.

## Out of scope (later phases)

- Multi-GP / multi-element CP-FEM with `tpt-fem-assembly`. **Phase 5.**
- Latent hardening (only the diagonal `q_αα = 1.0` part of Voce for now).
- Finite-strain `F` update — Phase 2 uses small-strain Cauchy
  elasticity at the GP. Switching to multiplicative `F = F^e F^p` is
  mechanical; see `RFC 0001-followup`.
- Texture evolution (handled by `tpt-mat-texture`).
- Taylor-factor calculation (handled by `tpt-mat-texture`).

## Verification

- Phase 2: Taylor factor for random FCC texture ≈ 3.06
  (Bhattacharyya et al., *Acta Metall.* 1991).
- Phase 2: single-crystal tension along `[001]` activates
  `(111)[1̄01]` first at `τ = g_0`, stress step Δσ = g_0 / 0.4082.
- Phase 5: Hill–Mandel condition satisfied across an RVE.

## References

- Peirce, Asaro, Needleman, *Acta Metall.* 30 (1982) 1087.
- Asaro, Needleman, *Acta Metall.* 33 (1985) 923.
- Bronkhorst, Hanson, Morris, *Int. J. Plast.* 8 (1992) 411.
- Roters et al., *Comput. Mater. Sci.* 53 (2012) 327
  (DAMASK CP-FEM review).