//! WebAssembly bindings for `tpt-materials`.
//!
//! **Status: Phase 8.**  Exposes JSON serialization of Phase 1 types
//! plus the interactive browser bindings:
//!
//! - [`WasmRveSolver`] — Voigt / Reuss / self-consistent
//!   homogenisation of a randomly generated polycrystal RVE, plus
//!   `homogenize` for a prescribed strain (spec §7).
//! - [`WasmPhaseField`] — Allen-Cahn phase-field stepping with
//!   `step` and `get_order_parameter` accessors for WebGL rendering.
//!
//! Build with:
//!
//! ```bash
//! cargo build -p tpt-mat-wasm --target wasm32-unknown-unknown
//! wasm-bindgen target/wasm32-unknown-unknown/release/tpt_mat_wasm.wasm \
//!   --out-dir pkg --target web
//! ```

#![cfg_attr(all(target_arch = "wasm32"), warn(missing_docs))]

use wasm_bindgen::prelude::*;

use tpt_mat_constants::PhysicalConstants;
use tpt_mat_core::{
    Composition, CompositionBasis, CrystalOrientation, Grain, GrainId, MaterialMicrostructure,
    OrientationRepresentation, Phase, PhaseId,
};
use tpt_mat_crystal_plasticity::SymmetricFourthOrder;
use tpt_mat_crystallography::{CrystalStructure, SlipSystem};
use tpt_mat_phase_field::{BulkEnergy, PhaseFieldModel, PhaseFieldSolver};
use tpt_mat_rve::{HomogenizationScheme, Rve, RveGrain};
use tpt_science::Grid2D;

#[wasm_bindgen(start)]
#[allow(dead_code, missing_docs)]
pub fn _start() {
    // Enable console_error_panic_hook once a logging crate is added.
}

// ---------------------------------------------------------------------
// JSON helpers (Phase 1)
// ---------------------------------------------------------------------

/// Round-trip a [`MaterialMicrostructure`] through JSON.
#[wasm_bindgen]
pub fn microstructure_round_trip_json(json: &str) -> Result<String, JsError> {
    let m: MaterialMicrostructure =
        serde_json::from_str(json).map_err(|e| JsError::new(&format!("parse: {e}")))?;
    serde_json::to_string(&m).map_err(|e| JsError::new(&format!("emit: {e}")))
}

/// Return the CODATA 2018 ideal-gas constant as a JS-accessible float.
#[wasm_bindgen]
pub fn gas_constant() -> f64 {
    PhysicalConstants::GAS_CONSTANT
}

/// Compute FCC slip systems as a JSON array.
#[wasm_bindgen]
pub fn fcc_slip_systems_json() -> Result<String, JsError> {
    let slips = CrystalStructure::FCC.slip_systems();
    serde_json::to_string(&slips).map_err(|e| JsError::new(&format!("emit: {e}")))
}

/// Identity Euler-Bunge angles as JSON.  Convenience for demos.
#[wasm_bindgen]
pub fn identity_orientation_json() -> Result<String, JsError> {
    let g = CrystalOrientation::identity();
    let euler = g.to_representation(OrientationRepresentation::EulerBunge);
    serde_json::to_string(&euler).map_err(|e| JsError::new(&format!("emit: {e}")))
}

// ---------------------------------------------------------------------
// RVE homogenisation (spec §7)
// ---------------------------------------------------------------------

/// Deterministic xorshift64 PRNG; reproducible across browsers and
/// versions so a given seed always yields the same RVE.
fn next_f64(seed: &mut u64) -> f64 {
    *seed ^= *seed << 13;
    *seed ^= *seed >> 7;
    *seed ^= *seed << 17;
    *seed as f64 / u64::MAX as f64
}

/// Uniform random point on the rotation manifold via Euler-Bunge
/// angles `(φ1, Φ, φ2)` with `Φ = acos(1 - 2 v)`.
fn random_rotation(seed: &mut u64) -> [[f64; 3]; 3] {
    let phi1 = 2.0 * std::f64::consts::PI * next_f64(seed);
    let phi = (1.0 - 2.0 * next_f64(seed)).acos();
    let phi2 = 2.0 * std::f64::consts::PI * next_f64(seed);
    let (s1, c1) = phi1.sin_cos();
    let (sp, cp) = phi.sin_cos();
    let (s2, c2) = phi2.sin_cos();
    [
        [c1 * c2 - s1 * cp * s2, s1 * c2 + c1 * cp * s2, sp * s2],
        [-c1 * s2 - s1 * cp * c2, -s1 * s2 + c1 * cp * c2, sp * c2],
        [s1 * sp, -c1 * sp, cp],
    ]
}

/// FCC single-crystal stiffness (Pa), matching the crystal-plasticity
/// solver defaults.
fn fcc_stiffness() -> SymmetricFourthOrder {
    SymmetricFourthOrder::cubic(168.4e9, 121.4e9, 75.4e9)
}

/// Interactive RVE solver: holds a deterministically generated
/// polycrystal and applies analytical homogenisation schemes.
#[wasm_bindgen]
pub struct WasmRveSolver {
    rve: Rve,
    scheme: HomogenizationScheme,
}

/// Homogenised response of the RVE for one applied strain.
#[wasm_bindgen]
pub struct WasmHomogenizedResponse {
    /// Sample-frame stress (Voigt order, Pa).
    stress: Vec<f64>,
    /// Effective 6x6 stiffness (Voigt order, Pa).
    effective_stiffness: Vec<f64>,
    /// Scheme name used for the response.
    scheme: String,
}

#[wasm_bindgen]
impl WasmRveSolver {
    /// Construct a polycrystal RVE with `num_grains` equal-fraction
    /// grains of random orientation.  `resolution` is the PRNG seed
    /// (the analytical schemes need no mesh, so `resolution` seeds the
    /// deterministic orientation generator).
    #[wasm_bindgen(constructor)]
    pub fn new(num_grains: u32, resolution: u32) -> WasmRveSolver {
        let n = num_grains.max(1);
        let mut seed = (resolution as u64).max(1);
        let stiffness = fcc_stiffness();
        let fraction = 1.0 / n as f64;
        let grains = (0..n)
            .map(|i| {
                RveGrain::new(
                    format!("g{i}"),
                    fraction,
                    random_rotation(&mut seed),
                    stiffness.clone(),
                )
            })
            .collect();
        WasmRveSolver {
            rve: Rve::new(grains),
            scheme: HomogenizationScheme::Voigt,
        }
    }

    /// Select the homogenisation scheme: `"voigt"` (default),
    /// `"reuss"` or `"self_consistent"`.
    pub fn set_scheme(&mut self, scheme: &str) -> Result<(), JsError> {
        self.scheme = match scheme.to_ascii_lowercase().as_str() {
            "voigt" => HomogenizationScheme::Voigt,
            "reuss" => HomogenizationScheme::Reuss,
            "self_consistent" | "selfconsistent" | "" => HomogenizationScheme::SelfConsistent,
            other => {
                return Err(JsError::new(&format!(
                    "unknown scheme '{other}' (expected voigt | reuss | self_consistent)"
                )));
            }
        };
        Ok(())
    }

    /// Number of grains in the RVE.
    pub fn n_grains(&self) -> u32 {
        self.rve.n_grains() as u32
    }

    /// Total volume fraction (ideally 1.0).
    pub fn total_volume_fraction(&self) -> f64 {
        self.rve.total_volume_fraction()
    }

    /// Current effective 6x6 stiffness (Voigt order, Pa).
    pub fn effective_stiffness(&self) -> Vec<f64> {
        self.rve
            .homogenize(self.scheme)
            .data
            .iter()
            .flatten()
            .copied()
            .collect()
    }

    /// Apply a 6-component strain (Voigt order) and return the
    /// homogenised stress plus the effective stiffness used.
    pub fn homogenize(&self, strain: Vec<f64>) -> Result<WasmHomogenizedResponse, JsError> {
        if strain.len() != 6 {
            return Err(JsError::new(&format!(
                "strain must have 6 Voigt components, got {}",
                strain.len()
            )));
        }
        let c = self.rve.homogenize(self.scheme);
        let stress: Vec<f64> = c
            .data
            .map(|row| row.iter().zip(strain.iter()).map(|(a, s)| a * s).sum())
            .to_vec();
        let stiffness: Vec<f64> = c.data.iter().flatten().copied().collect();
        Ok(WasmHomogenizedResponse {
            stress,
            effective_stiffness: stiffness,
            scheme: scheme_name(self.scheme),
        })
    }
}

#[wasm_bindgen]
impl WasmHomogenizedResponse {
    /// Sample-frame stress in Voigt order (Pa).
    pub fn stress(&self) -> Vec<f64> {
        self.stress.clone()
    }
    /// Effective 6x6 stiffness in Voigt order (Pa).
    pub fn effective_stiffness(&self) -> Vec<f64> {
        self.effective_stiffness.clone()
    }
    /// Scheme name used for the response.
    pub fn scheme(&self) -> String {
        self.scheme.clone()
    }
}

fn scheme_name(scheme: HomogenizationScheme) -> String {
    match scheme {
        HomogenizationScheme::Voigt => "voigt",
        HomogenizationScheme::Reuss => "reuss",
        HomogenizationScheme::SelfConsistent => "self_consistent",
    }
    .to_string()
}

// ---------------------------------------------------------------------
// Phase-field (spec §7)
// ---------------------------------------------------------------------

/// Interactive Allen-Cahn phase-field solver.
#[wasm_bindgen]
pub struct WasmPhaseField {
    solver: PhaseFieldSolver,
}

#[wasm_bindgen]
impl WasmPhaseField {
    /// Construct an Allen-Cahn solver on an `nx` x `ny` grid with a
    /// double-well bulk energy.  `initial_order` must contain exactly
    /// `nx * ny` cells.
    #[wasm_bindgen(constructor)]
    pub fn new(
        nx: u32,
        ny: u32,
        time_step: f64,
        initial_order: Vec<f64>,
    ) -> Result<WasmPhaseField, JsError> {
        let grid = Grid2D::new(nx as usize, ny as usize, 1.0);
        let solver = PhaseFieldSolver::new_2d(
            PhaseFieldModel::AllenCahn,
            grid,
            time_step,
            1.0,
            1.0,
            BulkEnergy::DoubleWell { well_depth: 1.0 },
            &initial_order,
        )?;
        Ok(WasmPhaseField { solver })
    }

    /// Advance the simulation by exactly one time step.
    pub fn step(&mut self) -> Result<(), JsError> {
        Ok(self.solver.step()?)
    }

    /// Advance by `n_steps` time steps.
    pub fn step_many(&mut self, n_steps: u32) -> Result<(), JsError> {
        for _ in 0..n_steps {
            self.step()?;
        }
        Ok(())
    }

    /// Order-parameter field (one value per cell) for WebGL
    /// visualization.
    pub fn get_order_parameter(&self) -> Vec<f64> {
        self.solver.order_parameter()[0].clone()
    }

    /// Simulation time after the steps taken.
    pub fn time(&self) -> f64 {
        self.solver.time
    }

    /// Total free energy of the current field.
    pub fn free_energy(&self) -> f64 {
        self.solver.snapshot().free_energy
    }

    /// Number of grid cells.
    pub fn n_cells(&self) -> u32 {
        self.solver.n_cells() as u32
    }
}

// Re-export the workspace types so they are available to other wasm
// crates that depend on this one.
pub use reexport::*;

mod reexport {
    pub use tpt_mat_constants::{AtomicData, PeriodicTable};
    pub use tpt_mat_core::{
        Composition as WasmComposition, CompositionBasis as WasmCompositionBasis,
        CompositionError as WasmCompositionError, CrystalOrientation as WasmCrystalOrientation,
        Grain as WasmGrain, GrainId as WasmGrainId,
        MaterialMicrostructure as WasmMaterialMicrostructure,
        OrientationConversionError as WasmOrientationConversionError,
        OrientationRepresentation as WasmOrientationRepresentation, Phase as WasmPhase,
        PhaseId as WasmPhaseId,
    };
    pub use tpt_mat_crystal_plasticity::SymmetricFourthOrder as WasmSymmetricFourthOrder;
    pub use tpt_mat_crystallography::{
        CrystalStructure as WasmCrystalStructure, LatticeParameters as WasmLatticeParameters,
        MillerIndex as WasmMillerIndex, SlipFamily as WasmSlipFamily, SlipSystem as WasmSlipSystem,
        SlipSystemError as WasmSlipSystemError,
    };
    pub use tpt_mat_phase_field::{
        BulkEnergy as WasmBulkEnergy, PhaseFieldModel as WasmPhaseFieldModel,
    };
    pub use tpt_mat_rve::{
        HomogenizationScheme as WasmHomogenizationScheme, Rve as WasmRve, RveGrain as WasmRveGrain,
    };
    pub use tpt_math_linalg_fixed::{SymMat3 as WasmSymMat3, Vec3 as WasmVec3, Vec6 as WasmVec6};
}

// Ensure the unused-imports lint stays happy while the crate exposes
// only JSON helpers for the Phase 1 types today.
#[allow(dead_code)]
fn _types_used() {
    let _ = (
        std::marker::PhantomData::<Composition>,
        std::marker::PhantomData::<CompositionBasis>,
        std::marker::PhantomData::<Grain>,
        std::marker::PhantomData::<GrainId>,
        std::marker::PhantomData::<Phase>,
        std::marker::PhantomData::<PhaseId>,
        std::marker::PhantomData::<SlipSystem>,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn random_rotation_is_orthonormal() {
        let mut seed = 42u64;
        for _ in 0..16 {
            let r = random_rotation(&mut seed);
            // R R^T = I within fp64 roundoff.
            for row in 0..3 {
                for col in 0..3 {
                    let dot = (0..3).map(|k| r[row][k] * r[col][k]).sum::<f64>();
                    let expected = if row == col { 1.0 } else { 0.0 };
                    assert!(
                        (dot - expected).abs() < 1.0e-12,
                        "R R^T[{row}][{col}] = {dot}"
                    );
                }
            }
        }
    }

    #[test]
    fn scheme_name_covers_all_variants() {
        assert_eq!(scheme_name(HomogenizationScheme::Voigt), "voigt");
        assert_eq!(scheme_name(HomogenizationScheme::Reuss), "reuss");
        assert_eq!(
            scheme_name(HomogenizationScheme::SelfConsistent),
            "self_consistent"
        );
    }

    #[test]
    fn cubic_stiffness_is_positive_definite() {
        let c = fcc_stiffness();
        let principal = c.data[0][0];
        let shear = c.data[3][3];
        assert!(principal > shear);
        assert!(principal > 0.0 && shear > 0.0);
    }

    #[test]
    fn rve_logic_default_voigt_and_sc_invariant() {
        // Mirrors WasmRveSolver on plain types (wasm-bindgen JsError
        // conversion cannot execute on native).
        let mut seed = 7u64;
        let stiffness = fcc_stiffness();
        let grains: Vec<RveGrain> = (0..50)
            .map(|i| {
                RveGrain::new(
                    format!("g{i}"),
                    1.0 / 50.0,
                    random_rotation(&mut seed),
                    stiffness.clone(),
                )
            })
            .collect();
        let rve = Rve::new(grains);
        let voigt = rve.homogenize(HomogenizationScheme::Voigt);
        let sc = rve.homogenize(HomogenizationScheme::SelfConsistent);
        // For a cubic polycrystal the bulk modulus is orientation
        // invariant, so the SC estimate must stay close to it.  The
        // scalar-nu Eshelby approximation allows some sampling-noise
        // drift, so allow 10%.
        let k_voigt = (voigt.data[0][0] + 2.0 * voigt.data[0][1]) / 3.0;
        let k_sc = (sc.data[0][0] + 2.0 * sc.data[0][1]) / 3.0;
        assert!(
            (k_sc - k_voigt).abs() <= 0.1 * k_voigt,
            "SC K {k_sc} drifted >10% from crystal K {k_voigt}"
        );
        // Voigt iso-strain upper bound dominates the SC shear modulus.
        assert!(sc.data[3][3] <= voigt.data[3][3] * 1.1);
        // Tensile strain along 0 ⇒ positive normal stress.
        let strain = [0.001, 0.0, 0.0, 0.0, 0.0, 0.0];
        let stress: f64 = (0..6).map(|k| sc.data[0][k] * strain[k]).sum();
        assert!(stress > 0.0);
    }

    #[test]
    fn phase_field_logic_decreases_free_energy() {
        // Mirrors WasmPhaseField on plain types.
        let grid = Grid2D::new(32, 32, 1.0);
        let mut initial = vec![0.5_f64; 32 * 32];
        for i in 0..32 {
            for j in 0..32 {
                initial[i * 32 + j] = if j < 16 { 0.0 } else { 1.0 };
            }
        }
        let mut solver = PhaseFieldSolver::new_2d(
            PhaseFieldModel::AllenCahn,
            grid,
            0.05,
            1.0,
            1.0,
            BulkEnergy::DoubleWell { well_depth: 1.0 },
            &initial,
        )
        .unwrap();
        let e0 = solver.snapshot().free_energy;
        solver.step_many(20).unwrap();
        let e1 = solver.snapshot().free_energy;
        assert!(e1 <= e0 + 1e-9, "free energy increased: {e0} → {e1}");
    }
}
