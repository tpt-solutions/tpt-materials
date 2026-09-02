//! WebAssembly bindings for `tpt-materials`.
//!
//! **Status: Phase 1 scaffold.**  Only exposes JSON serialization of
//! Phase 1 types today.  The full `WasmRveSolver` / `WasmPhaseField`
//! bindings arrive in Phase 8.
//!
//! Build with:
//!
//! ```bash
//! cargo build -p tpt-mat-wasm --target wasm32-unknown-unknown
//! wasm-bindgen ... --out-dir pkg
//! ```

#![cfg_attr(all(target_arch = "wasm32"), warn(missing_docs))]

use wasm_bindgen::prelude::*;

use tpt_mat_constants::PhysicalConstants;
use tpt_mat_core::{
    Composition, CompositionBasis, CrystalOrientation, Grain, GrainId, MaterialMicrostructure,
    OrientationRepresentation, Phase, PhaseId,
};
use tpt_mat_crystallography::{CrystalStructure, SlipSystem};

#[wasm_bindgen(start)]
#[allow(dead_code, missing_docs)]
pub fn _start() {
    // Enable console_error_panic_hook once a logging crate is added in Phase 8.
}

/// Round-trip a [`MaterialMicrostructure`] through JSON.  Useful in
/// integration tests until the typed bindings land in Phase 8.
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

/// Compute FCC slip systems as a JSON array.  Demonstrates that the
/// Phase 1 crystallography crate is callable from WASM.
#[wasm_bindgen]
pub fn fcc_slip_systems_json() -> Result<String, JsError> {
    let slips = CrystalStructure::FCC.slip_systems();
    serde_json::to_string(&slips).map_err(|e| JsError::new(&format!("emit: {e}")))
}

/// Identity Euler-Bunge angles as JSON.  Convenience for Phase 8 demos.
#[wasm_bindgen]
pub fn identity_orientation_json() -> Result<String, JsError> {
    let g = CrystalOrientation::identity();
    let euler = g.to_representation(OrientationRepresentation::EulerBunge);
    serde_json::to_string(&euler).map_err(|e| JsError::new(&format!("emit: {e}")))
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
    pub use tpt_mat_crystallography::{
        CrystalStructure as WasmCrystalStructure, LatticeParameters as WasmLatticeParameters,
        MillerIndex as WasmMillerIndex, SlipFamily as WasmSlipFamily, SlipSystem as WasmSlipSystem,
        SlipSystemError as WasmSlipSystemError,
    };
    pub use tpt_math_linalg_fixed::{SymMat3 as WasmSymMat3, Vec3 as WasmVec3, Vec6 as WasmVec6};
}

// Ensure the unused-imports lint stays happy while the crate is a scaffold.
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
