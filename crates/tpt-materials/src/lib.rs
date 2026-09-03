//! `tpt-materials` — umbrella / facade crate.
//!
//! Every spec §5 domain crate is re-exported behind an optional
//! feature so callers can pull in only what they need:
//!
//! ```toml
//! [dependencies]
//! tpt-materials = { version = "0.1", features = ["full"] }
//! ```
//!
//! or selectively:
//!
//! ```toml
//! [dependencies]
//! tpt-materials = { version = "0.1", features = ["crystal-plasticity", "homogenization"] }
//! ```

#![warn(missing_docs)]

#[cfg(feature = "crystallography")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_crystallography;
#[cfg(feature = "crystal-plasticity")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_crystal_plasticity;
#[cfg(feature = "hardening")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_hardening;
#[cfg(feature = "texture")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_texture;
#[cfg(feature = "phase-field")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_phase_field;
#[cfg(feature = "grain-growth")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_grain_growth;
#[cfg(feature = "solidification")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_solidification;
#[cfg(feature = "diffusion")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_diffusion;
#[cfg(feature = "phase-transform")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_phase_transform;
#[cfg(feature = "calphad")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_calphad;
#[cfg(feature = "homogenization")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_homogenization;
#[cfg(feature = "rve")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_rve;
#[cfg(feature = "composite-micro")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_composite_micro;
#[cfg(feature = "damage")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_damage;
#[cfg(feature = "fatigue")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_fatigue;
#[cfg(feature = "creep")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_creep;
#[cfg(feature = "heat-treatment")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_heat_treatment;
#[cfg(feature = "database")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_database;
#[cfg(feature = "wasm")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_wasm;

#[cfg(feature = "crystallography")]
pub mod crystallography {
    //! Re-export of `tpt-mat-crystallography`.
    pub use tpt_mat_crystallography::*;
}

#[cfg(feature = "crystal-plasticity")]
pub mod crystal_plasticity {
    //! Re-export of `tpt-mat-crystal-plasticity`.
    pub use tpt_mat_crystal_plasticity::*;
}

#[cfg(feature = "hardening")]
pub mod hardening {
    //! Re-export of `tpt-mat-hardening`.
    pub use tpt_mat_hardening::*;
}

#[cfg(feature = "texture")]
pub mod texture {
    //! Re-export of `tpt-mat-texture`.
    pub use tpt_mat_texture::*;
}

#[cfg(feature = "phase-field")]
pub mod phase_field {
    //! Re-export of `tpt-mat-phase-field`.
    pub use tpt_mat_phase_field::*;
}

#[cfg(feature = "grain-growth")]
pub mod grain_growth {
    //! Re-export of `tpt-mat-grain-growth`.
    pub use tpt_mat_grain_growth::*;
}

#[cfg(feature = "solidification")]
pub mod solidification {
    //! Re-export of `tpt-mat-solidification`.
    pub use tpt_mat_solidification::*;
}

#[cfg(feature = "diffusion")]
pub mod diffusion {
    //! Re-export of `tpt-mat-diffusion`.
    pub use tpt_mat_diffusion::*;
}

#[cfg(feature = "phase-transform")]
pub mod phase_transform {
    //! Re-export of `tpt-mat-phase-transform`.
    pub use tpt_mat_phase_transform::*;
}

#[cfg(feature = "calphad")]
pub mod calphad {
    //! Re-export of `tpt-mat-calphad`.
    pub use tpt_mat_calphad::*;
}

#[cfg(feature = "homogenization")]
pub mod homogenization {
    //! Re-export of `tpt-mat-homogenization`.
    pub use tpt_mat_homogenization::*;
}

#[cfg(feature = "rve")]
pub mod rve {
    //! Re-export of `tpt-mat-rve`.
    pub use tpt_mat_rve::*;
}

#[cfg(feature = "composite-micro")]
pub mod composite_micro {
    //! Re-export of `tpt-mat-composite-micro`.
    pub use tpt_mat_composite_micro::*;
}

#[cfg(feature = "damage")]
pub mod damage {
    //! Re-export of `tpt-mat-damage`.
    pub use tpt_mat_damage::*;
}

#[cfg(feature = "fatigue")]
pub mod fatigue {
    //! Re-export of `tpt-mat-fatigue`.
    pub use tpt_mat_fatigue::*;
}

#[cfg(feature = "creep")]
pub mod creep {
    //! Re-export of `tpt-mat-creep`.
    pub use tpt_mat_creep::*;
}

#[cfg(feature = "heat-treatment")]
pub mod heat_treatment {
    //! Re-export of `tpt-mat-heat-treatment`.
    pub use tpt_mat_heat_treatment::*;
}

#[cfg(feature = "database")]
pub mod database {
    //! Re-export of `tpt-mat-database`.
    pub use tpt_mat_database::*;
}

#[cfg(feature = "wasm")]
pub mod wasm_bindings {
    //! Re-export of `tpt-mat-wasm`.
    pub use tpt_mat_wasm::*;
}

/// Crate version string.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_non_empty() {
        assert!(!VERSION.is_empty());
    }
}