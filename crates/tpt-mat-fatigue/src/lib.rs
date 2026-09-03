//! High- and low-cycle fatigue models.
//!
//! This crate delivers the classical textbook fatigue-life and
//! crack-growth models required by Phase 6:
//!
//! - **S-N curves (Basquin 1910)**: high-cycle fatigue
//!   `σ_a = σ_f' (2 N_f)^b`, where `σ_f'` is the fatigue strength
//!   coefficient and `b` is the fatigue strength exponent (typically
//!   `b ∈ [-0.12, -0.05]` for steels).
//! - **Coffin–Manson (low-cycle fatigue)**: the plastic strain
//!   amplitude `Δε_p / 2 = ε_f' (2 N_f)^c`, where `ε_f'` is the
//!   fatigue ductility coefficient and `c ∈ [-0.7, -0.5]`.
//! - **Strain-life (Manson–Coffin–Basquin)** combination: the total
//!   strain amplitude `Δε / 2 = (σ_f' / E) (2 N_f)^b + ε_f' (2 N_f)^c`.
//! - **Walker mean-stress correction**: `σ_ar = σ_a / (1 - R)^γ`
//!   for the equivalent fully-reversed stress amplitude.
//! - **Paris law (crack growth)**: `da/dN = C (ΔK)^m` with the
//!   Forman / Walker / Nasgro variants for the threshold and
//!   critical stress-intensity range.
//! - **Rainflow counting**: extract cycles from a strain / load
//!   history using the classical ASTM E1049 four-point algorithm.
//!
//! # Conventions
//!
//! - `R = σ_min / σ_max` is the load ratio.  `R = -1` is fully
//!   reversed; `R = 0` is zero-to-tension; `R → 1` is mean-loaded.
//! - All cyclic quantities are *amplitudes* (half the peak-to-peak
//!   range) unless explicitly stated as a range (Δ).
//!
//! # References
//!
//! - Basquin, O. H. (1910).  "The exponential law of endurance
//!   tests."  Proc. ASTM 10, 625–630.
//! - Coffin, L. F. (1954).  "A study of the effects of cyclic
//!   thermal stresses on a ductile metal."  Trans. ASME 76,
//!   931–950.
//! - Manson, S. S. (1953).  "Behavior of materials under thermal
//!   stress."  NACA TN-2933.
//! - Paris, P. C. & Erdogan, F. (1963).  "A critical analysis of
//!   crack propagation laws."  Trans. ASME 85, 528–534.
//! - Walker, K. (1970).  "The effect of stress ratio during crack
//!   growth and fatigue."  ASTM STP 462, 1–14.

#![warn(missing_docs)]

mod basquin;
mod coffin_manson;
mod paris;
mod rainflow;
mod strain_life;
mod walker;

pub use basquin::{cycles_to_failure_basquin, fatigue_strength_at, BasquinParams};
pub use coffin_manson::{cycles_to_failure_coffin_manson, plastic_strain_at, CoffinMansonParams};
pub use paris::{crack_growth_rate, critical_crack_length, paris_lifetime, ParisParams};
pub use rainflow::{rainflow_count, rainflow_count_from_pairs, Cycle, RainflowResult};
pub use strain_life::{cycles_to_failure_strain_life, StrainLifeParams};
pub use walker::{equivalent_amplitude_walker, walker_gamma};
