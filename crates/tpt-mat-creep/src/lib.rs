//! Creep deformation and rupture models.
//!
//! This crate implements the classical analytical models for
//! time-dependent inelastic deformation under constant stress and
//! temperature:
//!
//! - **Norton–Bailey (1929/1934) power-law creep**: the steady-state
//!   creep rate `ε̇_ss = A σ^n exp(-Q / (R T))`.  For `n ∈ [3, 8]`
//!   the law is the canonical `power-law creep` regime; for
//!   `n > 8` it enters the `power-law breakdown` regime.
//! - **θ-projection (Wilshire & Burt, 1982)**: a four-parameter
//!   description of the *full* creep curve,
//!   `ε(t) = θ_1 (1 - exp(-θ_2 t)) + θ_3 (exp(θ_4 t) - 1)`.  The
//!   first term captures primary (decelerating) creep, the second
//!   tertiary (accelerating) creep.
//! - **Monkman–Grant (1956)**: an empirical relation between
//!   steady-state creep rate and rupture time,
//!   `t_r · ε̇_ss^m = C_mg` (typically `m ∈ [0.7, 1.0]`).
//! - **Larson–Miller (1952) parameter**: the time-temperature
//!   parameter `LMP = T (log t_r + C_LM)` (typically `C_LM ≈ 20`)
//!   used to extrapolate rupture data across temperatures.
//! - **Sherby–Dorn (1958) parameter**: a temperature-compensated
//!   time `θ_D = t_r · exp(-Q / (R T))`.
//!
//! # Conventions
//!
//! - Stress in MPa, strain dimensionless, time in seconds,
//!   temperature in Kelvin.
//! - The Norton–Bailey law here is the *secondary* (steady-state)
//!   creep rate, not a transient; for primary / tertiary use
//!   `theta_projection_creep_strain`.
//!
//! # References
//!
//! - Norton, F. H. (1929).  "The creep of steel at high
//!   temperatures."  Trans. ASME 51, 139–146.
//! - Bailey, R. W. (1934).  "The utilisation of creep test data
//!   in engineering design."  Proc. IME 131, 131–349.
//! - Monkman, F. C. & Grant, N. J. (1956).  "An empirical
//!   relationship between rupture life and minimum creep rate in
//!   creep-rupture tests."  Proc. ASTM 56, 593–620.
//! - Larson, F. R. & Miller, J. (1952).  "A time-temperature
//!   relationship for rupture and creep stresses."  Trans. ASME
//!   74, 765–775.
//! - Wilshire, B. & Burt, H. (1982).  "The θ-projection
//!   creep-constitutive equation."  Proc. 2nd Int. Conf. on
//!   Creep and Fracture of Engineering Materials and Structures.

#![warn(missing_docs)]

mod larson_miller;
mod monkman_grant;
mod norton;
mod sherby_dorn;
mod theta_projection;

pub use larson_miller::{larson_miller_parameter, rupture_time_from_lmp, LarsonMillerParams};
pub use monkman_grant::{monkman_grant_product, rupture_time_monkman_grant, MonkmanGrantParams};
pub use norton::{norton_creep_rate, rupture_time_norton, NortonBaileyParams};
pub use sherby_dorn::{rupture_time_from_dorn, sherby_dorn_parameter, SherbyDornParams};
pub use theta_projection::{
    theta_projection_creep_rate, theta_projection_creep_strain, ThetaProjectionParams,
};
