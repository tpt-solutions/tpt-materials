//! Voce, power-law, and Kocks-Mecking self-hardening parameters.

use serde::{Deserialize, Serialize};

/// Voce saturating hardening parameters (per slip system, but the
/// `VoceHardening` law applies them uniformly).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct VoceParams {
    /// Initial CRSS `τ_0` (MPa).
    pub tau_0: f64,
    /// Saturation CRSS `τ_s` (MPa).
    pub tau_s: f64,
    /// Initial linear hardening slope `θ_0` (MPa).
    pub theta_0: f64,
    /// Reference strain `γ_c` controlling saturation rate.
    pub gamma_c: f64,
}

impl Default for VoceParams {
    fn default() -> Self {
        Self {
            tau_0: 10.0,
            tau_s: 100.0,
            theta_0: 500.0,
            gamma_c: 0.05,
        }
    }
}

/// Power-law (Hutchinson) rate-dependent hardening parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PowerLawParams {
    /// Reference hardening modulus `h_0` (MPa).
    pub h_0: f64,
    /// Saturation CRSS `τ_s` (MPa).
    pub tau_s: f64,
}

impl Default for PowerLawParams {
    fn default() -> Self {
        Self {
            h_0: 1000.0,
            tau_s: 50.0,
        }
    }
}

/// Kocks-Mecking evolution parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct KocksMeckingParams {
    /// Initial CRSS `τ_0` (MPa).
    pub tau_0: f64,
    /// Saturation CRSS `τ_s` (MPa).
    pub tau_s: f64,
    /// Stage-II hardening slope `θ_0` (MPa).
    pub theta_0: f64,
}

impl Default for KocksMeckingParams {
    fn default() -> Self {
        Self {
            tau_0: 20.0,
            tau_s: 80.0,
            theta_0: 500.0,
        }
    }
}

/// Combined law: an underlying self-hardening modulus `h_0` plus a
/// latent-hardening factor `q_latent` (`h_{αβ} = h_0 (q_latent for α≠β, 1 for α=β)`).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CombinedParams {
    /// Self-hardening modulus `h_0` (MPa).
    pub h_0: f64,
    /// Latent / self ratio (`q ≥ 1` typically).
    pub q_latent: f64,
}

impl Default for CombinedParams {
    fn default() -> Self {
        Self {
            h_0: 1000.0,
            q_latent: 1.4,
        }
    }
}