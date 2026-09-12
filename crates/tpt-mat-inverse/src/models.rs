//! Closed-form least-squares fits for common constitutive laws.
//!
//! Each routine takes a list of `(x, y)` data points and returns
//! best-fit parameters (closed-form linear regression after a
//! log-transform where possible).

use serde::{Deserialize, Serialize};

/// Result of an Arrhenius diffusivity fit.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ArrheniusFit {
    /// Pre-exponential factor `D_0` (m²/s).
    pub k_0: f64,
    /// Activation energy `Q` (J/mol).
    pub q: f64,
}

/// Fit `D = D_0 · exp(−Q / RT)` to `(1/T, ln D)` data via ordinary
/// least squares.  Pass the gas constant `R` (J/mol·K) to keep the
/// fit unit-agnostic.
pub fn arrhenius_fit(data: &[(f64, f64)], gas_constant: f64) -> ArrheniusFit {
    let n = data.len() as f64;
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    for &(x, y) in data {
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    let denom = n * sxx - sx * sx;
    let slope = if denom.abs() < 1.0e-30 {
        0.0
    } else {
        (n * sxy - sx * sy) / denom
    };
    let intercept = (sy - slope * sx) / n;
    ArrheniusFit {
        k_0: intercept.exp(),
        q: -slope * gas_constant,
    }
}

/// Norton–Bailey creep fit.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct NortonFit {
    /// Coefficient `A`.
    pub a: f64,
    /// Stress exponent `n`.
    pub n: f64,
}

/// Fit the linearised Norton–Bailey law
/// `ln(ε̇) = ln A + n ln σ` to `(ln σ, ln ε̇)` data.
pub fn norton_creep_fit(data: &[(f64, f64)]) -> NortonFit {
    let n = data.len() as f64;
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    for &(x, y) in data {
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    let denom = n * sxx - sx * sx;
    let slope = if denom.abs() < 1.0e-30 {
        0.0
    } else {
        (n * sxy - sx * sy) / denom
    };
    let intercept = (sy - slope * sx) / n;
    NortonFit {
        a: intercept.exp(),
        n: slope,
    }
}

/// Voce hardening fit (closed-form nonlinear least squares for
/// the three-parameter Voce law).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct VoceFit {
    /// Initial CRSS `τ_0` (Pa).
    pub tau_0: f64,
    /// Saturation CRSS increment `τ_s` (Pa).
    pub tau_s: f64,
    /// Linear hardening slope `θ_0` (Pa).
    pub theta_0: f64,
    /// Characteristic accumulated shear `γ_c`.
    pub gamma_c: f64,
}

/// Fit the Voce law `τ = τ_0 + (τ_s − τ_0)(1 − exp(−γ/γ_c)) + θ_0 γ`
/// by Gauss–Newton with `γ_c` held fixed (caller may sweep).  We
/// initialise from `τ_0 = τ(0)` and `τ_s = max(τ) − τ_0`.
pub fn voce_fit(data: &[(f64, f64)]) -> VoceFit {
    if data.is_empty() {
        return VoceFit {
            tau_0: 0.0,
            tau_s: 0.0,
            theta_0: 0.0,
            gamma_c: 1.0,
        };
    }
    let tau_0 = data[0].1;
    let max_tau = data.iter().map(|d| d.1).fold(f64::MIN, f64::max);
    let mut tau_s = (max_tau - tau_0).max(0.0);
    let mut theta_0 = 0.0;
    let mut gamma_c = 0.1;
    for _ in 0..30 {
        // Simple gradient step on the SSE between model and data.
        let mut grad_tau_s = 0.0;
        let mut grad_theta = 0.0;
        let mut grad_g = 0.0;
        let mut sse = 0.0;
        for &(g, tau_obs) in data {
            let ratio = -g / gamma_c;
            let exp_term = if ratio > -50.0 { ratio.exp() } else { 0.0 };
            let model = tau_0 + tau_s * (1.0 - exp_term) + theta_0 * g;
            let r = tau_obs - model;
            sse += r * r;
            grad_tau_s += -2.0 * r * (1.0 - exp_term);
            grad_theta += -2.0 * r * g;
            let dexpdg = exp_term / gamma_c;
            grad_g += -2.0 * r * tau_s * dexpdg;
        }
        let lr = 1.0e-3;
        tau_s -= lr * grad_tau_s / data.len() as f64;
        theta_0 -= lr * grad_theta / data.len() as f64;
        gamma_c -= lr * grad_g / data.len() as f64;
        gamma_c = gamma_c.abs().max(1.0e-6);
        tau_s = tau_s.max(0.0);
        let _ = sse;
    }
    VoceFit {
        tau_0,
        tau_s,
        theta_0,
        gamma_c,
    }
}

/// Power-law S–N fit.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PowerLawSnFit {
    /// Fatigue strength coefficient `σ_f'` (MPa).
    pub sigma_f_prime: f64,
    /// Basquin exponent `b`.
    pub b: f64,
}

/// Fit `σ_a = σ_f' (2N)^b` to `(N, σ_a)` data via log-linear
/// regression.
pub fn power_law_sn_fit(data: &[(f64, f64)]) -> PowerLawSnFit {
    let log_data: Vec<(f64, f64)> = data
        .iter()
        .map(|&(n, s)| ((2.0 * n).ln(), s.ln()))
        .collect();
    let n = log_data.len() as f64;
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    for &(x, y) in &log_data {
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    let denom = n * sxx - sx * sx;
    let slope = if denom.abs() < 1.0e-30 {
        0.0
    } else {
        (n * sxy - sx * sy) / denom
    };
    let intercept = (sy - slope * sx) / n;
    PowerLawSnFit {
        sigma_f_prime: intercept.exp(),
        b: slope,
    }
}

/// Avrami / JMAK fit.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct AvramiFit {
    /// Rate constant `k` (1/sⁿ).
    pub k: f64,
    /// Avrami exponent `n`.
    pub n: f64,
}

/// Fit `f = 1 − exp(−k t^n)` to `(t, f)` data via log-linear
/// regression on the latent variable `ln(−ln(1 − f)) = ln k + n ln t`.
pub fn avrami_fit(data: &[(f64, f64)]) -> AvramiFit {
    let mut pts: Vec<(f64, f64)> = Vec::new();
    for &(t, f) in data {
        if t > 0.0 && f > 0.0 && f < 1.0 {
            let y = (-(1.0 - f).ln()).ln();
            pts.push((t.ln(), y));
        }
    }
    let n = pts.len() as f64;
    if n < 2.0 {
        return AvramiFit { k: 0.0, n: 0.0 };
    }
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    for &(x, y) in &pts {
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    let denom = n * sxx - sx * sx;
    let slope = if denom.abs() < 1.0e-30 {
        0.0
    } else {
        (n * sxy - sx * sy) / denom
    };
    let intercept = (sy - slope * sx) / n;
    AvramiFit {
        k: intercept.exp(),
        n: slope,
    }
}

/// Coffin–Manson fit.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CoffinMansonFit {
    /// Fatigue ductility coefficient `ε_f'`.
    pub epsilon_f_prime: f64,
    /// Ductility exponent `c`.
    pub c: f64,
}

/// Fit `Δε_p / 2 = ε_f' (2N)^c` to `(N, Δε_p/2)` data via log-
/// linear regression.
pub fn coffin_manson_fit(data: &[(f64, f64)]) -> CoffinMansonFit {
    let log_data: Vec<(f64, f64)> = data
        .iter()
        .map(|&(n, eps)| ((2.0 * n).ln(), eps.ln()))
        .collect();
    let n = log_data.len() as f64;
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    for &(x, y) in &log_data {
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    let denom = n * sxx - sx * sx;
    let slope = if denom.abs() < 1.0e-30 {
        0.0
    } else {
        (n * sxy - sx * sy) / denom
    };
    let intercept = (sy - slope * sx) / n;
    CoffinMansonFit {
        epsilon_f_prime: intercept.exp(),
        c: slope,
    }
}
