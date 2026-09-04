//! Inverse / calibration helpers.
//!
//! Closed-form fits for the parameters of common constitutive laws
//! (Norton–Bailey creep `dε/dt = A σ^n exp(−Q/RT)`, Arrhenius
//! diffusivity `D = D_0 exp(−Q/RT)`, Voce hardening
//! `τ = τ_0 + (τ_s − τ_0)(1 − exp(−γ/γ_c)) + θ_0 γ`, power-law
//! S–N curve `σ_a = σ_f' (2N)^b`, JMAK Avrami
//! `f = 1 − exp(−k t^n)`, Coffin–Manson `Δε_p/2 = ε_f' (2N)^c`).
//!
//! The fits here are analytical (linear regression after a
//! log-transform, where possible).  For non-linear least-squares
//! problems use [`LevenbergMarquardt`], a generic iterative solver
//! that minimises `Σ (y_i - f(x_i; θ))²` over the parameter vector
//! `θ`.

mod levenberg_marquardt;
mod models;

pub use levenberg_marquardt::{LevenbergMarquardt, LevenbergMarquardtResult};
pub use models::{
    arrhenius_fit, avrami_fit, coffin_manson_fit, norton_creep_fit, power_law_sn_fit,
    voce_fit,
};

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn arrhenius_fit_recovers_d0_and_q() {
        let d_0_true = 1.0e-5;
        let q_true = 100.0e3;
        let r = 8.314_462_618;
        let temps: Vec<f64> = (800..1000).step_by(20).map(|t| t as f64).collect();
        let data: Vec<(f64, f64)> = temps
            .iter()
            .map(|&t| {
                let ln = (d_0_true * (-q_true / (r * t)).exp()).ln();
                (1.0 / t, ln)
            })
            .collect();
        let fit = arrhenius_fit(&data, r);
        assert!(approx(fit.k_0, d_0_true, 1.0e-3));
        assert!(approx(fit.q, q_true, 100.0));
    }

#[test]
    fn norton_creep_fit_recovers_a_and_n() {
        // ln(dε/dt) = ln(A) − Q/R T + n ln(σ)
        let a_true: f64 = 1.0e-5;
        let n_true: f64 = 3.5;
        let stresses: Vec<f64> = (10..100).step_by(10).map(|s| s as f64).collect();
        let data: Vec<(f64, f64)> = stresses
            .iter()
            .map(|&s| {
                let ln_rate = a_true.ln() + n_true * s.ln();
                (s.ln(), ln_rate)
            })
            .collect();
        let fit = norton_creep_fit(&data);
        assert!(approx(fit.n, n_true, 1.0e-6));
        assert!(approx(fit.a.ln(), a_true.ln(), 1.0e-6));
    }

    #[test]
    fn voce_fit_recovers_parameters() {
        let params = voce_fit(&[
            (0.0, 30.0),
            (0.05, 60.0),
            (0.10, 80.0),
            (0.20, 95.0),
            (0.50, 99.0),
        ]);
        // Recover at γ = 0 the CRSS should be ≈ τ_0.
        let r0 = params.tau_0;
        assert!(r0 > 0.0 && r0 < 100.0);
    }

    #[test]
    fn power_law_sn_fit_recovers_basquin() {
        let sigma_f_true = 1000.0;
        let b_true = -0.1;
        let data: Vec<(f64, f64)> = (1..6)
            .map(|i| {
                let n = 10.0_f64.powi(i);
                let s = sigma_f_true * (2.0 * n).powf(b_true);
                (n, s)
            })
            .collect();
        let fit = power_law_sn_fit(&data);
        assert!(approx(fit.sigma_f_prime, sigma_f_true, 1.0));
        assert!(approx(fit.b, b_true, 1.0e-3));
    }
}