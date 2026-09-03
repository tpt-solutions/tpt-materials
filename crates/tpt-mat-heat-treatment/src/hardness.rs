//! Rule-of-mixtures and Martensite–Carbon hardness models.

use serde::{Deserialize, Serialize};

/// Representative Vickers hardness values for the major
/// microstructural constituents of steels.
///
/// These are order-of-magnitude calibration values used by
/// the Maynier-style regression tradition
/// (`HV = Σ f_i HV_i + ΔHV_age`).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PhaseHardness {
    /// Martensite Vickers hardness.
    pub martensite: f64,
    /// Bainite Vickers hardness.
    pub bainite: f64,
    /// Pearlite Vickers hardness.
    pub pearlite: f64,
    /// Ferrite Vickers hardness.
    pub ferrite: f64,
    /// Retained austenite Vickers hardness.
    pub austenite: f64,
}

impl Default for PhaseHardness {
    fn default() -> Self {
        Self {
            martensite: 800.0,
            bainite: 450.0,
            pearlite: 250.0,
            ferrite: 100.0,
            austenite: 200.0,
        }
    }
}

/// Compute the rule-of-mixtures hardness from a vector of
/// phase fractions and the per-phase hardness table.
///
/// Phase fractions are normalised before use, so callers may
/// pass `Σ f_i ≈ 1` without strict equality.
///
/// # Arguments
///
/// - `fractions` — `&[f_martensite, f_bainite, f_pearlite, f_ferrite, f_austenite]`
/// - `table`     — per-phase hardness values
/// - `extra`     — additional contribution from precipitation
///   hardening (aging) in HV
pub fn hardness_from_fractions(fractions: &[f64; 5], table: &PhaseHardness, extra: f64) -> f64 {
    let sum: f64 = fractions.iter().sum();
    if sum <= 0.0 {
        return table.ferrite;
    }
    let hvs = [
        table.martensite,
        table.bainite,
        table.pearlite,
        table.ferrite,
        table.austenite,
    ];
    let mut weighted = 0.0;
    for (f, h) in fractions.iter().zip(hvs.iter()) {
        weighted += (f / sum) * h;
    }
    weighted + extra.max(0.0)
}

/// Martensite-Carbon hardness regression.
///
/// `HV_M ≈ 127 + 949 (wt% C) + 49 (wt% C)^2`
///
/// This is the classic Equation 6 from Krauss (2015) and
/// recovers `HV ≈ 925` at `C = 0.8 wt%`.
///
/// # Arguments
///
/// - `carbon_wt_pct` — carbon content in weight percent
pub fn hardness_martensite(carbon_wt_pct: f64) -> f64 {
    127.0 + 949.0 * carbon_wt_pct + 49.0 * carbon_wt_pct * carbon_wt_pct
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn rule_of_mixtures_returns_weighted_sum() {
        let t = PhaseHardness::default();
        // 50% martensite + 50% ferrite -> (800+100)/2 = 450.
        let f = [0.5, 0.0, 0.0, 0.5, 0.0];
        let h = hardness_from_fractions(&f, &t, 0.0);
        assert!(approx(h, 450.0, 1.0e-6));
    }

    #[test]
    fn rule_of_mixtures_adds_age_hardening() {
        let t = PhaseHardness::default();
        let f = [1.0, 0.0, 0.0, 0.0, 0.0];
        let h0 = hardness_from_fractions(&f, &t, 0.0);
        let h_age = hardness_from_fractions(&f, &t, 50.0);
        assert!(approx(h_age - h0, 50.0, 1.0e-6));
    }

    #[test]
    fn rule_of_mixtures_normalises_unnormalised_fractions() {
        let t = PhaseHardness::default();
        let f = [1.0, 0.0, 0.0, 1.0, 0.0];
        let h = hardness_from_fractions(&f, &t, 0.0);
        // Should be same as [0.5, 0, 0, 0.5, 0].
        assert!(approx(h, 450.0, 1.0e-6));
    }

    #[test]
    fn martensite_carbon_regression_recovers_classic_value() {
        // 0.4 wt% C -> HV ≈ 514.
        let hv = hardness_martensite(0.4);
        assert!(approx(hv, 514.8, 0.5));
        // 0.8 wt% C -> HV ≈ 917.
        let hv = hardness_martensite(0.8);
        assert!(approx(hv, 917.6, 0.5));
    }

    #[test]
    fn martensite_carbon_regression_is_monotonic() {
        let mut last = 0.0;
        for c in [0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0] {
            let hv = hardness_martensite(c);
            assert!(hv > last, "non-monotonic at c={c}");
            last = hv;
        }
    }
}
