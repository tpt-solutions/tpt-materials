//! HEDE / HELP embrittlement indicators.

use serde::{Deserialize, Serialize};

/// HEDE (Hydrogen-Enhanced Decohesion) model parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct HedeParams {
    /// Critical hydrogen concentration at the fracture site
    /// (mol/m³).
    pub c_critical: f64,
}

/// HELP (Hydrogen-Enhanced Localised Plasticity) parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct HelpParams {
    /// Critical local hydrogen concentration for HELP activation.
    pub c_critical: f64,
    /// Critical stress triaxiality for HELP activation.
    pub triaxiality_critical: f64,
}

/// HEDE indicator: returns `true` if the local concentration
/// meets or exceeds the decohesion threshold.
pub fn hede_threshold(c_local: f64, params: &HedeParams) -> bool {
    c_local >= params.c_critical
}

/// HELP indicator: returns `true` if both the local H
/// concentration and the stress triaxiality meet their
/// respective thresholds.
pub fn help_threshold(c_local: f64, triaxiality: f64, params: &HelpParams) -> bool {
    c_local >= params.c_critical && triaxiality >= params.triaxiality_critical
}

/// Combined susceptibility index `C × T` (concentration times
/// triaxiality), normalised by the critical values.  Values
/// above 1.0 indicate high risk.
pub fn susceptibility_index(
    c_local: f64,
    triaxiality: f64,
    hede: &HedeParams,
    help: &HelpParams,
) -> f64 {
    let nc = (c_local / hede.c_critical).max(0.0);
    let nt = (triaxiality / help.triaxiality_critical).max(0.0);
    nc * nt
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn hede_above_threshold_returns_true() {
        let p = HedeParams { c_critical: 1.0 };
        assert!(hede_threshold(2.0, &p));
    }

    #[test]
    fn hede_below_threshold_returns_false() {
        let p = HedeParams { c_critical: 1.0 };
        assert!(!hede_threshold(0.5, &p));
    }

    #[test]
    fn help_requires_both_thresholds() {
        let p = HelpParams {
            c_critical: 1.0,
            triaxiality_critical: 0.5,
        };
        assert!(!help_threshold(2.0, 0.1, &p));
        assert!(!help_threshold(0.1, 2.0, &p));
        assert!(help_threshold(2.0, 2.0, &p));
    }

    #[test]
    fn susceptibility_index_above_one_indicates_risk() {
        let hede = HedeParams { c_critical: 1.0 };
        let help = HelpParams {
            c_critical: 1.0,
            triaxiality_critical: 0.5,
        };
        // c = 2 (2x), T = 1.0 (2x critical) -> 2 * 2 = 4.
        let s = susceptibility_index(2.0, 1.0, &hede, &help);
        assert!(approx(s, 4.0, 1.0e-9));
        assert!(s > 1.0);
    }
}