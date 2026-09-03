//! Polarization curves: scan the potential axis and return the net
//! current density `i(E) = i_a(E) + i_c(E)` (with sign convention
//! anodic-positive).

use serde::{Deserialize, Serialize};

use crate::butler_volmer::butler_volmer_current_density;
use crate::rate::CorrosionModel;

/// Polarization branch identifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum PolarizationBranch {
    /// Net current (anodic − |cathodic|).
    Net,
    /// Anodic current only.
    Anodic,
    /// Cathodic current only (magnitude).
    Cathodic,
}

/// Sampled polarization curve.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PolarizationCurve {
    /// Potential values (V vs SHE).
    pub potentials: Vec<f64>,
    /// Current densities (A/m²), one per branch (matching
    /// `branch` selector).
    pub currents: Vec<f64>,
    /// Which branch the `currents` field contains.
    pub branch: PolarizationBranch,
}

/// Sample the polarization curve over `potential_range = (E_min, E_max)`.
///
/// `num_points` — number of samples (≥ 2).
pub fn polarization_curve(
    model: &CorrosionModel,
    potential_range: (f64, f64),
    num_points: usize,
    branch: PolarizationBranch,
) -> PolarizationCurve {
    assert!(num_points >= 2);
    let (e_min, e_max) = potential_range;
    let n = num_points;
    let mut potentials = Vec::with_capacity(n);
    let mut currents = Vec::with_capacity(n);
    for i in 0..n {
        let frac = i as f64 / (n - 1) as f64;
        let e = e_min + frac * (e_max - e_min);
        let i_a = butler_volmer_current_density(&model.anode, e);
        let i_c = butler_volmer_current_density(&model.cathode, e);
        let i_out = match branch {
            PolarizationBranch::Net => i_a + i_c,
            PolarizationBranch::Anodic => i_a,
            PolarizationBranch::Cathodic => -i_c,
        };
        potentials.push(e);
        currents.push(i_out);
    }
    PolarizationCurve {
        potentials,
        currents,
        branch,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::butler_volmer::ElectrodeKinetics;
    use crate::rate::CorrosionModel;

    fn fe_model() -> CorrosionModel {
        let anode = ElectrodeKinetics::from_alphas(-0.44, 1.0e-3, 0.5, 0.5, 2.0, 298.0);
        let cathode = ElectrodeKinetics::from_alphas(-0.41, 1.0e-1, 0.5, 0.5, 2.0, 298.0);
        CorrosionModel::new(anode, cathode, 7.0, 298.0)
    }

    #[test]
    fn polarization_curve_returns_n_points() {
        let pc = polarization_curve(&fe_model(), (-0.8, 0.0), 50, PolarizationBranch::Net);
        assert_eq!(pc.potentials.len(), 50);
        assert_eq!(pc.currents.len(), 50);
    }

    #[test]
    fn polarization_curve_anodic_branch_positive_far_anode() {
        let pc = polarization_curve(&fe_model(), (-0.8, 0.5), 100, PolarizationBranch::Anodic);
        // Far above E_eq(anode) the anodic branch dominates.
        let last = pc.currents.last().copied().unwrap();
        assert!(last > 0.0);
    }

    #[test]
    fn polarization_curve_cathodic_branch_positive_far_cathode() {
        let pc = polarization_curve(&fe_model(), (-1.0, 0.0), 100, PolarizationBranch::Cathodic);
        let first = pc.currents.first().copied().unwrap();
        assert!(first > 0.0);
    }

    #[test]
    fn polarization_curve_potentials_monotonic() {
        let pc = polarization_curve(&fe_model(), (-0.5, 0.0), 20, PolarizationBranch::Net);
        for i in 1..pc.potentials.len() {
            assert!(pc.potentials[i] > pc.potentials[i - 1]);
        }
    }
}
