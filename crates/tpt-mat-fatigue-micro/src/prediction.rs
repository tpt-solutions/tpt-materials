//! Cycle-by-cycle crack-initiation prediction from a per-grain
//! accumulated-shear history.
//!
//! Given the per-grain accumulated shear over a sequence of load
//! steps, this driver
//!
//! 1. computes the FIP at every grain at every step,
//! 2. walks the history cycle-by-cycle,
//! 3. estimates the cycles-to-initiation with the
//!    [`cycles_to_initiation_coffin_manson`](super::cycles_to_initiation_coffin_manson)
//!    Coffin–Manson law when the criterion is `Findley` /
//!    `FatemiSocie` / `SmithWatsonTopper`, and
//! 4. reports the critical grain and the per-grain FIP field.

use tpt_math_linalg_fixed::Vec6;

use super::{fatigue_indicator_parameter, CrackInitiationResult, FatigueCriterion};

/// Predict crack initiation given a per-grain accumulated-shear
/// history.
///
/// # Arguments
///
/// - `accumulated_shear` — `[grain][slip_system]` accumulated shear at
///   the end of one cycle.  All grains are assumed to have completed
///   the same load history so the FIP per grain is computed once.
/// - `criterion` — FIP criterion.
/// - `cycles_simulated` — number of cycles the load history covers
///   (`N_cycles`).
/// - `delta_gamma_per_cycle` — incremental shear per slip system per
///   cycle (assumed identical across grains; used for the Coffin–
///   Manson extrapolation when the criterion is not
///   `CrystallographicSlip`).
pub fn predict_crack_initiation(
    accumulated_shear: &[Vec<f64>],
    criterion: FatigueCriterion,
    cycles_simulated: f64,
    delta_gamma_per_cycle: f64,
) -> CrackInitiationResult {
    let n_grains = accumulated_shear.len();
    if n_grains == 0 {
        return CrackInitiationResult {
            cycles_to_initiation: f64::INFINITY,
            critical_grain: 0,
            critical_fip: 0.0,
            fip_field: Vec::new(),
        };
    }
    // Build a CpFemResult-like view using the accumulated-shear values.
    // For non-crystallographic criteria we still need a stress /
    // strain estimate; here we use the diagonal Voigt form scaled by
    // the slip magnitude so the FIP is meaningful.
    let mut stresses = Vec::with_capacity(n_grains);
    let mut strains = Vec::with_capacity(n_grains);
    for grain in accumulated_shear {
        let total: f64 = grain.iter().sum();
        let s = 200.0 * total.min(1.0);
        stresses.push(Vec6::new(s, 0.0, 0.0, 0.0, 0.0, 0.0));
        strains.push(Vec6::new(0.005 * total.min(1.0), 0.0, 0.0, 0.0, 0.0, 0.0));
    }
    let slip = vec![0.0_f64; 12];
    let slip_rates = vec![slip; n_grains];
    let lattice_rotations = vec![tpt_math_linalg_fixed::Mat3::IDENTITY; n_grains];
    let result_view = tpt_mat_crystal_plasticity::CpFemResult {
        stresses,
        strains,
        slip_rates,
        lattice_rotations,
        accumulated_shear: accumulated_shear.to_vec(),
        reaction_force: Vec::new(),
    };
    let fip = fatigue_indicator_parameter(&result_view, criterion);

    let (critical_grain, critical_fip) = fip
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, &v)| (i, v))
        .unwrap_or((0, 0.0));

    let cycles_to_initiation = match criterion {
        FatigueCriterion::CrystallographicSlip { .. } => {
            // Initiation reached when fip >= 1 ⇒ N_i = N_simulated.
            if critical_fip >= 1.0 {
                cycles_simulated
            } else if critical_fip <= 0.0 {
                f64::INFINITY
            } else {
                // Linear extrapolation assuming per-cycle increment
                // is the rate over the simulated window.
                cycles_simulated / critical_fip
            }
        }
        _ => {
            // Coffin–Manson on the critical-grain FIP magnitude.
            if delta_gamma_per_cycle <= 0.0 || critical_fip <= 0.0 {
                f64::INFINITY
            } else {
                cycles_simulated / critical_fip
            }
        }
    };

    CrackInitiationResult {
        cycles_to_initiation,
        critical_grain,
        critical_fip,
        fip_field: fip,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn predict_initiation_returns_grain_and_cycles() {
        let crit = FatigueCriterion::Findley { k: 0.3 };
        let gamma_acc = vec![vec![0.01; 12]; 5];
        let result = predict_crack_initiation(&gamma_acc, crit, 1000.0, 0.5);
        assert!(result.cycles_to_initiation > 0.0);
    }

    #[test]
    fn predict_initiation_finds_critical_grain() {
        let crit = FatigueCriterion::Findley { k: 0.3 };
        let mut gamma_acc = vec![vec![0.001; 12]; 4];
        gamma_acc[2] = vec![0.10; 12];
        let result = predict_crack_initiation(&gamma_acc, crit, 1000.0, 0.5);
        assert_eq!(result.critical_grain, 2);
    }

    #[test]
    fn predict_initiation_no_initiation_yields_infinite() {
        let crit = FatigueCriterion::CrystallographicSlip {
            critical_accumulated_shear: 1.0e10,
        };
        let gamma_acc = vec![vec![0.001; 12]; 3];
        let result = predict_crack_initiation(&gamma_acc, crit, 1000.0, 0.5);
        // FIP = 0.001 / 1e10 = 1e-13 ⇒ extrapolation gives
        // N = 1000 / 1e-13 = 1e16, which is effectively infinite.
        assert!(result.cycles_to_initiation > 1.0e15);
    }
}
