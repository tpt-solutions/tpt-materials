//! Miner (1945) linear damage accumulation.
//!
//! `D = Σ n_i / N_i`
//!
//! with failure at `D = 1` (Miner's rule).  The method is widely
//! used for variable-amplitude high-cycle fatigue: it sums the
//! partial damage from each stress cycle and predicts failure
//! when the cumulative damage reaches 1.

use serde::{Deserialize, Serialize};

/// A single `(n_i, N_i)` pair: `n_i` cycles at a stress amplitude
/// for which the fatigue-life curve predicts `N_i` cycles to
/// failure.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MinerCyclicInputs {
    /// Number of applied cycles at this stress level.
    pub n_i: f64,
    /// Number of cycles to failure at this stress level (from
    /// e.g. `BasquinParams::cycles_to_failure_basquin`).
    pub n_failure_i: f64,
}

impl MinerCyclicInputs {
    /// Damage contribution `n_i / N_i` from this single block.
    pub fn damage_contribution(&self) -> f64 {
        if self.n_failure_i <= 0.0 {
            return 0.0;
        }
        self.n_i / self.n_failure_i
    }
}

/// Cumulative Miner damage `Σ n_i / N_i` for the supplied cycle
/// blocks.
pub fn miner_damage_accumulation(blocks: &[MinerCyclicInputs]) -> f64 {
    blocks.iter().map(|b| b.damage_contribution()).sum()
}

/// Remaining-life fraction: how many *additional* cycles are
/// allowed at the *current* stress level before Miner damage
/// reaches 1.  Negative means the structure has already failed.
pub fn miner_remaining_life(
    blocks: &[MinerCyclicInputs],
    current_d_per_cycle: f64,
) -> f64 {
    if current_d_per_cycle <= 0.0 {
        return f64::INFINITY;
    }
    let d_current = miner_damage_accumulation(blocks);
    (1.0 - d_current) / current_d_per_cycle
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn miner_damage_zero_for_empty() {
        assert_eq!(miner_damage_accumulation(&[]), 0.0);
    }

    #[test]
    fn miner_damage_sums_individual_contributions() {
        let blocks = vec![
            MinerCyclicInputs { n_i: 100.0, n_failure_i: 1000.0 },
            MinerCyclicInputs { n_i: 500.0, n_failure_i: 1000.0 },
        ];
        let d = miner_damage_accumulation(&blocks);
        assert_relative_eq!(d, 0.1 + 0.5, epsilon = 1e-12);
    }

    #[test]
    fn miner_damage_ignores_blocks_with_zero_n_failure() {
        let blocks = vec![MinerCyclicInputs { n_i: 100.0, n_failure_i: 0.0 }];
        assert_eq!(miner_damage_accumulation(&blocks), 0.0);
    }

    #[test]
    fn remaining_life_decreases_with_existing_damage() {
        let blocks = vec![MinerCyclicInputs { n_i: 500.0, n_failure_i: 1000.0 }];
        let n_remaining = miner_remaining_life(&blocks, 0.001);
        // Damage 0.5 + 500 more cycles at 0.001/cycle ⇒ 1000 left.
        assert!((n_remaining - 500.0).abs() < 1e-6);
    }

    #[test]
    fn remaining_life_is_infinite_for_zero_rate() {
        let blocks = vec![MinerCyclicInputs { n_i: 100.0, n_failure_i: 1000.0 }];
        assert!(miner_remaining_life(&blocks, 0.0).is_infinite());
    }
}