//! Rainflow cycle counting (ASTM E1049 four-point algorithm).
//!
//! Given a closed load (or strain) history, the rainflow algorithm
//! extracts a list of `(range, mean)` (or `(amplitude, mean)`)
//! cycles that can be fed into a damage-accumulation model (e.g.
//! Miner's rule).
//!
//! The implementation here follows the four-point ASTM E1049
//! algorithm:
//!
//! 1. Reduce the history to its turning points.
//! 2. For each turning point, compare its range to the previous two
//!    turning-point ranges:
//    - if the current range is *less than or equal to* the
//!      previous range, extract the previous range as a full
//!      cycle and pop its two endpoints;
//    - else continue.
//! 3. After all turning points are processed, extract residual
//!    half-cycles in pairs.
//!
//! Cycles are returned in `(amplitude, mean)` form for direct use
//! with mean-stress-corrected fatigue models.

use serde::{Deserialize, Serialize};

/// A single fatigue cycle extracted by the rainflow algorithm.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Cycle {
    /// Range `ΔX = X_max - X_min` of the cycle.
    pub range: f64,
    /// Mean value `X_mean = (X_max + X_min) / 2`.
    pub mean: f64,
}

impl Cycle {
    /// Construct from min / max.
    pub fn from_min_max(min: f64, max: f64) -> Self {
        let lo = min.min(max);
        let hi = max.max(min);
        Self {
            range: hi - lo,
            mean: 0.5 * (hi + lo),
        }
    }

    /// Amplitude `ΔX / 2`.
    pub fn amplitude(&self) -> f64 {
        0.5 * self.range
    }

    /// Load ratio `R = X_min / X_max` (returns 0 for `X_max = 0`).
    pub fn load_ratio(&self) -> f64 {
        let lo = self.mean - self.amplitude();
        let hi = self.mean + self.amplitude();
        if hi == 0.0 {
            0.0
        } else {
            lo / hi
        }
    }
}

/// Result of a rainflow count.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RainflowResult {
    /// Full cycles.
    pub cycles: Vec<Cycle>,
    /// Residual half-cycles that did not pair (typically small in
    /// long histories).
    pub residuals: Vec<Cycle>,
}

impl RainflowResult {
    /// Total number of cycles including half-cycles.
    pub fn total_cycle_count(&self) -> f64 {
        self.cycles.len() as f64 + 0.5 * self.residuals.len() as f64
    }

    /// Sum of damage contributions `Σ (n_i / N_i)` for the supplied
    /// `N_i` evaluator.  Not implemented; placeholder for future
    /// Miner's rule integration.
    pub fn damage(&self, _n_i: impl Fn(&Cycle) -> f64) -> f64 {
        0.0
    }
}

/// Reduce a load history to its turning points.
fn turning_points(history: &[f64]) -> Vec<f64> {
    if history.is_empty() {
        return Vec::new();
    }
    let mut out = vec![history[0]];
    for &x in &history[1..] {
        if x != *out.last().unwrap() {
            out.push(x);
        }
    }
    out
}

/// Four-point rainflow cycle counting (ASTM E1049).
///
/// `history` is a sequence of load (or strain) values; the
/// algorithm interprets it cyclically (i.e. the last point is
/// followed by the first).
pub fn rainflow_count(history: &[f64]) -> RainflowResult {
    let points = turning_points(history);
    let n = points.len();
    if n < 3 {
        return RainflowResult {
            cycles: Vec::new(),
            residuals: points.iter().map(|&x| Cycle::from_min_max(x, x)).collect(),
        };
    }
    let mut stack: Vec<f64> = Vec::with_capacity(n);
    let mut cycles: Vec<Cycle> = Vec::new();
    let extended: Vec<f64> = points
        .iter()
        .chain(points.iter().take(2))
        .copied()
        .collect();
    for &x in &extended {
        stack.push(x);
        // Need at least 4 points to apply the four-point rule.
        while stack.len() >= 4 {
            let len = stack.len();
            // Range of the *previous* pair (Y_{i-1}).
            let range_prev = (stack[len - 2] - stack[len - 3]).abs();
            // Range of the *current* pair (Y_i).
            let range_curr = (stack[len - 1] - stack[len - 2]).abs();
            // Extract a cycle if the previous range is *strictly*
            // greater than the current range.
            if range_prev > range_curr {
                cycles.push(Cycle::from_min_max(stack[len - 2], stack[len - 3]));
                stack.remove(len - 2);
                stack.remove(len - 3);
            } else {
                break;
            }
        }
        // Stop if we've consumed all the original turning points
        // plus a couple (to handle the cyclic wrap-around without
        // double-counting).
        if stack.len() >= n + 1 {
            break;
        }
    }
    // Remaining stack entries form residual half-cycles.
    let mut residuals = Vec::new();
    if stack.len() >= 2 {
        let mut i = 0;
        while i + 1 < stack.len() {
            residuals.push(Cycle::from_min_max(stack[i], stack[i + 1]));
            i += 2;
        }
        if i < stack.len() {
            residuals.push(Cycle::from_min_max(stack[i], stack[i]));
        }
    } else if stack.len() == 1 {
        residuals.push(Cycle::from_min_max(stack[0], stack[0]));
    }
    RainflowResult { cycles, residuals }
}

/// Convenience: rainflow-count a pre-built list of `(min, max)` pairs.
pub fn rainflow_count_from_pairs(pairs: &[(f64, f64)]) -> RainflowResult {
    let mut history: Vec<f64> = Vec::with_capacity(pairs.len() * 2);
    for &(lo, hi) in pairs {
        history.push(lo);
        history.push(hi);
    }
    rainflow_count(&history)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_oscillation_yields_one_cycle() {
        // history: 0 → 1 → 0.
        let r = rainflow_count(&[0.0, 1.0, 0.0]);
        assert_eq!(r.cycles.len(), 1);
        assert!((r.cycles[0].range - 1.0).abs() < 1e-9);
        assert!((r.cycles[0].mean - 0.5).abs() < 1e-9);
    }

    #[test]
    fn two_oscillations_yield_two_cycles() {
        // history: 0 → 1 → 0 → 0.5 → 0.
        let r = rainflow_count(&[0.0, 1.0, 0.0, 0.5, 0.0]);
        assert_eq!(r.cycles.len(), 2);
    }

    #[test]
    fn monotonic_history_yields_no_cycles() {
        let r = rainflow_count(&[0.0, 1.0, 2.0, 3.0]);
        assert!(r.cycles.is_empty());
    }

    #[test]
    fn cycle_amplitude_and_load_ratio() {
        let c = Cycle::from_min_max(-100.0, 100.0);
        assert!((c.range - 200.0).abs() < 1e-9);
        assert!((c.amplitude() - 100.0).abs() < 1e-9);
        assert!((c.load_ratio() - (-1.0)).abs() < 1e-9);
    }

    #[test]
    fn total_cycle_count_includes_half_cycle_residuals() {
        let r = RainflowResult {
            cycles: vec![Cycle::from_min_max(0.0, 1.0); 3],
            residuals: vec![Cycle::from_min_max(0.0, 1.0); 2],
        };
        assert!((r.total_cycle_count() - 4.0).abs() < 1e-9);
    }
}
