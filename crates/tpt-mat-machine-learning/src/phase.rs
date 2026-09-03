//! Multi-class softmax phase predictor.

use serde::{Deserialize, Serialize};

/// Soft-class prediction: probabilities over a discrete set of phase
/// labels plus the most likely label.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PhasePrediction {
    /// Probabilities (one per phase label).
    pub probabilities: Vec<f64>,
    /// Index of the most likely phase.
    pub most_likely: usize,
}

/// Multi-class softmax phase predictor.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PhasePredictor {
    /// Per-class weight matrix (`n_classes × n_features`).
    pub weights: Vec<Vec<f64>>,
    /// Per-class bias.
    pub bias: Vec<f64>,
    /// Phase labels (one per class).
    pub labels: Vec<String>,
}

impl PhasePredictor {
    /// Softmax activation `p_i = exp(z_i) / Σ_j exp(z_j)`.
    pub fn predict(&self, x: &[f64]) -> PhasePrediction {
        let n_classes = self.weights.len();
        let n_features = self.weights[0].len();
        debug_assert_eq!(x.len(), n_features);
        let mut z = vec![0.0; n_classes];
        for i in 0..n_classes {
            let mut acc = self.bias[i];
            for j in 0..n_features {
                acc += self.weights[i][j] * x[j];
            }
            z[i] = acc;
        }
        let max_z = z.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let mut exp_z = vec![0.0; n_classes];
        for i in 0..n_classes {
            exp_z[i] = (z[i] - max_z).exp();
        }
        let sum: f64 = exp_z.iter().sum();
        for v in exp_z.iter_mut() {
            *v /= sum;
        }
        let mut most_likely = 0;
        let mut best = f64::NEG_INFINITY;
        for (i, &p) in exp_z.iter().enumerate() {
            if p > best {
                best = p;
                most_likely = i;
            }
        }
        PhasePrediction {
            probabilities: exp_z,
            most_likely,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn softmax_probabilities_sum_to_one() {
        let p = PhasePredictor {
            weights: vec![vec![1.0, 0.5], vec![0.5, 1.0]],
            bias: vec![0.0, 0.0],
            labels: vec!["FCC".to_string(), "BCC".to_string()],
        };
        let pred = p.predict(&[0.5, 0.5]);
        let sum: f64 = pred.probabilities.iter().sum();
        assert!(approx(sum, 1.0, 1.0e-9));
    }

    #[test]
    fn softmax_picks_max_weight() {
        let p = PhasePredictor {
            weights: vec![vec![0.0, 0.0], vec![10.0, 10.0]],
            bias: vec![0.0, 0.0],
            labels: vec!["a".to_string(), "b".to_string()],
        };
        let pred = p.predict(&[1.0, 1.0]);
        assert_eq!(pred.most_likely, 1);
    }
}
