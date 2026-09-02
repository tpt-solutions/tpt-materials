//! Latent-hardening matrix `h_{αβ}` (cross-hardening between slip
//! systems).
//!
//! Latent hardening is the increase of CRSS on system `α` due to slip
//! activity on system `β` (`α ≠ β`).  The standard
//! [`LatentHardeningMatrix::diagonal`] returns a matrix that only
//! self-hardens; [`LatentHardeningMatrix::with_latent`] returns
//! `h_{αβ} = h_0` if `α = β` and `h_{αβ} = h_0 q` if `α ≠ β`, where
//! `q ≥ 1` is the latent-to-self ratio (Taylor, 1934; Kocks & Brown,
//! 1966).

use serde::{Deserialize, Serialize};

/// `n_slip × n_slip` hardening matrix.
///
/// Stored row-major for cache-friendly inner-loop access during the
/// `Δτ_c^α = Σ_β h_{αβ} Δγ^β` accumulation in
/// [`CombinedHardening`](crate::CombinedHardening).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LatentHardeningMatrix {
    /// Flat row-major data, length `n × n`.
    pub data: Vec<f64>,
    /// Dimension.
    pub n: usize,
}

impl LatentHardeningMatrix {
    /// Diagonal matrix `h_{αβ} = h_0 δ_{αβ}` (no cross-hardening) for
    /// `n` slip systems.
    pub fn diagonal(n: usize, h_0: f64) -> Self {
        let mut data = vec![0.0; n * n];
        for i in 0..n {
            data[i * n + i] = h_0;
        }
        Self { data, n }
    }

    /// Latent-hardening matrix `h_{αβ} = h_0 (q for α ≠ β, 1 for α = β)`
    /// with the same dimension as `self`.
    pub fn with_latent(&self, q: f64, h_0: f64) -> Self {
        let n = self.n;
        let mut data = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                data[i * n + j] = if i == j { h_0 } else { h_0 * q };
            }
        }
        Self { data, n }
    }

    /// Entry `h_{αβ}`.
    #[inline]
    pub fn get(&self, alpha: usize, beta: usize) -> f64 {
        self.data[alpha * self.n + beta]
    }

    /// Set entry `h_{αβ}`.
    #[inline]
    pub fn set(&mut self, alpha: usize, beta: usize, value: f64) {
        self.data[alpha * self.n + beta] = value;
    }

    /// Construct an `n × n` zero matrix (no self- or latent-hardening).
    pub fn zero(n: usize) -> Self {
        Self {
            data: vec![0.0; n * n],
            n,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diagonal_has_only_self_terms() {
        let m = LatentHardeningMatrix::diagonal(3, 100.0);
        assert_eq!(m.n, 3);
        assert_eq!(m.get(0, 0), 100.0);
        assert_eq!(m.get(1, 1), 100.0);
        assert_eq!(m.get(2, 2), 100.0);
        assert_eq!(m.get(0, 1), 0.0);
    }

    #[test]
    fn with_latent_sets_off_diagonal() {
        let m = LatentHardeningMatrix::diagonal(3, 50.0).with_latent(1.4, 50.0);
        for i in 0..m.n {
            for j in 0..m.n {
                let expected = if i == j { 50.0 } else { 70.0 };
                assert_eq!(m.get(i, j), expected);
            }
        }
    }
}