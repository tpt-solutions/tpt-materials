//! Generic non-linear least-squares via Levenberg–Marquardt.
//!
//! Solves `min_θ Σ_i r_i(θ)²` where the residuals are supplied as
//! a closure `r(θ) -> Vec<f64>`.  An optional analytical or
//! numerical Jacobian is used; if no Jacobian is supplied the
//! solver uses forward finite differences.
//!
//! This is a small, dependency-free Levenberg–Marquardt routine —
//! adequate for fitting constitutive-model parameters with `n ≤ 50`
//! parameters.

use serde::{Deserialize, Serialize};

/// Result of a Levenberg–Marquardt solve.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LevenbergMarquardtResult {
    /// Best-fit parameters.
    pub theta: Vec<f64>,
    /// Final sum of squared residuals.
    pub sse: f64,
    /// Number of iterations executed.
    pub iterations: usize,
    /// `true` if the iteration converged before the iteration cap.
    pub converged: bool,
}

/// Levenberg–Marquardt driver.
pub struct LevenbergMarquardt {
    /// Initial damping factor.
    pub initial_damping: f64,
    /// Damping scale factor on success/failure (typ. 10).
    pub damping_scale: f64,
    /// Convergence tolerance on the parameter step norm.
    pub param_tol: f64,
    /// Convergence tolerance on the SSE change.
    pub sse_tol: f64,
    /// Maximum iterations.
    pub max_iterations: usize,
}

impl Default for LevenbergMarquardt {
    fn default() -> Self {
        Self {
            initial_damping: 1.0e-3,
            damping_scale: 10.0,
            param_tol: 1.0e-8,
            sse_tol: 1.0e-12,
            max_iterations: 200,
        }
    }
}

impl LevenbergMarquardt {
    /// Construct with defaults.
    pub fn new() -> Self {
        Self::default()
    }

    /// Solve `min_θ Σ r_i(θ)²` for the given closure.  `r` returns
    /// the residual at `theta`; `theta_init` is the starting
    /// parameter vector.
    pub fn solve<F>(&self, theta_init: &[f64], r: F) -> LevenbergMarquardtResult
    where
        F: Fn(&[f64]) -> Vec<f64>,
    {
        let n = theta_init.len();
        let mut theta = theta_init.to_vec();
        let mut damping = self.initial_damping.max(1.0e-30);
        let mut sse = sum_sq(&r(&theta));
        let mut converged = false;
        let mut iterations = 0;
        let mut prev_sse = sse;
        for it in 0..self.max_iterations {
            iterations = it + 1;
            let residuals = r(&theta);
            let m = residuals.len();
            if m == 0 || n == 0 {
                break;
            }
            // Numerical Jacobian (forward differences).
            let mut j = vec![vec![0.0_f64; n]; m];
            let eps = 1.0e-7;
            for j_col in 0..n {
                let mut t_p = theta.clone();
                t_p[j_col] += eps * theta[j_col].abs().max(1.0e-3);
                let r_p = r(&t_p);
                for i in 0..m {
                    j[i][j_col] = (r_p[i] - residuals[i]) / (eps * theta[j_col].abs().max(1.0e-3));
                }
            }
            // J^T J + λ I.
            let mut jtj = vec![vec![0.0_f64; n]; n];
            for i in 0..n {
                for k in 0..n {
                    let mut acc = 0.0;
                    for row in 0..m {
                        acc += j[row][i] * j[row][k];
                    }
                    jtj[i][k] = acc;
                }
            }
            for i in 0..n {
                jtj[i][i] += damping * (1.0 + jtj[i][i].abs());
            }
            // J^T r.
            let mut jtr = vec![0.0_f64; n];
            for i in 0..n {
                let mut acc = 0.0;
                for row in 0..m {
                    acc += j[row][i] * residuals[row];
                }
                jtr[i] = acc;
            }
            // Solve (J^T J + λ I) δθ = -J^T r.
            let neg_jtr: Vec<f64> = jtr.iter().map(|v| -v).collect();
            let Some(delta) = solve_linear(&jtj, &neg_jtr) else {
                damping *= self.damping_scale;
                continue;
            };
            let mut theta_new = theta.clone();
            for i in 0..n {
                theta_new[i] += delta[i];
            }
            let new_sse = sum_sq(&r(&theta_new));
            if new_sse < sse {
                let ratio = if prev_sse > 0.0 {
                    (prev_sse - new_sse) / prev_sse
                } else {
                    0.0
                };
                theta = theta_new;
                prev_sse = sse;
                sse = new_sse;
                damping = (damping / self.damping_scale).max(1.0e-30);
                if ratio < self.sse_tol {
                    converged = true;
                    break;
                }
            } else {
                damping *= self.damping_scale;
                if damping > 1.0e10 {
                    break;
                }
            }
            let step_norm: f64 = delta.iter().map(|d| d * d).sum::<f64>().sqrt();
            if step_norm < self.param_tol {
                converged = true;
                break;
            }
        }
        LevenbergMarquardtResult {
            theta,
            sse,
            iterations,
            converged,
        }
    }
}

fn sum_sq(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum()
}

fn solve_linear(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = a.len();
    if b.len() != n {
        return None;
    }
    let mut m = vec![vec![0.0_f64; n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            m[i][j] = a[i][j];
        }
        m[i][n] = b[i];
    }
    for i in 0..n {
        let mut piv = i;
        for k in (i + 1)..n {
            if m[k][i].abs() > m[piv][i].abs() {
                piv = k;
            }
        }
        if m[piv][i].abs() < 1.0e-15 {
            return None;
        }
        m.swap(i, piv);
        let inv_piv = 1.0 / m[i][i];
        for j in 0..=n {
            m[i][j] *= inv_piv;
        }
        for k in 0..n {
            if k != i {
                let f = m[k][i];
                for j in 0..=n {
                    m[k][j] -= f * m[i][j];
                }
            }
        }
    }
    Some((0..n).map(|i| m[i][n]).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recovers_a_and_b_in_y_a_x_plus_b() {
        let lm = LevenbergMarquardt::new();
        let result = lm.solve(&[0.0_f64, 0.0], |theta| {
            let a = theta[0];
            let b = theta[1];
            let xs = [0.0, 1.0, 2.0, 3.0, 4.0];
            xs.iter()
                .map(|&x| ((3.0 * x + 2.0) - (a * x + b)))
                .collect()
        });
        assert!(result.converged);
        assert!((result.theta[0] - 3.0).abs() < 1.0e-4);
        assert!((result.theta[1] - 2.0).abs() < 1.0e-4);
    }

    #[test]
    fn recovers_exponential_decay_rate() {
        let lm = LevenbergMarquardt::new();
        let result = lm.solve(&[1.0_f64, 0.0], |theta| {
            let (a, k) = (theta[0], theta[1]);
            let ts = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0];
            ts.iter()
                .map(|&t| ((5.0_f64 * (-1.5_f64 * t).exp()) - (a * (-k * t).exp())))
                .collect()
        });
        assert!(result.converged);
        assert!((result.theta[0] - 5.0).abs() < 1.0e-3);
        assert!((result.theta[1] - 1.5).abs() < 1.0e-3);
    }
}