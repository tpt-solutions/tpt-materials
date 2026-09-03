//! Kernel ridge regression (RBF or polynomial kernel).

use serde::{Deserialize, Serialize};

use super::TrainingResult;

/// Kernel ridge regression model.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct KernelRidgeModel {
    /// Training inputs (support points).
    pub x_train: Vec<Vec<f64>>,
    /// Dual coefficients (one per support point).
    pub alpha: Vec<f64>,
    /// Kernel bandwidth `γ` for the RBF.
    pub gamma: f64,
    /// Regularisation strength `λ`.
    pub lambda: f64,
}

impl KernelRidgeModel {
    /// RBF kernel value `exp(−γ ‖x − y‖²)`.
    fn rbf(&self, x: &[f64], y: &[f64]) -> f64 {
        let mut acc = 0.0;
        for (xi, yi) in x.iter().zip(y.iter()) {
            let d = xi - yi;
            acc += d * d;
        }
        (-self.gamma * acc).exp()
    }

    /// Predict `y` for a single feature row.
    pub fn predict(&self, x: &[f64]) -> f64 {
        let mut acc = 0.0;
        for (xi, &ai) in self.x_train.iter().zip(self.alpha.iter()) {
            acc += ai * self.rbf(x, xi);
        }
        acc
    }
}

/// Train a kernel ridge regression model with an RBF kernel.
///
/// - `gamma` — kernel bandwidth.
/// - `lambda` — regularisation strength (Tikhonov).
pub fn kernel_ridge_regression(
    x: &[Vec<f64>],
    y: &[f64],
    gamma: f64,
    lambda: f64,
) -> (KernelRidgeModel, TrainingResult) {
    let n = x.len();
    assert_eq!(x.len(), y.len());
    // K_{ij} = exp(−γ ‖x_i − x_j‖²)
    let mut k = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            let mut acc = 0.0;
            for (a, b) in x[i].iter().zip(x[j].iter()) {
                let d = a - b;
                acc += d * d;
            }
            k[i][j] = (-gamma * acc).exp();
        }
    }
    // (K + λ I) α = y
    for i in 0..n {
        k[i][i] += lambda;
    }
    let alpha = solve_dense(&k, y);
    let model = KernelRidgeModel {
        x_train: x.to_vec(),
        alpha: alpha.clone(),
        gamma,
        lambda,
    };
    let mut sse = 0.0;
    let mut sst = 0.0;
    let y_mean = y.iter().sum::<f64>() / n as f64;
    for (row, &yi) in x.iter().zip(y.iter()) {
        let pred = model.predict(row);
        sse += (pred - yi) * (pred - yi);
        sst += (yi - y_mean) * (yi - y_mean);
    }
    let r2 = if sst > 0.0 { 1.0 - sse / sst } else { 0.0 };
    (
        model,
        TrainingResult {
            train_mse: sse / n as f64,
            train_r2: r2,
            num_samples: n,
        },
    )
}

fn solve_dense(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = a.len();
    let mut aug = vec![vec![0.0; n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = a[i][j];
        }
        aug[i][n] = b[i];
    }
    for i in 0..n {
        let mut piv = i;
        for k in (i + 1)..n {
            if aug[k][i].abs() > aug[piv][i].abs() {
                piv = k;
            }
        }
        if aug[piv][i].abs() < 1.0e-15 {
            return vec![0.0; n];
        }
        aug.swap(i, piv);
        let inv = 1.0 / aug[i][i];
        for j in 0..=n {
            aug[i][j] *= inv;
        }
        for k in 0..n {
            if k != i {
                let f = aug[k][i];
                for j in 0..=n {
                    aug[k][j] -= f * aug[i][j];
                }
            }
        }
    }
    let mut out = vec![0.0; n];
    for i in 0..n {
        out[i] = aug[i][n];
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn krr_fits_sine() {
        // y = sin(x) on a 0..2π grid.
        let n = 30;
        let mut x = Vec::with_capacity(n);
        let mut y = Vec::with_capacity(n);
        for i in 0..n {
            let xi = 2.0 * core::f64::consts::PI * (i as f64) / (n as f64 - 1.0);
            x.push(vec![xi]);
            y.push(xi.sin());
        }
        let (model, res) = kernel_ridge_regression(&x, &y, 1.0, 1.0e-3);
        assert!(res.train_r2 > 0.95);
        let p = model.predict(&[core::f64::consts::FRAC_PI_2]);
        assert!((p - 1.0).abs() < 0.1);
    }

    #[test]
    fn krr_two_features() {
        let x = vec![
            vec![0.0, 0.0],
            vec![1.0, 1.0],
            vec![2.0, 2.0],
            vec![3.0, 3.0],
        ];
        let y = vec![0.0, 2.0, 4.0, 6.0];
        let (m, _) = kernel_ridge_regression(&x, &y, 1.0, 1.0e-2);
        let p = m.predict(&[1.5, 1.5]);
        assert!((p - 3.0).abs() < 1.0);
    }
}
