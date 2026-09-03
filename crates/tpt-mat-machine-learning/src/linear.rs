//! Ordinary / ridge linear regression (closed-form).

use serde::{Deserialize, Serialize};

use super::TrainingResult;

/// Fitted linear model `y = X · w + b`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LinearModel {
    /// Per-feature weights.
    pub weights: Vec<f64>,
    /// Bias (intercept).
    pub bias: f64,
}

impl LinearModel {
    /// Predict `y` for a single feature row.
    pub fn predict(&self, x: &[f64]) -> f64 {
        debug_assert_eq!(x.len(), self.weights.len());
        let mut acc = self.bias;
        for (w, xi) in self.weights.iter().zip(x.iter()) {
            acc += w * xi;
        }
        acc
    }
}

/// Train an ordinary least-squares model.
///
/// `x` — design matrix (`n_samples × n_features`), row-major.
/// `y` — targets (`n_samples`).
///
/// Returns the fitted [`LinearModel`].
pub fn linear_regression(x: &[Vec<f64>], y: &[f64]) -> LinearModel {
    ridge_regression(x, y, 0.0)
}

/// Train a ridge regression model with regularisation strength `λ`.
///
/// Solves `(Xᵀ X + λ I) w = Xᵀ (y − b)` with the bias fitted
/// separately to keep the penalty off the intercept term.
pub fn ridge_regression(x: &[Vec<f64>], y: &[f64], lambda: f64) -> LinearModel {
    assert!(!x.is_empty());
    assert_eq!(x.len(), y.len());
    let n = x.len();
    let d = x[0].len();
    // Mean-impute the targets for the bias.
    let mut y_mean = 0.0;
    for &yi in y {
        y_mean += yi;
    }
    y_mean /= n as f64;
    let mut x_mean = vec![0.0; d];
    for row in x {
        for j in 0..d {
            x_mean[j] += row[j];
        }
    }
    for j in 0..d {
        x_mean[j] /= n as f64;
    }
    let mut xtx = vec![vec![0.0; d]; d];
    let mut xty = vec![0.0; d];
    for (row, &yi) in x.iter().zip(y.iter()) {
        let xc: Vec<f64> = (0..d).map(|j| row[j] - x_mean[j]).collect();
        let yc = yi - y_mean;
        for i in 0..d {
            xty[i] += xc[i] * yc;
            for j in 0..d {
                xtx[i][j] += xc[i] * xc[j];
            }
        }
    }
    for i in 0..d {
        xtx[i][i] += lambda;
    }
    let w = solve_symmetric(&xtx, &xty).unwrap_or_else(|| vec![0.0; d]);
    let bias = y_mean - (0..d).map(|j| w[j] * x_mean[j]).sum::<f64>();
    LinearModel { weights: w, bias }
}

/// Train and report the training MSE / R².
pub fn train(x: &[Vec<f64>], y: &[f64], lambda: f64) -> (LinearModel, TrainingResult) {
    let model = ridge_regression(x, y, lambda);
    let n = x.len();
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

/// Solve `A w = b` for a symmetric positive-definite `A` via
/// Cholesky decomposition.  Returns `None` if the decomposition
/// fails.
fn solve_symmetric(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = a.len();
    // Cholesky: A = L Lᵀ
    let mut l = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i][j];
            for k in 0..j {
                sum -= l[i][k] * l[j][k];
            }
            if i == j {
                if sum <= 0.0 {
                    return None;
                }
                l[i][j] = sum.sqrt();
            } else {
                l[i][j] = sum / l[j][j];
            }
        }
    }
    // Solve L z = b
    let mut z = vec![0.0; n];
    for i in 0..n {
        let mut sum = b[i];
        for j in 0..i {
            sum -= l[i][j] * z[j];
        }
        z[i] = sum / l[i][i];
    }
    // Solve Lᵀ w = z
    let mut w = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = z[i];
        for j in (i + 1)..n {
            sum -= l[j][i] * w[j];
        }
        w[i] = sum / l[i][i];
    }
    Some(w)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn linear_fit_perfect_line() {
        let x = vec![vec![1.0], vec![2.0], vec![3.0], vec![4.0]];
        let y = vec![2.0, 4.0, 6.0, 8.0];
        let (model, res) = train(&x, &y, 0.0);
        assert!(approx(res.train_mse, 0.0, 1.0e-9));
        assert!(approx(model.weights[0], 2.0, 1.0e-6));
        assert!(approx(model.bias, 0.0, 1.0e-6));
    }

    #[test]
    fn linear_fit_noisy_data_reasonable_r2() {
        let x = vec![vec![0.0], vec![1.0], vec![2.0], vec![3.0], vec![4.0]];
        let y = vec![0.1, 1.1, 1.9, 3.1, 3.9];
        let (_, res) = train(&x, &y, 0.0);
        assert!(res.train_r2 > 0.99);
    }

    #[test]
    fn ridge_reduces_overfitting() {
        let x = vec![vec![0.0, 1.0], vec![1.0, 0.0], vec![1.0, 1.0]];
        let y = vec![1.0, 1.0, 0.5];
        let (m_ols, _) = train(&x, &y, 0.0);
        let (m_ridge, _) = train(&x, &y, 1.0);
        // Ridge weights should be smaller in magnitude than OLS.
        let ols_norm: f64 = m_ols.weights.iter().map(|w| w * w).sum();
        let ridge_norm: f64 = m_ridge.weights.iter().map(|w| w * w).sum();
        assert!(ridge_norm < ols_norm);
    }

    #[test]
    fn predict_multi_feature() {
        let x = vec![vec![1.0, 2.0], vec![3.0, 4.0], vec![5.0, 6.0]];
        let y = vec![1.0, 2.0, 3.0];
        let (m, _) = train(&x, &y, 0.0);
        let p = m.predict(&[7.0, 8.0]);
        assert!(p.is_finite());
    }
}
