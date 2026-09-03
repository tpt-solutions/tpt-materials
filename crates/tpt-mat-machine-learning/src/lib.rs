//! Materials-informatics machine-learning surrogates.
//!
//! A small collection of analytic regression / classification models
//! suitable for surrogate-modelling a CP-FEM sweep, a homogenization
//! map, or a CALPHAD parameter fit.
//!
//! # Models
//!
//! - [`LinearRegression`] — ordinary least squares with optional
//!   ridge (Tikhonov) regularisation.
//! - [`KernelRidgeRegression`] — kernel ridge regression with
//!   RBF or polynomial kernels.
//! - [`PolynomialFeatures`] — augment a design matrix with
//!   cross-products up to a chosen degree.
//! - [`PhasePredictor`] — multi-label softmax classification for
//!   phase-fident predictions.
//!
//! The models are deliberately simple (closed-form, no GPU) so they
//! can serve as the **baseline** that any deep-learning surrogate in
//! the wider `tpt-mat-*` ecosystem is benchmarked against.

#![warn(missing_docs)]

mod features;
mod linear;
mod phase;
mod rbf;

use serde::{Deserialize, Serialize};

pub use features::polynomial_features;
pub use linear::{linear_regression, ridge_regression, LinearModel};
pub use phase::{PhasePrediction, PhasePredictor};
pub use rbf::{kernel_ridge_regression, KernelRidgeModel};

/// Composition descriptor (in atomic fraction / wt pct).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Composition {
    /// Element symbol or label.
    pub components: Vec<(String, f64)>,
}

impl Composition {
    /// Convert to a feature vector (just the component fractions in
    /// the supplied order).
    pub fn to_features(&self) -> Vec<f64> {
        self.components.iter().map(|(_, x)| *x).collect()
    }
}

/// A machine-learning model category.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum MlModel {
    /// Regression model (single scalar output).
    PropertyPredictor,
    /// Multi-class classifier for phase / category prediction.
    PhasePredictor,
    /// General surrogate for a CP-FEM / homogenisation sweep.
    SurrogateModel,
}

/// Result of a model training run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TrainingResult {
    /// Mean-squared error on the training set.
    pub train_mse: f64,
    /// R² score on the training set.
    pub train_r2: f64,
    /// Number of training samples.
    pub num_samples: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn composition_to_features() {
        let c = Composition {
            components: vec![("Al".to_string(), 0.95), ("Cu".to_string(), 0.05)],
        };
        let f = c.to_features();
        assert_eq!(f, vec![0.95, 0.05]);
    }
}
