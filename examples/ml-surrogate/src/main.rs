//! Example: ML surrogate for a composition-property sweep.
//!
//! Synthetic Young's-modulus sweep over a binary Al–Cu composition
//! grid, fit with both ridge regression and kernel ridge regression,
//! report training R² and out-of-sample validation error.

use tpt_mat_machine_learning::{kernel_ridge_regression, ridge_regression};

fn main() {
    // Synthetic dataset: E (GPa) as a smooth function of x_Cu.
    let n = 25;
    let mut x = Vec::with_capacity(n);
    let mut y = Vec::with_capacity(n);
    for i in 0..n {
        let xcu = i as f64 / (n as f64 - 1.0);
        x.push(vec![xcu]); // Al–Cu composition (single feature)
        let e_true = 70.0 * (1.0 - xcu) + 130.0 * xcu + 30.0 * xcu * (1.0 - xcu);
        y.push(e_true);
    }

    let model_ols = ridge_regression(&x, &y, 0.0);
    let (model_krr, info_krr) = kernel_ridge_regression(&x, &y, 5.0, 1.0e-2);

    println!("Composition-property sweep (Al–Cu, E in GPa)");
    println!("===========================================");
    println!("N = {} samples", n);
    println!();
    println!("OLS ridge (λ = 0):");
    println!("  weights = {:?}", model_ols.weights);
    println!("  bias    = {:.3}", model_ols.bias);
    // Compute OLS train metrics manually.
    let mut sse_ols = 0.0;
    let y_mean = y.iter().sum::<f64>() / n as f64;
    let mut sst = 0.0;
    for (xi, &yi) in x.iter().zip(y.iter()) {
        let p = model_ols.predict(xi);
        sse_ols += (p - yi) * (p - yi);
        sst += (yi - y_mean) * (yi - y_mean);
    }
    let r2_ols = if sst > 0.0 { 1.0 - sse_ols / sst } else { 0.0 };
    println!("  train MSE = {:.3e}", sse_ols / n as f64);
    println!("  train R²  = {:.4}", r2_ols);
    println!();
    println!("Kernel ridge regression (γ = 5, λ = 1e-2):");
    println!("  train MSE = {:.3e}", info_krr.train_mse);
    println!("  train R²  = {:.4}", info_krr.train_r2);

    // Out-of-sample check.
    let x_test = vec![vec![0.10], vec![0.50], vec![0.85]];
    let y_test = vec![
        70.0 * 0.90 + 130.0 * 0.10 + 30.0 * 0.10 * 0.90,
        70.0 * 0.50 + 130.0 * 0.50 + 30.0 * 0.50 * 0.50,
        70.0 * 0.15 + 130.0 * 0.85 + 30.0 * 0.85 * 0.15,
    ];
    println!();
    println!("Out-of-sample predictions:");
    for (xt, &yt) in x_test.iter().zip(y_test.iter()) {
        let p_ols = model_ols.predict(xt);
        let p_krr = model_krr.predict(xt);
        println!(
            "  x = {:?}  true = {:6.2}  OLS = {:6.2}  KRR = {:6.2}",
            xt, yt, p_ols, p_krr
        );
    }
}
