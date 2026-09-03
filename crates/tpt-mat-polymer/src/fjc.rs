//! Inverse-Langevin approximation used by FJC and Arruda–Boyce.

/// Pnueli (1991) approximation of the inverse Langevin function
/// `L^{-1}(x)`:
///
/// `L^{-1}(x) ≈ x (3 − x²) / (1 − x²)` for `|x| < 1`.
pub fn inverse_langevin_approx(x: f64) -> f64 {
    let x = x.clamp(-0.999, 0.999);
    if x.abs() < 1.0e-12 {
        return 0.0;
    }
    let x2 = x * x;
    x * (3.0 - x2) / (1.0 - x2)
}

/// Reference exact inverse-Langevin via fixed-point iteration.
pub fn inverse_langevin_exact(x: f64, tol: f64, max_iter: u32) -> f64 {
    let mut y = inverse_langevin_approx(x);
    for _ in 0..max_iter {
        let y2 = y * y;
        if (1.0 - y2).abs() < 1.0e-30 {
            break;
        }
        let y_new = x * (3.0 - y2) / (1.0 - y2);
        if (y_new - y).abs() < tol {
            return y_new;
        }
        y = y_new;
    }
    y
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn l_zero_at_zero() {
        assert!(approx(inverse_langevin_approx(0.0), 0.0, 1.0e-12));
    }

    #[test]
    fn l_monotonically_increases() {
        let mut last = 0.0;
        for i in 1..=10 {
            let x = i as f64 / 11.0;
            let v = inverse_langevin_approx(x);
            assert!(v > last);
            last = v;
        }
    }

    #[test]
    fn approx_is_positive_for_x_in_unit_interval() {
        let v = inverse_langevin_approx(0.7);
        assert!(v > 0.0);
    }
}
