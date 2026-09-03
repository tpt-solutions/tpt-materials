//! Polynomial feature augmentation.

/// Append polynomial features (degree `d`) to a 1-D design matrix.
///
/// Produces monomials with non-decreasing index sequences so that
/// `poly_features([2, 3], 2) = [1, 2, 3, 4, 6, 9]`.
pub fn polynomial_features(x: &[f64], degree: usize) -> Vec<f64> {
    let d = x.len();
    let mut out = Vec::new();
    out.push(1.0);
    // Monomials of increasing total degree.
    let mut indices: Vec<usize> = Vec::new();
    fn emit(out: &mut Vec<f64>, x: &[f64], indices: &[usize]) {
        let mut v = 1.0;
        for &i in indices {
            v *= x[i];
        }
        out.push(v);
    }
    fn recurse(
        out: &mut Vec<f64>,
        x: &[f64],
        indices: &mut Vec<usize>,
        start: usize,
        remaining: usize,
    ) {
        if remaining == 0 {
            emit(out, x, indices);
            return;
        }
        for i in start..x.len() {
            indices.push(i);
            recurse(out, x, indices, i, remaining - 1);
            indices.pop();
        }
    }
    for deg in 1..=degree {
        recurse(&mut out, x, &mut indices, 0, deg);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn poly_features_constant() {
        let f = polynomial_features(&[1.0, 2.0], 1);
        assert_eq!(f, vec![1.0, 1.0, 2.0]);
    }

    #[test]
    fn poly_features_degree_two() {
        let f = polynomial_features(&[2.0, 3.0], 2);
        // [1, x1, x2, x1^2, x1*x2, x2^2] = [1, 2, 3, 4, 6, 9]
        assert_eq!(f, vec![1.0, 2.0, 3.0, 4.0, 6.0, 9.0]);
    }

    #[test]
    fn poly_features_degree_three() {
        let f = polynomial_features(&[2.0], 3);
        // [1, x, x², x³] = [1, 2, 4, 8]
        assert_eq!(f, vec![1.0, 2.0, 4.0, 8.0]);
    }
}
