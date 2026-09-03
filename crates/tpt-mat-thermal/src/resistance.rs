//! Kapitza interfacial thermal-resistance correction.

/// Effective conductivity of a laminate composite with
/// per-interface Kapitza resistance `r_i` (m²·K/W).
///
/// `k_eff = L_total / (Σ L_j / k_j + Σ r_i)`
/// where `L_j` is the thickness of layer `j`, and the
/// Kapitza resistances appear between every pair of adjacent
/// layers (so `r_i` is counted `n-1` times for `n` layers).
pub fn kapitza_correction(k_layers: &[f64], l_layers: &[f64], r_interfaces: &[f64]) -> f64 {
    assert_eq!(k_layers.len(), l_layers.len(), "k and L must match");
    let n = k_layers.len();
    let l_total: f64 = l_layers.iter().sum();
    if l_total <= 0.0 {
        return 0.0;
    }
    let bulk_resistance: f64 = l_layers
        .iter()
        .zip(k_layers.iter())
        .map(|(l, k)| l / k)
        .sum();
    let interface_resistance: f64 = r_interfaces.iter().sum();
    let total = bulk_resistance + interface_resistance;
    if total <= 0.0 {
        return 0.0;
    }
    l_total / total
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn kapitza_correction_zero_interface_recovers_voigt_average() {
        // Two equal-thickness layers, no interface resistance
        // -> arithmetic mean.
        let k = kapitza_correction(&[10.0, 100.0], &[0.5, 0.5], &[]);
        // L_total = 1, resistance = 0.5/10 + 0.5/100 = 0.055
        // k_eff = 1 / 0.055 = 18.18...
        let expected = 1.0 / (0.5 / 10.0 + 0.5 / 100.0);
        assert!(approx(k, expected, 1.0e-9));
    }

    #[test]
    fn kapitza_correction_reduces_effective_conductivity() {
        let k_clean = kapitza_correction(&[10.0, 100.0], &[0.5, 0.5], &[]);
        let k_dirty = kapitza_correction(&[10.0, 100.0], &[0.5, 0.5], &[0.01]);
        assert!(k_dirty < k_clean);
    }
}