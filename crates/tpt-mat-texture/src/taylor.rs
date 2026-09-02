//! Taylor factor: the polycrystal-averaged ratio of macroscopic
//! stress to slip-system shear stress under fully-constrained
//! (Taylor, 1938) plasticity.
//!
//! The classical result for an FCC aggregate with random texture is
//! `M ≈ 3.06` (Taylor's original upper-bound calculation), and for
//! BCC it is `M ≈ 2.75`.
//!
//! **Approximation in this crate.**  A fully-rigorous Taylor factor
//! requires solving a small linear complementarity problem at every
//! strain state (Bishop & Hill, 1951; Kocks, 1970) — the macroscopic
//! stress is the one that activates the 5 independent slip systems
//! that minimise `Σ |γ^α|`.  This crate ships a simple upper-bound
//! proxy that uses the per-grain maximum Schmid factor:
//!
//! `M ≈ ⟨ 1 / max_schmid_factor ⟩`
//!
//! which converges to `M ≈ 2.0–2.3` for FCC (rather than 3.06)
//! because single-system activation is over-generous.  A full
//! `TaylorFactorSolver` (LCP, Bishop-Hill) lands in Phase 5 alongside
//! the RVE / homogenization stack.

use tpt_mat_crystallography::{CrystalStructure, SlipSystem};

/// Compute the Taylor factor for an FCC polycrystal with the supplied
/// `n_samples` random tensile axes.  Returns the average
/// `M = ⟨ 1 / max_schmid ⟩` over the sample directions.
///
/// Pass `seed` for deterministic results.
pub fn taylor_factor(n_samples: usize, seed: u64) -> f64 {
    let slips = CrystalStructure::FCC.slip_systems();
    taylor_factor_with_slips(&slips, n_samples, seed)
}

/// Taylor factor with an explicit slip-system list.
pub fn taylor_factor_with_slips(
    slips: &[SlipSystem],
    n_samples: usize,
    seed: u64,
) -> f64 {
    let mut rng = Lcg::new(seed.max(1));
    let mut acc = 0.0;
    let mut count = 0;
    for _ in 0..n_samples {
        let axis = random_unit_vector(&mut rng);
        let m = SlipSystem::max_schmid_factor(slips, axis);
        if m > 0.0 {
            acc += 1.0 / m;
            count += 1;
        }
    }
    if count == 0 {
        0.0
    } else {
        acc / count as f64
    }
}

fn random_unit_vector(rng: &mut Lcg) -> [f64; 3] {
    // Marsaglia (1972) method.
    loop {
        let u = 2.0 * rng.next_f64() - 1.0;
        let v = 2.0 * rng.next_f64() - 1.0;
        let s = u * u + v * v;
        if s < 1.0 && s > 0.0 {
            let factor = 2.0 * (1.0 - s).sqrt();
            return [u * factor, v * factor, 1.0 - 2.0 * s];
        }
    }
}

struct Lcg(u64);
impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed)
    }
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0
    }
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn taylor_factor_proxy_for_random_fcc_is_about_2() {
        // Single-Schmid proxy converges to ~2.0–2.3 for FCC.  The
        // full Taylor factor of 3.06 requires the Bishop–Hill LCP
        // solver (Phase 5).
        let m = taylor_factor(4096, 0xDEADBEEF);
        assert!((m - 2.15).abs() < 0.30, "M = {m}");
    }

    #[test]
    fn taylor_factor_is_deterministic() {
        let a = taylor_factor(256, 42);
        let b = taylor_factor(256, 42);
        assert_eq!(a, b);
    }
}