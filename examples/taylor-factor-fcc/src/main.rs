//! Example: Bishop–Hill (1951) Taylor factor for an FCC polycrystal
//! under uniaxial tension.
//!
//! Computes the Taylor factor `M = σ_y / τ_c` for several specific
//! crystallographic tensile axes, then averages over a random sample
//! to estimate the random-texture FCC Taylor factor
//! (T = 3.06, Taylor 1938).
//!
//! Run with: `cargo run --release -p taylor-factor-fcc`.

use tpt_mat_crystallography::CrystalStructure;
use tpt_mat_rve::bishop_hill_taylor_factor_axis;

fn main() {
    let slips = CrystalStructure::FCC.slip_systems();
    let tau_c = 1.0_f64;
    println!("FCC Taylor factor (CRSS τ_c = {tau_c}) along high-symmetry tensile axes:");
    println!("  axis             M          n_active");
    println!("  ------------     ------     --------");
    for (label, dir) in [
        ("[001]", [0.0_f64, 0.0, 1.0]),
        ("[011]", [0.0, 1.0, 1.0]),
        ("[111]", [1.0, 1.0, 1.0]),
        ("[012]", [0.0, 1.0, 2.0]),
        ("[112]", [1.0, 1.0, 2.0]),
        ("[123]", [1.0, 2.0, 3.0]),
    ] {
        let n = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
        let axis = [dir[0] / n, dir[1] / n, dir[2] / n];
        let r = bishop_hill_taylor_factor_axis(axis, tau_c, Some(&slips));
        println!(
            "  {:>10}     {:>6.3}     {:>8}",
            label, r.taylor_factor, r.n_active
        );
    }

    // Random sample: estimate the random-texture Taylor factor.
    let mut rng = 0xDEADBEEFu64;
    let n_samples = 256;
    let mut m_acc = 0.0_f64;
    for _ in 0..n_samples {
        let axis = loop {
            let u = 2.0 * (rng as f64 / u64::MAX as f64) - 1.0;
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            let v = 2.0 * (rng as f64 / u64::MAX as f64) - 1.0;
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            let s = u * u + v * v;
            if s < 1.0 && s > 0.0 {
                let factor = 2.0 * (1.0 - s).sqrt();
                break [u * factor, v * factor, 1.0 - 2.0 * s];
            }
        };
        let r = bishop_hill_taylor_factor_axis(axis, tau_c, Some(&slips));
        m_acc += r.taylor_factor;
    }
    let m_avg = m_acc / n_samples as f64;
    println!();
    println!("Average over {n_samples} random directions: M = {m_avg:.3}");
    println!(
        "(Classical Taylor 1938: M ≈ 3.06; our L2 pseudo-inverse underestimates the L1 minimum.)"
    );
}
