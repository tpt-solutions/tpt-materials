//! Cross-crate verification tests called out in `todo.md` Phase 1.
//!
//! These are deliberately written to exercise the public API of
//! multiple crates together; per-crate unit tests live in the
//! individual crates.

use tpt_mat_crystallography::{CrystalStructure, SlipSystem};
use tpt_math_linalg_fixed::{SymMat3, Vec3};

#[test]
fn fcc_has_exactly_12_slip_systems() {
    assert_eq!(CrystalStructure::FCC.slip_systems().len(), 12);
}

#[test]
fn bcc_has_exactly_12_slip_systems() {
    assert_eq!(CrystalStructure::BCC.slip_systems().len(), 12);
}

#[test]
fn hcp_has_full_slip_inventory() {
    let s = CrystalStructure::HCP.slip_systems();
    assert_eq!(s.len(), 24, "HCP must have 3+3+6+12 = 24 systems");
}

#[test]
fn schmid_tensor_is_symmetric_for_random_unit_vectors() {
    // Property-based: for any unit s and any unit n, P = sym(s ⊗ n) is
    // a symmetric matrix.  We sample a handful of representative
    // directions.
    let cases: [(Vec3, Vec3); 6] = [
        (Vec3::X, Vec3::Z),
        (Vec3::Y, Vec3::X),
        (Vec3::Z, Vec3::Y),
        (
            Vec3::new(1.0, 1.0, 0.0).normalized(),
            Vec3::new(1.0, -1.0, 1.0).normalized(),
        ),
        (
            Vec3::new(2.0, -1.0, 1.0).normalized(),
            Vec3::new(-1.0, 0.5, 1.5).normalized(),
        ),
        (
            Vec3::new(-3.0, 0.7, 2.0).normalized(),
            Vec3::new(0.1, -0.4, 1.0).normalized(),
        ),
    ];
    for (s, n) in cases {
        let p: SymMat3 = SlipSystem::schmid_tensor(s, n);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (p.get(i, j) - p.get(j, i)).abs() < 1e-12,
                    "Schmid tensor not symmetric at ({i}, {j}) for s={s:?}, n={n:?}"
                );
            }
        }
    }
}

#[test]
fn schmid_factor_max_for_uniaxial_tension_is_cos_phi_cos_lambda() {
    // Single slip system: s = (1,0,0), n = (0,1,0).  Uniaxial tension
    // along d = (1,1,1)/√3.  Schmid factor should be
    // |cos(angle(s, d)) · cos(angle(n, d))| = |1/√3 · 1/√3| = 1/3.
    let s = Vec3::X;
    let n = Vec3::Y;
    let d: [f64; 3] = [1.0, 1.0, 1.0];
    let slips =
        vec![
            SlipSystem::from_vectors(s, n, 1.0, tpt_mat_crystallography::SlipFamily::Fcc110)
                .unwrap(),
        ];
    let m = SlipSystem::max_schmid_factor(&slips, d);
    assert!((m - 1.0 / 3.0).abs() < 1e-9, "got {m}");
}
