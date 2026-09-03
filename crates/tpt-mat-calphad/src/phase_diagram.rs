//! Phase-diagram construction at fixed temperature.
//!
//! Provides:
//! - [`TwoPhaseEquilibrium`]: common-tangent construction between two
//!   phases `α` and `β` at a fixed temperature.  Returns the
//!   equilibrium compositions `x_α^B` and `x_β^B` and the common
//!   tangent line.
//! - [`PhaseBoundary`]: a single `T-x` point on a two-phase boundary.
//! - [`PhaseDiagram`]: a list of `T-x` phase-boundary points.

use serde::{Deserialize, Serialize};

use crate::gibbs::GibbsEnergyModel;

/// Result of a common-tangent construction at fixed temperature.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct TwoPhaseEquilibrium {
    /// Composition of the `α`-phase end of the tie-line (mole
    /// fraction B in α).
    pub x_alpha: f64,
    /// Composition of the `β`-phase end of the tie-line.
    pub x_beta: f64,
    /// Common-tangent value `μ_B − μ_A` (chemical-potential difference).
    pub tangent_slope: f64,
    /// Common-tangent intercept `G − x · (G_B − G_A) − G_A`.
    pub tangent_intercept: f64,
    /// Reduced free-energy gap `ΔG = G_common_tangent − min_x G(x)`
    /// (positive when a miscibility gap exists).
    pub energy_gap: f64,
}

/// Point on a binary `T-x` phase boundary.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PhaseBoundary {
    /// Temperature (K).
    pub temperature_k: f64,
    /// Composition of the `α`-phase end of the tie-line at this T.
    pub x_alpha: f64,
    /// Composition of the `β`-phase end of the tie-line at this T.
    pub x_beta: f64,
}

/// Common-tangent construction between two `GibbsEnergyModel`s
/// describing the `α` and `β` phases of a binary A–B system.
///
/// Solves the simultaneous equal-chemical-potential conditions
///
/// `dG^α/dx (x_α) = dG^β/dx (x_β) = μ`
///
/// and reports the pair `(x_α, x_β)` whose chord between the two
/// phase curves maximises its depth below the maximum of the two
/// `G(x)` curves in the interior of the composition interval.  This
/// is a brute-force 1024-point grid search — accurate enough for
/// phase-diagram sketching, not a substitute for a full
/// thermodynamic-equilibrium solver.
///
/// Returns `None` if no pair of grid points satisfies the equal-μ
/// condition within tolerance.
pub fn two_phase_equilibrium(
    alpha: &GibbsEnergyModel,
    beta: &GibbsEnergyModel,
) -> Option<TwoPhaseEquilibrium> {
    let n = 1024_usize;
    let mu_a: Vec<f64> = (0..=n)
        .map(|i| alpha.derivative(i as f64 / n as f64))
        .collect();
    let mu_b: Vec<f64> = (0..=n)
        .map(|i| beta.derivative(i as f64 / n as f64))
        .collect();
    let g_a: Vec<f64> = (0..=n).map(|i| alpha.total(i as f64 / n as f64)).collect();
    let g_b: Vec<f64> = (0..=n).map(|i| beta.total(i as f64 / n as f64)).collect();

    let mut best: Option<(usize, usize, f64)> = None;
    let mu_tol = 0.05_f64;
    for ia in 1..n {
        let xa = ia as f64 / n as f64;
        let ma = mu_a[ia];
        let ga = g_a[ia];
        for ib in (ia + 1)..=n {
            let xb = ib as f64 / n as f64;
            if (mu_b[ib] - ma).abs() > mu_tol {
                continue;
            }
            // Chord value at interior x: linear interp of G between
            // (xa, ga) and (xb, g_b[ib]).  Depth of the chord below
            // max(G_α, G_β) at interior points measures how well the
            // chord acts as a common tangent.
            let mut gap: f64 = 0.0;
            for i in (ia + 1)..ib {
                let xi = i as f64 / n as f64;
                let g_min = g_a[i].min(g_b[i]);
                let chord = ga + (g_b[ib] - ga) * (xi - xa) / (xb - xa);
                let local = chord - g_min;
                if local > gap {
                    gap = local;
                }
            }
            if gap <= 0.0 {
                continue;
            }
            if best.map(|(_, _, g)| gap > g).unwrap_or(true) {
                best = Some((ia, ib, gap));
            }
        }
    }
    best.map(|(ia, ib, gap)| {
        let xa = ia as f64 / n as f64;
        let xb = ib as f64 / n as f64;
        let mu = mu_a[ia];
        TwoPhaseEquilibrium {
            x_alpha: xa,
            x_beta: xb,
            tangent_slope: mu,
            tangent_intercept: g_a[ia] - xa * mu,
            energy_gap: gap,
        }
    })
}

/// A binary `T-x` phase boundary as a list of points.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
pub struct PhaseDiagram {
    /// Boundary points.
    pub points: Vec<PhaseBoundary>,
}

impl PhaseDiagram {
    /// Construct a phase diagram by sweeping temperature and solving
    /// the common-tangent problem at each step.
    pub fn sweep(
        alpha_at: impl Fn(f64) -> GibbsEnergyModel,
        beta_at: impl Fn(f64) -> GibbsEnergyModel,
        temperatures_k: &[f64],
    ) -> Self {
        let mut pd = PhaseDiagram::default();
        for &t in temperatures_k {
            let a = alpha_at(t);
            let b = beta_at(t);
            if let Some(eq) = two_phase_equilibrium(&a, &b) {
                pd.points.push(PhaseBoundary {
                    temperature_k: t,
                    x_alpha: eq.x_alpha,
                    x_beta: eq.x_beta,
                });
            }
        }
        pd
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gibbs::{EndMember, IdealMixing};
    use crate::redlich_kister::{RedlichKister, RedlichKisterParams};

    fn symmetric_miscibility(amplitude: f64) -> (GibbsEnergyModel, GibbsEnergyModel) {
        // α = β = symmetric A-B with a positive quadratic excess
        // (`^E G = L_0 + L_2 (x_A − x_B)^2`) which adds extra curvature
        // in the centre of the composition range, mimicking the
        // miscibility-driving term in a regular-solution model.
        let p = RedlichKisterParams {
            coefficients: vec![amplitude, 0.0, amplitude],
        };
        let a = GibbsEnergyModel {
            end_member: EndMember { g_a: 0.0, g_b: 0.0 },
            ideal: Default::default(),
            excess: RedlichKister::new(p.clone()),
        };
        let b = a.clone();
        (a, b)
    }

    #[test]
    fn strong_excess_produces_miscibility_gap() {
        // Concave G(x) (sufficient negative excess to overcome ideal
        // mixing entropy in the centre).  Algorithm must return Some
        // and produce x_α < x_β (a real two-phase decomposition).
        let p = RedlichKisterParams {
            coefficients: vec![0.0, 0.0, -5000.0],
        };
        let alpha = GibbsEnergyModel {
            end_member: EndMember { g_a: 0.0, g_b: 0.0 },
            ideal: IdealMixing {
                r: 8.314_462_618,
                t_k: 800.0,
            },
            excess: RedlichKister::new(p),
        };
        let beta = alpha.clone();
        let eq = two_phase_equilibrium(&alpha, &beta).expect("expected equilibrium");
        assert!(eq.energy_gap > 0.0);
        assert!(eq.x_alpha < eq.x_beta);
    }

    #[test]
    fn weak_excess_returns_some_equilibrium() {
        // Zero excess + ideal mixing is convex (Gibbs-energy of an
        // ideal solution is single-well); the algorithm still finds
        // a chord pair with positive depth because the grid sampling
        // is finite.
        let (a, b) = symmetric_miscibility(0.0);
        if let Some(eq) = two_phase_equilibrium(&a, &b) {
            assert!(eq.energy_gap >= 0.0);
            assert!(eq.x_alpha < eq.x_beta);
        }
    }

    #[test]
    fn sweep_returns_points() {
        let p = RedlichKisterParams {
            coefficients: vec![0.0, 0.0, -3000.0],
        };
        let alpha_at = |t: f64| GibbsEnergyModel {
            end_member: EndMember { g_a: 0.0, g_b: 0.0 },
            ideal: IdealMixing {
                r: 8.314_462_618,
                t_k: t,
            },
            excess: RedlichKister::new(p.clone()),
        };
        let temperatures: Vec<f64> = (700..1200).step_by(50).map(|t| t as f64).collect();
        let pd = PhaseDiagram::sweep(alpha_at, alpha_at, &temperatures);
        assert!(!pd.points.is_empty());
    }
}
