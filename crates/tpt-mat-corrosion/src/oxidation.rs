//! High-temperature oxidation (Wagner parabolic scale growth).
//!
//! Models for the growth of a protective oxide scale on a metal at
//! elevated temperature.  In the parabolic regime (Wagner 1933) the
//! scale thickness `x` obeys
//!
//! `x² = 2 k_p t + x_0²`
//!
//! where `k_p` is the parabolic rate constant (typically Arrhenius
//! in temperature).  The rate constant has units of m²/s and is
//! related to the more familiar mass-gain parabolic constant `k_p_m`
//! by
//!
//! `k_p = k_p_m² / (ρ_scale²) · (M_O / M_scale)² · (z_e / z_e_O)²`
//!
//! where `M_O`/`M_scale` are the molar masses of oxygen and the
//! oxide, and `z_e`/`z_e_O` are the electron counts in the
//! metal-dissolution and oxygen-reduction half-reactions.
//!
//! Wagner–Hauffe semiconductor doping rules (p-type oxide growth
//! accelerated by acceptor cations, n-type accelerated by donor
//! cations) are encoded qualitatively through a `doping_effect`
//! factor that scales `k_p`.
//!
//! Breakaway oxidation (stringer 1965; payer et al. 1980) — the
//! transition from protective parabolic growth to linear or
//! catastrophic growth when the scale cracks or spalls — is exposed
//! as a [`BreakawayCriterion`] predicate over local stress/strain or
//! scale-thickness thresholds.

use serde::{Deserialize, Serialize};

/// Wagner parabolic rate constant `k_p` (m²/s) at temperature `T`
/// (K) and its Arrhenius parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ParabolicRateConstant {
    /// Pre-exponential factor `k_p^0` (m²/s).
    pub k_p_0: f64,
    /// Activation energy `Q` (J/mol).
    pub activation_energy: f64,
}

impl ParabolicRateConstant {
    /// Evaluate `k_p(T) = k_p_0 * exp(-Q / (R T))`.
    pub fn at(&self, temperature_k: f64, gas_constant: f64) -> f64 {
        self.k_p_0 * (-self.activation_energy / (gas_constant * temperature_k)).exp()
    }
}

/// Parabolic scale-growth state and predictive driver.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ScaleGrowth {
    /// Parabolic rate constant at the operating temperature (m²/s).
    pub k_p: f64,
    /// Initial scale thickness `x_0` (m).  Defaults to zero for a
    /// freshly exposed surface.
    pub initial_thickness: f64,
}

impl ScaleGrowth {
    /// Construct from explicit `k_p`.
    pub fn new(k_p: f64) -> Self {
        Self {
            k_p,
            initial_thickness: 0.0,
        }
    }

    /// Construct from a [`ParabolicRateConstant`] evaluated at `T`.
    pub fn at_temperature(
        params: ParabolicRateConstant,
        temperature_k: f64,
        gas_constant: f64,
    ) -> Self {
        Self {
            k_p: params.at(temperature_k, gas_constant),
            initial_thickness: 0.0,
        }
    }

    /// Thickness `x(t)` from `x² = 2 k_p t + x_0²`.
    pub fn thickness(&self, time_s: f64) -> f64 {
        (2.0 * self.k_p * time_s + self.initial_thickness.powi(2)).sqrt()
    }

    /// Mass gain per unit area `Δm/A` (kg/m²).  Uses
    /// `Δm/A = ρ_scale · x · (M_metal / (z_e · M_O))`.  Pass the
    /// oxide density (kg/m³) and the metal/oxygen molar masses
    /// (kg/mol) plus the metal valence.
    pub fn mass_gain(
        &self,
        time_s: f64,
        oxide_density_kg_per_m3: f64,
        molar_mass_metal_kg_per_mol: f64,
        molar_mass_oxygen_kg_per_mol: f64,
        metal_valence: f64,
    ) -> f64 {
        let x = self.thickness(time_s);
        let n_o = oxide_density_kg_per_m3 * x / molar_mass_oxygen_kg_per_mol;
        n_o * molar_mass_metal_kg_per_mol / metal_valence
    }

    /// Instantaneous scale-growth rate `dx/dt = k_p / x`.
    pub fn rate(&self, time_s: f64) -> f64 {
        let x = self.thickness(time_s);
        if x < 1.0e-30 {
            0.0
        } else {
            self.k_p / x
        }
    }
}

/// Wagner–Hauffe doping rule.  `factor` is a multiplier applied to
/// `k_p`: positive values accelerate oxidation (p-type oxide +
/// acceptor cations OR n-type oxide + donor cations), negative
/// values decelerate (the opposite pair).  Typical magnitudes are
/// `|factor| ≤ 3`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DopingEffect {
    /// Multiplicative correction to `k_p`.
    pub factor: f64,
}

impl DopingEffect {
    /// Apply to a base [`ParabolicRateConstant`] to produce a
    /// corrected one.
    pub fn apply(&self, base: ParabolicRateConstant) -> ParabolicRateConstant {
        ParabolicRateConstant {
            k_p_0: base.k_p_0 * self.factor,
            activation_energy: base.activation_energy,
        }
    }
}

/// Breakaway criterion: when the scale exceeds `critical_thickness`
/// OR the local strain exceeds `critical_strain`, parabolic growth
/// gives way to linear / catastrophic growth.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BreakawayCriterion {
    /// Critical thickness `x_c` (m).
    pub critical_thickness: f64,
    /// Critical in-plane strain (dimensionless).
    pub critical_strain: f64,
}

impl BreakawayCriterion {
    /// `true` if either threshold is exceeded at the given state.
    pub fn triggered(&self, thickness: f64, in_plane_strain: f64) -> bool {
        thickness >= self.critical_thickness || in_plane_strain.abs() >= self.critical_strain
    }
}

/// Linear / breakaway growth rate `dx/dt = k_l` (m/s) once the
/// protective scale has failed.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LinearBreakawayRate {
    /// Linear rate constant (m/s).
    pub k_l: f64,
}

impl LinearBreakawayRate {
    /// Construct.
    pub fn new(k_l: f64) -> Self {
        Self { k_l }
    }

    /// Thickness `x(t) = x_breakaway + k_l (t - t_breakaway)`.
    pub fn thickness_from_breakaway(
        &self,
        time_s: f64,
        time_breakaway_s: f64,
        thickness_breakaway: f64,
    ) -> f64 {
        thickness_breakaway + self.k_l * (time_s - time_breakaway_s).max(0.0)
    }
}

/// Pick the more conservative (thicker) of parabolic and linear
/// breakaway growth models — the classical Evans picture where the
/// oxide first grows parabolically, then transitions to linear
/// once the scale fails.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct OxidationModel {
    /// Parabolic growth parameters.
    pub parabolic: ScaleGrowth,
    /// Breakaway criterion.
    pub breakaway: BreakawayCriterion,
    /// Linear rate once breakaway triggers.
    pub linear: LinearBreakawayRate,
}

impl OxidationModel {
    /// `x(t)` with breakaway transition.
    pub fn thickness(&self, time_s: f64, in_plane_strain: f64) -> f64 {
        let x_par = self.parabolic.thickness(time_s);
        if !self.breakaway.triggered(x_par, in_plane_strain) {
            return x_par;
        }
        // Find the approximate breakaway time by bisection.
        let mut lo = 0.0_f64;
        let mut hi = time_s.max(1.0e-3);
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            let x_mid = self.parabolic.thickness(mid);
            if self.breakaway.triggered(x_mid, in_plane_strain) {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let t_b = 0.5 * (lo + hi);
        let x_b = self.parabolic.thickness(t_b);
        self.linear.thickness_from_breakaway(time_s, t_b, x_b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn parabolic_thickness_follows_x_squared_law() {
        let g = ScaleGrowth::new(1.0e-12);
        let x1 = g.thickness(1.0e3);
        let x4 = g.thickness(4.0e3);
        assert!(approx(x4 / x1, 2.0, 1.0e-9));
    }

    #[test]
    fn parabolic_rate_decreases_with_time() {
        let g = ScaleGrowth::new(1.0e-12);
        let r1 = g.rate(100.0);
        let r2 = g.rate(10_000.0);
        assert!(r1 > r2);
    }

    #[test]
    fn arrhenius_kp_doubles_over_reasonable_temperature_range() {
        let p = ParabolicRateConstant {
            k_p_0: 1.0e-6,
            activation_energy: 150.0e3,
        };
        let r = 8.314_462_618;
        let k800 = p.at(800.0, r);
        let k900 = p.at(900.0, r);
        // Activation energy of 150 kJ/mol gives ~25x increase per 100K
        assert!(k900 > k800 * 10.0);
    }

    #[test]
    fn mass_gain_grows_monotonically() {
        let g = ScaleGrowth::new(1.0e-12);
        let m1 = g.mass_gain(100.0, 5000.0, 0.055_845, 0.016, 2.0);
        let m2 = g.mass_gain(10_000.0, 5000.0, 0.055_845, 0.016, 2.0);
        assert!(m2 > m1);
    }

    #[test]
    fn doping_effect_scales_kp_by_factor() {
        let p = ParabolicRateConstant {
            k_p_0: 1.0e-6,
            activation_energy: 100.0e3,
        };
        let accel = DopingEffect { factor: 2.0 };
        let p2 = accel.apply(p);
        assert!(approx(p2.k_p_0, 2.0e-6, 1.0e-15));
    }

    #[test]
    fn breakaway_triggers_on_thickness_or_strain() {
        let c = BreakawayCriterion {
            critical_thickness: 1.0e-6,
            critical_strain: 0.01,
        };
        assert!(c.triggered(1.5e-6, 0.0));
        assert!(c.triggered(0.5e-6, 0.02));
        assert!(!c.triggered(0.5e-6, 0.005));
    }

    #[test]
    fn oxidation_model_grows_linearly_after_breakaway() {
        let model = OxidationModel {
            parabolic: ScaleGrowth::new(1.0e-12),
            breakaway: BreakawayCriterion {
                critical_thickness: 1.0e-6,
                critical_strain: 1.0,
            },
            linear: LinearBreakawayRate::new(1.0e-9),
        };
        let x_short = model.thickness(1.0e2, 0.0);
        let x_long = model.thickness(1.0e10, 0.0);
        let ratio = x_long / x_short;
        // Long-time growth should be far in excess of sqrt(4) from
        // pure parabolic.
        assert!(ratio > 100.0);
    }
}