//! Recrystallization kinetics.
//!
//! - [`JmakRecrystallization`] — JMAK-for-RX: `X(t) = 1 - exp(-f(t))`
// with `f(t) = (k_G · t)^n` for static RX (Avrami form).
//! - [`RxxKinetics`] — coupled dynamic recrystallization (DRX) with a
//! stored-energy-driven nucleation-and-growth driver.
//! - [`ZenerHollomon`] — `Z = ε̇ · exp(Q / (R T))` parameter and the
//! classic Sellars-Tegart constitutive relation for steady-state
//! stress.

use serde::{Deserialize, Serialize};

/// Static recrystallization kinetics following the JMAK (Avrami)
/// form `X(t) = 1 − exp(−k · t^n)` with Arrhenius temperature
/// dependence `k(T) = k_0 · exp(−Q / (R T))`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct JmakRecrystallization {
    /// Pre-exponential `k_0` (1/s).
    pub k_0: f64,
    /// Activation energy `Q` (J/mol).
    pub activation_energy: f64,
    /// Avrami exponent `n` (typically 2–4).
    pub avrami_n: f64,
}

impl Default for JmakRecrystallization {
    fn default() -> Self {
        Self {
            k_0: 1.0e8,
            activation_energy: 180.0e3,
            avrami_n: 2.0,
        }
    }
}

impl JmakRecrystallization {
    /// Construct.
    pub fn new(k_0: f64, activation_energy: f64, avrami_n: f64) -> Self {
        Self {
            k_0,
            activation_energy,
            avrami_n,
        }
    }

    /// Recrystallized volume fraction `X(t, T)`.
    pub fn recrystallized_fraction(
        &self,
        time_s: f64,
        temperature_k: f64,
        gas_constant: f64,
    ) -> f64 {
        let k_t = self.k_0
            * (-self.activation_energy / (gas_constant * temperature_k)).exp();
        1.0 - (-(k_t * time_s).powf(self.avrami_n)).exp()
    }

    /// Time to a given recrystallized fraction at temperature `T`.
    pub fn time_to_fraction(
        &self,
        x_target: f64,
        temperature_k: f64,
        gas_constant: f64,
    ) -> f64 {
        let k_t = self.k_0
            * (-self.activation_energy / (gas_constant * temperature_k)).exp();
        if k_t <= 0.0 || x_target <= 0.0 || x_target >= 1.0 {
            return f64::INFINITY;
        }
        ((-x_target.ln() / 1.0).log(std::f64::consts::E)
            / self.avrami_n)
        .exp()
            / k_t
    }

    /// t_{50} (time to 50 % RX).
    pub fn t_50(&self, temperature_k: f64, gas_constant: f64) -> f64 {
        self.time_to_fraction(0.5, temperature_k, gas_constant)
    }
}

/// Zener–Hollomon parameter `Z = ε̇ · exp(Q / (R T))` and the
/// Sellars–Tegart steady-state constitutive relation
/// `σ_s = (1 / α) · asinh((Z / A)^(1/n))`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ZenerHollomon {
    /// Strain rate (1/s).
    pub strain_rate: f64,
    /// Activation energy for hot working (J/mol).
    pub activation_energy: f64,
    /// Sellars-Tegart `α` (1/MPa).
    pub alpha: f64,
    /// Sellars-Tegart `n`.
    pub n: f64,
    /// Sellars-Tegart `A` (1/s).
    pub a: f64,
}

impl Default for ZenerHollomon {
    fn default() -> Self {
        Self {
            strain_rate: 1.0,
            activation_energy: 312.0e3,
            alpha: 0.014,
            n: 4.3,
            a: 1.0e8,
        }
    }
}

impl ZenerHollomon {
    /// `Z(ε̇, T)` (units of 1/s).
    pub fn z(&self, temperature_k: f64, gas_constant: f64) -> f64 {
        self.strain_rate * (self.activation_energy / (gas_constant * temperature_k)).exp()
    }

    /// Steady-state flow stress (MPa).
    pub fn steady_state_stress(self, temperature_k: f64, gas_constant: f64) -> f64 {
        let z = self.z(temperature_k, gas_constant);
        let arg = (z / self.a).powf(1.0 / self.n);
        (1.0 / self.alpha) * arg.asinh()
    }
}

/// Dynamic recrystallization driver: maintains a stored-energy
/// fraction that is consumed by nucleation of new grains, each of
/// which grows at rate `G(T)`.  The recrystallized volume follows
/// `X(t) = 1 − exp(−π · Ṅ · G³ · t⁴ / 3)` (Cahn–Hagel-style),
/// where `Ṅ(T)` is a temperature-dependent nucleation rate.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DrxKinetics {
    /// Nucleation rate pre-exponential `Ṅ_0` (m⁻³·s⁻¹).
    pub nucleation_rate_0: f64,
    /// Nucleation activation energy `Q_N` (J/mol).
    pub nucleation_q: f64,
    /// Growth-rate pre-exponential `G_0` (m/s).
    pub growth_rate_0: f64,
    /// Growth activation energy `Q_G` (J/mol).
    pub growth_q: f64,
    /// Initial stored-energy fraction (dimensionless, 0–1).
    pub stored_energy_initial: f64,
}

impl Default for DrxKinetics {
    fn default() -> Self {
        Self {
            nucleation_rate_0: 1.0e6,
            nucleation_q: 200.0e3,
            growth_rate_0: 1.0e-6,
            growth_q: 150.0e3,
            stored_energy_initial: 1.0,
        }
    }
}

impl DrxKinetics {
    /// `Ṅ(T)` (m⁻³·s⁻¹).
    pub fn nucleation_rate(&self, temperature_k: f64, gas_constant: f64) -> f64 {
        self.nucleation_rate_0
            * (-self.nucleation_q / (gas_constant * temperature_k)).exp()
    }

    /// `G(T)` (m/s).
    pub fn growth_rate(&self, temperature_k: f64, gas_constant: f64) -> f64 {
        self.growth_rate_0 * (-self.growth_q / (gas_constant * temperature_k)).exp()
    }

    /// Recrystallized volume fraction at time `t` (Cahn–Hagel form).
    pub fn recrystallized_fraction(
        &self,
        time_s: f64,
        temperature_k: f64,
        gas_constant: f64,
    ) -> f64 {
        let n = self.nucleation_rate(temperature_k, gas_constant);
        let g = self.growth_rate(temperature_k, gas_constant);
        let expo = std::f64::consts::PI * n * g.powi(3) * time_s.powi(4) / 3.0;
        1.0 - (-expo).exp()
    }

    /// Remaining stored energy fraction (`1 − X`).
    pub fn stored_energy_fraction(
        &self,
        time_s: f64,
        temperature_k: f64,
        gas_constant: f64,
    ) -> f64 {
        self.stored_energy_initial
            * (1.0 - self.recrystallized_fraction(time_s, temperature_k, gas_constant))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn jmak_fraction_monotone_in_time() {
        let r = JmakRecrystallization::default();
        let gas = 8.314_462_618;
        let t = 900.0;
        let x1 = r.recrystallized_fraction(10.0, t, gas);
        let x2 = r.recrystallized_fraction(1000.0, t, gas);
        assert!(x2 > x1);
    }

    #[test]
    fn jmak_higher_temperature_faster_rxn() {
        let r = JmakRecrystallization::default();
        let gas = 8.314_462_618;
        let x_lo = r.recrystallized_fraction(100.0, 700.0, gas);
        let x_hi = r.recrystallized_fraction(100.0, 900.0, gas);
        assert!(x_hi > x_lo);
    }

    #[test]
    fn jmak_t50_matches_half_fraction() {
        let r = JmakRecrystallization::default();
        let gas = 8.314_462_618;
        let t = 900.0;
        let t50 = r.t_50(t, gas);
        let x = r.recrystallized_fraction(t50, t, gas);
        assert!(approx(x, 0.5, 1.0e-6));
    }

    #[test]
    fn zener_hollomon_increases_with_strain_rate() {
        let mut z = ZenerHollomon::default();
        z.strain_rate = 0.1;
        let gas = 8.314_462_618;
        let z1 = z.z(1000.0, gas);
        z.strain_rate = 10.0;
        let z2 = z.z(1000.0, gas);
        assert!(z2 > z1);
    }

    #[test]
    fn drx_fraction_monotone_in_time() {
        let d = DrxKinetics {
            nucleation_rate_0: 1.0e20,
            nucleation_q: 100.0e3,
            growth_rate_0: 1.0e-3,
            growth_q: 80.0e3,
            stored_energy_initial: 1.0,
        };
        let gas = 8.314_462_618;
        let x1 = d.recrystallized_fraction(0.01, 1000.0, gas);
        let x2 = d.recrystallized_fraction(1.0, 1000.0, gas);
        assert!(x2 > x1);
    }

    #[test]
    fn drx_stored_energy_decreases_as_rxn_proceeds() {
        let d = DrxKinetics {
            nucleation_rate_0: 1.0e20,
            nucleation_q: 100.0e3,
            growth_rate_0: 1.0e-3,
            growth_q: 80.0e3,
            stored_energy_initial: 1.0,
        };
        let gas = 8.314_462_618;
        let s1 = d.stored_energy_fraction(0.01, 1000.0, gas);
        let s2 = d.stored_energy_fraction(1.0, 1000.0, gas);
        assert!(s2 < s1);
    }
}