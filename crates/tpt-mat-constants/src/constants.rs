//! CODATA 2018 physical constants used across `tpt-materials`.

use serde::{Deserialize, Serialize};

/// Physical constants used across `tpt-materials`.
///
/// Values are CODATA 2018 recommended values, SI units.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PhysicalConstants;

impl PhysicalConstants {
    /// Boltzmann constant `k_B`.  `1.380 649e-23 J/K` (exact, since 2019).
    pub const BOLTZMANN: f64 = 1.380_649e-23;

    /// Ideal gas constant `R = N_A k_B`.  `8.314 462 618 J/(mol·K)` (exact).
    pub const GAS_CONSTANT: f64 = 8.314_462_618;

    /// Avogadro constant `N_A`.  `6.022 140 76e23 /mol` (exact).
    pub const AVOGADRO: f64 = 6.022_140_76e23;

    /// Faraday constant `F = N_A e`.  CODATA 2018 value, `96_485.33212 C/mol`.
    pub const FARADAY: f64 = 96_485.332_12;

    /// Planck constant `h`.  `6.626 070 15e-34 J·s` (exact).
    pub const PLANCK: f64 = 6.626_070_15e-34;

    /// Reduced Planck constant `ℏ = h / 2π`.  `1.054 571 817e-34 J·s`.
    pub const PLANCK_REDUCED: f64 = 1.054_571_817e-34;

    /// Elementary charge `e`.  `1.602 176 634e-19 C` (exact).
    pub const ELEMENTARY_CHARGE: f64 = 1.602_176_634e-19;

    /// Electron mass `m_e`.  `9.109_383_701_5e-31 kg` (CODATA 2018).
    pub const ELECTRON_MASS: f64 = 9.109_383_701_5e-31;

    /// Atomic mass constant `u`.  `1.660_539_066_60e-27 kg`.
    pub const ATOMIC_MASS_UNIT: f64 = 1.660_539_066_60e-27;

    /// Vacuum permittivity `ε₀`.  `8.854_187_812_8e-12 F/m`.
    pub const VACUUM_PERMITTIVITY: f64 = 8.854_187_812_8e-12;

    /// Vacuum permeability `μ₀`.  `1.256_637_062_12e-6 N/A²`.
    pub const VACUUM_PERMEABILITY: f64 = 1.256_637_062_12e-6;

    /// Speed of light in vacuum `c`.  `299 792 458 m/s` (exact).
    pub const SPEED_OF_LIGHT: f64 = 299_792_458.0;

    /// Stefan–Boltzmann constant `σ`.  `5.670 374 419e-8 W/(m²·K⁴)`.
    pub const STEFAN_BOLTZMANN: f64 = 5.670_374_419e-8;
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_testkit::assert_relative_eq;

    #[test]
    fn gas_constant_equals_avogadro_times_boltzmann() {
        // R is stored rounded to f64 precision; the relationship
        // R = N_A * k_B holds to better than 1e-6.
        assert_relative_eq!(
            PhysicalConstants::GAS_CONSTANT,
            PhysicalConstants::AVOGADRO * PhysicalConstants::BOLTZMANN,
            max_relative = 1e-6
        );
    }

    #[test]
    fn faraday_equals_avogadro_times_elementary_charge() {
        assert_relative_eq!(
            PhysicalConstants::FARADAY,
            PhysicalConstants::AVOGADRO * PhysicalConstants::ELEMENTARY_CHARGE,
            max_relative = 1e-3
        );
    }
}
