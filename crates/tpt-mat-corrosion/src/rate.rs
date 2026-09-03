//! Mixed-potential theory and corrosion rate.
//!
//! Wagner–Traud (1938) decomposed the corrosion current `i_corr` of
//! a coupled anode + cathode system into the intersection of the
//! anodic current density of the dissolution reaction with the
//! magnitude of the cathodic current density.  In Tafel-regime
//! (large `|η|` relative to the Tafel slopes) the intersection is
//!
//! ```text
//! β = b_a · b_c / (b_a + b_c)
//! i_corr = (i_0_a)^(b_c / (b_a + b_c))
//!           · (i_0_c)^(b_a / (b_a + b_c))
//!           · exp( (E_eq_c − E_eq_a) / (b_a + b_c) )
//! E_corr = β · ln(i_corr / i_0_a) + E_eq_a
//! ```

use serde::{Deserialize, Serialize};

use crate::butler_volmer::{butler_volmer_current_density, ElectrodeKinetics};

/// A coupled corrosion model: an anode + cathode + electrolyte.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CorrosionModel {
    /// Anodic half-reaction (metal dissolution).
    pub anode: ElectrodeKinetics,
    /// Cathodic half-reaction (e.g. hydrogen evolution / oxygen
    /// reduction).
    pub cathode: ElectrodeKinetics,
    /// Electrolyte pH (informational; affects E_eqs upstream).
    pub ph: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl CorrosionModel {
    /// Construct from an anode + cathode + pH + temperature.
    pub fn new(
        anode: ElectrodeKinetics,
        cathode: ElectrodeKinetics,
        ph: f64,
        temperature: f64,
    ) -> Self {
        Self {
            anode,
            cathode,
            ph,
            temperature,
        }
    }
}

/// Result of a corrosion-rate solve.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CorrosionRate {
    /// Corrosion current density (A/m²).
    pub current_density: f64,
    /// Corrosion potential (V vs SHE).
    pub corrosion_potential: f64,
    /// Penetration rate (mm/yr).
    pub penetration_rate_mm_per_yr: f64,
    /// Mass-loss rate (g/m²·day).
    pub mass_loss_rate_g_per_m2_day: f64,
}

/// Result of a mixed-potential solve.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MixedPotentialResult {
    /// Corrosion current density (A/m²) — positive scalar.
    pub corrosion_current_density: f64,
    /// Corrosion potential (V vs SHE).
    pub corrosion_potential: f64,
}

/// Solve the mixed-potential problem numerically: find the potential
/// `E` where the anodic current `i_a(E)` of the metal dissolution
/// equals the (magnitude of the) cathodic current `i_c(E)`.
pub fn mixed_potential(model: &CorrosionModel) -> MixedPotentialResult {
    let lo = model
        .cathode
        .equilibrium_potential
        .min(model.anode.equilibrium_potential);
    let hi = model
        .cathode
        .equilibrium_potential
        .max(model.anode.equilibrium_potential);
    let span = (hi - lo).abs().max(1.0);
    let mut lo = lo - 0.5 * span;
    let mut hi = hi + 0.5 * span;
    let mut f_lo = anodic_minus_cathodic(model, lo);
    let mut f_hi = anodic_minus_cathodic(model, hi);
    // Expand bounds if necessary.
    let mut expand = 0;
    while f_lo * f_hi > 0.0 && expand < 20 {
        lo -= 0.5 * span;
        hi += 0.5 * span;
        f_lo = anodic_minus_cathodic(model, lo);
        f_hi = anodic_minus_cathodic(model, hi);
        expand += 1;
    }
    // Bisection.
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let f_mid = anodic_minus_cathodic(model, mid);
        if f_mid.abs() < 1.0e-12 {
            lo = mid;
            hi = mid;
            break;
        }
        if f_lo * f_mid < 0.0 {
            hi = mid;
            f_hi = f_mid;
        } else {
            lo = mid;
            f_lo = f_mid;
        }
    }
    let e_corr = 0.5 * (lo + hi);
    let i_corr = butler_volmer_current_density(&model.anode, e_corr).abs();
    MixedPotentialResult {
        corrosion_current_density: i_corr,
        corrosion_potential: e_corr,
    }
}

/// `f(E) = i_a(E) + i_c(E)`; root when `i_a(E) = |i_c(E)|`.
fn anodic_minus_cathodic(model: &CorrosionModel, e: f64) -> f64 {
    let i_a = butler_volmer_current_density(&model.anode, e);
    let i_c = butler_volmer_current_density(&model.cathode, e);
    i_a + i_c
}

/// Solve for the corrosion rate of a coupled model.
///
/// The penetration rate is computed from the anodic reaction
/// (`valence = n`, `atomic_mass = M`, density `ρ`):
///
/// ```text
/// CR (m/s)  =  i_corr · M / (n · F · ρ)
/// CR (mm/yr) = CR · 1000 · 365.25 · 86400
/// ```
pub fn corrosion_rate(
    model: &CorrosionModel,
    atomic_mass_kg_per_mol: f64,
    valence: f64,
    density_kg_per_m3: f64,
) -> CorrosionRate {
    let mp = mixed_potential(model);
    let cr_m_per_s = mp.corrosion_current_density * atomic_mass_kg_per_mol
        / (valence * crate::FARADAY * density_kg_per_m3);
    let penetration_rate_mm_per_yr = cr_m_per_s * 1000.0 * 365.25 * 86400.0;
    // mass loss rate (g/m²·day) = ρ · CR (kg/m²/s) · 86400 s/day · 1000 g/kg
    let mass_loss_g_per_m2_day = cr_m_per_s * density_kg_per_m3 * 1.0e3 * 86400.0;
    CorrosionRate {
        current_density: mp.corrosion_current_density,
        corrosion_potential: mp.corrosion_potential,
        penetration_rate_mm_per_yr,
        mass_loss_rate_g_per_m2_day: mass_loss_g_per_m2_day,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::butler_volmer::ElectrodeKinetics;

    fn fe_model() -> CorrosionModel {
        // Iron dissolution vs hydrogen evolution at pH 7, 298 K.
        let anode = ElectrodeKinetics::from_alphas(
            -0.44,  // Fe → Fe²⁺ + 2e⁻
            1.0e-3, // i_0
            0.5, 0.5, 2.0, 298.0,
        );
        let cathode = ElectrodeKinetics::from_alphas(
            -0.41, // 2H⁺ + 2e⁻ → H₂ (pH 7 → −0.41 V)
            1.0e-1, 0.5, 0.5, 2.0, 298.0,
        );
        CorrosionModel::new(anode, cathode, 7.0, 298.0)
    }

    #[test]
    fn mixed_potential_between_eq_potentials() {
        let mp = mixed_potential(&fe_model());
        // Should lie between E_eq(anode) = −0.44 and E_eq(cathode) = −0.41.
        assert!(mp.corrosion_potential > -0.44 && mp.corrosion_potential < -0.41);
        assert!(mp.corrosion_current_density > 0.0);
    }

    #[test]
    fn corrosion_rate_reasonable_for_iron() {
        let cr = corrosion_rate(&fe_model(), 0.055_845, 2.0, 7874.0);
        // i_corr ~ 1e-3 A/m² ⇒ penetration ~0.5 mm/yr (Jones 1996).
        assert!(cr.current_density > 0.0);
        assert!(cr.penetration_rate_mm_per_yr > 0.0);
        assert!(cr.penetration_rate_mm_per_yr < 5.0);
    }

    #[test]
    fn corrosion_rate_scales_with_current_density() {
        let cr1 = corrosion_rate(&fe_model(), 0.055_845, 2.0, 7874.0);
        let cr2 = corrosion_rate(
            &CorrosionModel {
                anode: fe_model().anode,
                cathode: ElectrodeKinetics::from_alphas(-0.41, 1.0e-2, 0.5, 0.5, 2.0, 298.0),
                ph: 7.0,
                temperature: 298.0,
            },
            0.055_845,
            2.0,
            7874.0,
        );
        // Lower i_0(cathode) ⇒ smaller i_corr.
        assert!(cr2.current_density < cr1.current_density);
    }
}
