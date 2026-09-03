//! Example: SiC / Al composite effective thermal properties.
//!
//! Sweeps the SiC inclusion fraction `f_SiC ∈ [0, 0.5]` and
//! computes the effective thermal conductivity and CTE under
//! several homogenization bounds.
//!
//! Material data (room-temperature):
//!   Al:  k = 237 W/(m·K),  CTE = 23.1e-6 /K
//!   SiC: k = 120 W/(m·K), CTE = 4.0e-6  /K
//!   Al:  K = 76 GPa,    Al: G = 26 GPa
//!   SiC: K = 220 GPa,   SiC: G = 180 GPa
//!
//! Run with `cargo run -p composite-thermal-properties`.

use tpt_mat_thermal::{effective_conductivity, effective_cte, ConductivityBound, CteBound};

const K_AL: f64 = 237.0;
const K_SIC: f64 = 120.0;
const CTE_AL: f64 = 23.1e-6;
const CTE_SIC: f64 = 4.0e-6;
const KMOD_AL: f64 = 76.0e9;
const GMOD_AL: f64 = 26.0e9;
const KMOD_SIC: f64 = 220.0e9;
const GMOD_SIC: f64 = 180.0e9;

fn main() {
    println!("SiC / Al composite effective thermal properties");
    println!("================================================");
    println!("  Al  : k = {K_AL:>5} W/m·K, CTE = {CTE_AL:.2e} /K, K = {KMOD_AL:.2e} Pa");
    println!("  SiC : k = {K_SIC:>5} W/m·K, CTE = {CTE_SIC:.2e} /K, K = {KMOD_SIC:.2e} Pa\n");

    println!("  f_SiC    Voigt      Reuss      VRH        HS_lo      HS_up      MG");
    println!("  -----    -------    -------    -------    -------    -------    -------");
    let fractions = [0.0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50];
    for &f in &fractions {
        let k_arr = [K_AL, K_SIC];
        let ff = [1.0 - f, f];
        let v = effective_conductivity(ConductivityBound::Voigt, &k_arr, &ff);
        let r = effective_conductivity(ConductivityBound::Reuss, &k_arr, &ff);
        let vrh = effective_conductivity(ConductivityBound::VoigtReussHill, &k_arr, &ff);
        let lo = effective_conductivity(ConductivityBound::HashinShtrikmanLower, &k_arr, &ff);
        let up = effective_conductivity(ConductivityBound::HashinShtrikmanUpper, &k_arr, &ff);
        let mg = effective_conductivity(ConductivityBound::MaxwellGarnett, &k_arr, &ff);
        println!(
            "  {f:>4.2}    {v:>7.2}    {r:>7.2}    {vrh:>7.2}    {lo:>7.2}    {up:>7.2}    {mg:>7.2}"
        );
    }

    println!("\nEffective CTE (1/K) vs f_SiC:");
    println!("  f_SiC    Turner     Kerner    RH_lower   RH_upper");
    println!("  -----    --------   --------  --------   --------");
    for &f in &fractions {
        let cte_arr = [CTE_AL, CTE_SIC];
        let k_arr = [KMOD_AL, KMOD_SIC];
        let g_arr = [GMOD_AL, GMOD_SIC];
        let ff = [1.0 - f, f];
        let turner = effective_cte(CteBound::Turner, &cte_arr, &k_arr, &g_arr, &ff);
        let kerner = effective_cte(CteBound::Kerner, &cte_arr, &k_arr, &g_arr, &ff);
        let rh_l = effective_cte(CteBound::RosenHashinLower, &cte_arr, &k_arr, &g_arr, &ff);
        let rh_u = effective_cte(CteBound::RosenHashinUpper, &cte_arr, &k_arr, &g_arr, &ff);
        println!("  {f:>4.2}    {turner:.3e}   {kerner:.3e}   {rh_l:.3e}   {rh_u:.3e}");
    }

    println!("\nNotes:");
    println!("  - k_eff decreases as SiC (lower k than Al) is added.");
    println!("  - CTE drops sharply with SiC fraction (low CTE phase dominates the bulk-weighted Turner bound).");
    println!("  - Hashin-Shtrikman bounds enclose the Maxwell-Garnett estimate and lie inside Voigt/Reuss.");
}
