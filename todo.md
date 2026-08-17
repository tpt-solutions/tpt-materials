# tpt-materials — Project Todo

Organization: **TPT Solutions** · License: **MIT OR Apache-2.0** (dual)

---

## Phase 0 — Repo Scaffolding

- [ ] Workspace `Cargo.toml`
- [ ] `.gitignore`
- [ ] `LICENSE-MIT` (TPT Solutions copyright)
- [ ] `LICENSE-APACHE` (TPT Solutions copyright)
- [ ] `README.md` (spec §13 template)
- [ ] `CONTRIBUTING.md` (fork → branch → code+tests → fmt/clippy/test → `cargo deny check licenses` → PR + DCO sign-off → RFC discussion → 2 approvals)
- [ ] `SECURITY.md` (private disclosure process)
- [ ] `CODE_OF_CONDUCT.md`
- [ ] `CHANGELOG.md`
- [ ] `deny.toml` (allow MIT/Apache-2.0/BSD-2/BSD-3/ISC/Zlib/Unicode-3.0; `copyleft = "deny"`; `unlicensed = "deny"`)
- [ ] `rustfmt.toml`
- [ ] `clippy.toml`
- [ ] `.github/workflows/ci.yml`
- [ ] `.github/workflows/license.yml`
- [ ] `.github/workflows/benchmark.yml`
- [ ] `.github/workflows/docs.yml`
- [ ] `.github/workflows/release.yml`
- [ ] `.github/ISSUE_TEMPLATE/`
- [ ] `.github/PULL_REQUEST_TEMPLATE.md`
- [ ] Directory skeleton: `crates/`, `examples/`, `test-data/{ebsd,crystal-structures,phase-diagrams,golden}/`, `benches/`, `docs/{book,rfc,api}/`, `rfcs/`
- [ ] Public GitHub Projects roadmap board
- [ ] Trademark note: "TPT Materials" name reserved by TPT

---

## Phase 1 — Foundation (Months 1-3)

**Crates:** `tpt-mat-core`, `tpt-mat-constants`, `tpt-mat-crystallography`, `tpt-mat-wasm` (scaffold)
**Substrate:** `tpt-math-linalg-fixed` (Schmid tensors, elasticity tensors)

- [ ] `tpt-mat-core`
  - [ ] `MaterialMicrostructure` (id, name, phases, grains, volume_element, temperature)
  - [ ] `Phase` (id, name, crystal_structure, composition, volume_fraction, properties)
  - [ ] `Grain` (id, phase, orientation, centroid, equivalent_radius, neighbors)
  - [ ] `Composition` + `CompositionBasis` (Atomic/Weight/Mole fraction)
  - [ ] `CrystalOrientation` + `OrientationRepresentation` (EulerBunge, Quaternion, RotationMatrix, Rodrigues, AxisAngle)
- [ ] `tpt-mat-crystallography`
  - [ ] `CrystalStructure` enum (FCC, BCC, HCP, Diamond, SimpleCubic, BCT, Custom) + `LatticeParameters`
  - [ ] `MillerIndex`
  - [ ] `SlipSystem` (slip_direction, slip_plane_normal, critical_resolved_shear_stress)
  - [ ] `CrystalStructure::slip_systems()` for FCC (12: {111}<110>), BCC (12: {110}<111>), HCP (basal/prismatic/pyramidal)
  - [ ] `SlipSystem::schmid_tensor()`
  - [ ] `SlipSystem::resolved_shear_stress()`
- [ ] `tpt-mat-constants`
  - [ ] `PhysicalConstants` (Boltzmann, gas constant, Avogadro, Faraday, Planck)
  - [ ] `PeriodicTable::atomic_mass()`
  - [ ] `PeriodicTable::atomic_radius()`
- [ ] `tpt-mat-wasm` — crate scaffold only (full bindings in later phases)
- [ ] Verification test: FCC slip system count == 12
- [ ] Verification test: Schmid tensor symmetry
- [ ] Example scaffold: `fcc-single-crystal-tension` (data only, solver comes Phase 2)

**Milestone:** Calculate Schmid tensors for FCC/BCC/HCP

---

## Phase 2 — Crystal Plasticity (Months 4-6)

**Crates:** `tpt-mat-crystal-plasticity`, `tpt-mat-hardening`, `tpt-mat-texture`
**Substrate:** `tpt-fem` (mesh, elements, periodic BCs, assembly, Newton-Raphson solver)

- [ ] `tpt-mat-crystal-plasticity`
  - [ ] `CrystalPlasticityModel` (crystal_structure, slip_systems, hardening_law, rate_sensitivity, reference_strain_rate, elastic_tensor)
  - [ ] `HardeningLaw` enum: Voce, PowerLaw, KocksMecking, Combined (latent hardening matrix)
  - [ ] `CrystalPlasticityModel::plastic_velocity_gradient()` (viscoplastic power-law flow rule)
  - [ ] `CrystalPlasticityModel::update_hardening()`
  - [ ] `PlasticIncrement` (velocity_gradient, slip_rates)
  - [ ] `CpFemSolver` (mesh, material, boundary_conditions)
  - [ ] `CpFemSolver::solve_increment()` (assemble stiffness → Newton-Raphson → update stress/hardening/lattice rotation)
  - [ ] `CpFemResult` (stresses, strains, slip_rates, lattice_rotations, accumulated_shear)
- [ ] `tpt-mat-hardening` — Voce and power-law hardening implementations wired to `HardeningLaw`
- [ ] `tpt-mat-texture`
  - [ ] `TextureAnalyzer` (orientations, weights)
  - [ ] `TextureAnalyzer::from_ebsd()`
  - [ ] `TextureAnalyzer::pole_figure()` → `PoleFigure`
  - [ ] `TextureAnalyzer::orientation_distribution_function()` (ODF via KDE)
  - [ ] `TextureAnalyzer::taylor_factor()`
- [ ] RFC 0001: `rfcs/0001-crystal-plasticity-fem.md`
- [ ] Golden test data: `test-data/golden/crystal-plasticity/fcc-single-crystal-tension.json`
- [ ] Golden test data: `test-data/golden/crystal-plasticity/bcc-polycrystal-rve.json`
- [ ] Golden test data: `test-data/golden/crystal-plasticity/taylor-factor-random.json`
- [ ] Verification test: FCC random-texture Taylor factor ≈ 3.06
- [ ] Example: `examples/fcc-single-crystal-tension/`

**Milestone:** Single crystal tension test with correct slip activation

---

## Phase 3 — Phase-Field (Months 7-9)

**Crates:** `tpt-mat-phase-field`, `tpt-mat-grain-growth`, `tpt-mat-solidification`
**Substrate:** `tpt-science` (Allen-Cahn, Cahn-Hilliard equations)

- [ ] `tpt-mat-phase-field`
  - [ ] `PhaseFieldSolver` (grid, model, time_step, mobility, gradient_coefficient)
  - [ ] `PhaseFieldModel` enum: AllenCahn, CahnHilliard, Kobayashi (latent heat), MultiPhase (num_grains)
  - [ ] `FreeEnergyFunctional` + `BulkEnergy` (DoubleWell, Polynomial, RegularSolution)
  - [ ] `PhaseFieldSolver::step_allen_cahn()`
  - [ ] `PhaseFieldSolver::step_cahn_hilliard()`
  - [ ] `PhaseFieldSolver::compute_laplacian()`
  - [ ] `PhaseFieldResult` (order_parameter, concentration, free_energy, interface_area)
- [ ] `tpt-mat-grain-growth`
  - [ ] `GrainGrowthSolver` (phase_field, mobilities)
  - [ ] `GrainBoundaryMobility` (base_mobility, misorientation_dependent)
  - [ ] `GrainGrowthSolver::simulate_growth()`
  - [ ] `GrainGrowthSolver::grain_size_distribution()` → `GrainSizeDistribution`
- [ ] `tpt-mat-solidification`
  - [ ] `SolidificationSolver` (phase_field, thermal_field, anisotropy)
  - [ ] `AnisotropyModel` (strength, mode)
  - [ ] `SolidificationSolver::simulate_dendrite()`
  - [ ] `SolidificationSolver::secondary_arm_spacing()`
  - [ ] `SolidificationResult` (solid_fraction, tip_velocity, primary/secondary arm spacing, microstructure)
- [ ] RFC 0002: `rfcs/0002-phase-field-framework.md`
- [ ] Golden test data: `allen-cahn-grain-growth.json`, `cahn-hilliard-spinodal.json`, `kobayashi-dendrite.json`
- [ ] Verification test: free energy monotonically decreases under Allen-Cahn step
- [ ] Example: `examples/spinodal-decomposition/`
- [ ] Example: `examples/dendritic-solidification/`

**Milestone:** Simulate spinodal decomposition and dendritic solidification

---

## Phase 4 — Diffusion & Transformation (Months 10-12)

**Crates:** `tpt-mat-diffusion`, `tpt-mat-phase-transform`, `tpt-mat-calphad`
**Substrate:** `tpt-science` (Fick's laws, Darken's equation)

- [ ] `tpt-mat-diffusion`
  - [ ] `DiffusionSolver` (grid, diffusivity, boundary_conditions)
  - [ ] `DiffusivityModel` enum: Constant, Arrhenius (D0, activation energy), CompositionDependent
  - [ ] `DiffusivityModel::at_temperature()`
  - [ ] `DiffusionSolver::solve_transient()` (Fick's second law)
  - [ ] `DiffusionSolver::solve_steady_state()`
  - [ ] `DiffusionResult` (concentration, flux, total_diffused)
- [ ] `tpt-mat-phase-transform`
  - [ ] `PhaseTransformation` (ttt_diagram, cct_diagram, kinetics)
  - [ ] `TransformationKinetics` enum: JohnsonMehlAvrami, KoistinenMarburger
  - [ ] `PhaseTransformation::transformed_fraction()`
  - [ ] `TttDiagram`, `CctDiagram`
- [ ] `tpt-mat-calphad`
  - [ ] `CalphadDatabase` (phases, elements)
  - [ ] `CalphadPhase` + `GibbsEnergyModel` + `InteractionParameter`
  - [ ] `CalphadDatabase::gibbs_energy()`
  - [ ] `CalphadDatabase::equilibrium()` → `PhaseEquilibrium`
  - [ ] `CalphadDatabase::phase_diagram()` → `BinaryPhaseDiagram`
- [ ] Golden test data: `fickian-diffusion-couple.json`, `arrhenius-temperature-dependence.json`
- [ ] Benchmark: `benches/diffusion-couple.rs`

**Milestone:** Simulate diffusion couple and phase transformation

---

## Phase 5 — Micro-Mechanics (Months 13-15)

**Crates:** `tpt-mat-rve`, `tpt-mat-homogenization`, `tpt-mat-composite-micro`
**Substrate:** `tpt-fem` (Voxel/hex meshing, periodic BCs)

- [ ] `tpt-mat-rve`
  - [ ] `RepresentativeVolumeElement` (mesh, grains, boundary_conditions, material)
  - [ ] `RveBoundaryCondition` enum: Periodic, KinematicUniform (Taylor), StaticUniform (Sachs)
  - [ ] `RepresentativeVolumeElement::generate_voronoi()`
  - [ ] `RepresentativeVolumeElement::generate_from_ebsd()`
  - [ ] `RepresentativeVolumeElement::homogenize()` (Hill-Mandel condition)
  - [ ] `HomogenizedResponse` (average_stress, average_strain, tangent_modulus, local_stresses, local_slip_rates)
- [ ] `tpt-mat-homogenization`
  - [ ] `HomogenizationMethod` enum: Voigt, Reuss, Hill, MoriTanaka, SelfConsistent, FftBased, FemBased
  - [ ] `Homogenizer::effective_stiffness()` — Voigt average
  - [ ] `Homogenizer::effective_stiffness()` — Reuss average
  - [ ] `Homogenizer::effective_stiffness()` — Mori-Tanaka
  - [ ] `Homogenizer::effective_stiffness()` — self-consistent
  - [ ] `Homogenizer::effective_stiffness()` — FFT-based (Moulinec-Suquet)
- [ ] `tpt-mat-composite-micro`
  - [ ] `CompositeMicroMechanics` (fiber, matrix, fiber_volume_fraction, fiber_arrangement)
  - [ ] `FiberArrangement` enum: Unidirectional, Woven, Random
  - [ ] `CompositeMicroMechanics::rule_of_mixtures_stiffness()`
  - [ ] `CompositeMicroMechanics::halpin_tsai()`
  - [ ] `CompositeMicroMechanics::fiber_misorientation_effect()`
- [ ] RFC 0003: `rfcs/0003-rve-periodic-boundary.md`
- [ ] Golden test data: `voronoi-rve-homogenization.json`, `mori-tanaka-stiffness.json`, `composite-rule-of-mixtures.json`
- [ ] Verification test: Hill-Mandel condition (σ̄:ε̄ == volume-averaged σ:ε)
- [ ] Benchmark: `benches/rve-homogenization.rs`, `benches/crystal-plasticity-fem.rs`
- [ ] Example: `examples/polycrystal-rve-homogenization/`

**Milestone:** Polycrystal RVE homogenization with Hill-Mandel verification

---

## Phase 6 — Degradation (Months 16-18)

**Crates:** `tpt-mat-damage`, `tpt-mat-fatigue-micro`, `tpt-mat-corrosion`
**Substrate:** `tpt-science` (electrochemistry, corrosion kinetics)

- [ ] `tpt-mat-damage`
  - [ ] `GursonTvergaardNeedleman` (f_0, f_c, f_f, q_1, q_2, q_3, f_n, s_n, e_n)
  - [ ] `GursonTvergaardNeedleman::yield_function()`
  - [ ] `GursonTvergaardNeedleman::update_porosity()` (void growth + nucleation)
- [ ] `tpt-mat-fatigue-micro`
  - [ ] `MicrostructuralFatigue` (rve, criterion)
  - [ ] `FatigueCriterion` enum: Findley, FatemiSocie, SmithWatsonTopper, CrystallographicSlip
  - [ ] `MicrostructuralFatigue::predict_crack_initiation()` → `CrackInitiationResult`
  - [ ] `MicrostructuralFatigue::fatigue_indicator_parameter()`
- [ ] `tpt-mat-corrosion`
  - [ ] `CorrosionModel` (anode, cathode, electrolyte)
  - [ ] `ElectrodeKinetics` (exchange_current_density, tafel_slope, equilibrium_potential)
  - [ ] `CorrosionModel::corrosion_rate()` (Butler-Volmer)
  - [ ] `CorrosionModel::polarization_curve()`
- [ ] Golden test data: `gtn-void-growth.json`, `fatigue-crack-initiation.json`, `corrosion-polarization.json`

**Milestone:** Predict fatigue crack initiation site in polycrystal

---

## Phase 7 — Energy Materials & Manufacturing (Months 19-21)

**Crates:** `tpt-mat-battery`, `tpt-mat-hydrogen-storage`, `tpt-mat-polymer`, `tpt-mat-hydrogel`, `tpt-mat-additive`, `tpt-mat-welding`, `tpt-mat-heat-treatment`

- [ ] `tpt-mat-battery`
  - [ ] `BatteryElectrodeModel` (active_material, degradation_mechanisms)
  - [ ] `ActiveMaterial` (chemistry, particle_radius, diffusion_coefficient, partial_molar_volume)
  - [ ] `DegradationMechanism` enum: SeiGrowth, ParticleCracking, LithiumPlating, TransitionMetalDissolution
  - [ ] `BatteryElectrodeModel::simulate_diffusion_stress()`
  - [ ] `BatteryElectrodeModel::capacity_fade_curve()` → `DegradationCurve`
- [ ] `tpt-mat-hydrogen-storage`
  - [ ] `HydrogenStorageMaterial` (hydride_type, storage_capacity_wt_pct, absorption_kinetics)
  - [ ] `HydrideType` enum: MetalHydride, ChemicalHydride, PorousMaterial
  - [ ] `HydrogenStorageMaterial::pvc_isotherm()`
- [ ] `tpt-mat-polymer`
  - [ ] `PolymerModel` (chain_model, crosslink_density)
  - [ ] `ChainModel` enum: FreelyJointedChain, WormLikeChain, ArrudaBoyce
  - [ ] `PolymerModel::stress_strain()` (Arruda-Boyce 8-chain, inverse Langevin function)
- [ ] `tpt-mat-hydrogel` — soft-matter hydrogel swelling/mechanics model
- [ ] `tpt-mat-additive`
  - [ ] `AdditiveManufacturingModel` (process, material)
  - [ ] `AmProcess` enum: LaserPowderBedFusion, DirectedEnergyDeposition, ElectronBeamMelting
  - [ ] `AdditiveManufacturingModel::thermal_history()`
  - [ ] `AdditiveManufacturingModel::predict_microstructure()` → `PredictedMicrostructure`
  - [ ] `AdditiveManufacturingModel::residual_stress()` → `ResidualStressField`
- [ ] `tpt-mat-welding`
  - [ ] `WeldModel` (base_metal, filler_metal, process)
  - [ ] `WeldModel::heat_affected_zone()` → `HazResult`
  - [ ] `WeldModel::predict_haz_microstructure()` → `HazMicrostructure`
- [ ] `tpt-mat-heat-treatment`
  - [ ] `HeatTreatmentSimulator` (alloy, phase_diagram)
  - [ ] `HeatTreatmentProcess` enum: Annealing, Quenching, Tempering, Aging, SolutionTreatment
  - [ ] `HeatTreatmentSimulator::simulate()` → `HeatTreatmentResult`
  - [ ] `HeatTreatmentSimulator::predict_hardness()`
- [ ] RFC 0004: `rfcs/0004-battery-degradation-model.md`
- [ ] RFC 0005: `rfcs/0005-additive-manufacturing.md`
- [ ] Golden test data: `battery-sei-growth.json`, `electrode-capacity-fade.json`
- [ ] Example: `examples/battery-electrode-degradation/`
- [ ] Example: `examples/additive-manufacturing-microstructure/`

**Milestone:** Generate battery capacity fade curve for `tpt-energy`

---

## Phase 8 — Informatics & Ecosystem (Months 22-24)

**Crates:** `tpt-mat-database`, `tpt-mat-machine-learning`

- [ ] `tpt-mat-database`
  - [ ] `MaterialsDatabase` (materials: HashMap<MaterialId, MaterialRecord>)
  - [ ] `MaterialRecord` (name, composition, mechanical, thermal, electrical, sources)
  - [ ] `MaterialsDatabase::search_by_property()`
  - [ ] `MaterialsDatabase::load_builtin()`
- [ ] `tpt-mat-machine-learning`
  - [ ] `MaterialsMl` (model)
  - [ ] `MlModel` enum: PropertyPredictor, PhasePredictor, SurrogateModel
  - [ ] `MaterialsMl::train()`
  - [ ] `MaterialsMl::predict()`
- [ ] RFC 0006: `rfcs/0006-calphad-integration.md` (if not already closed in Phase 4)
- [ ] `tpt-mat-wasm` — full bindings
  - [ ] `WasmRveSolver` (`new()`, `homogenize()`)
  - [ ] `WasmPhaseField` (`step()`, `get_order_parameter()`)
  - [ ] Interactive RVE explorer demo (browser)
  - [ ] Browser-based phase-field demo
- [ ] Cross-repo integration: export battery degradation curves to `tpt-energy`
- [ ] Cross-repo integration: export composite fatigue / alloy creep to `tpt-transport`
- [ ] Cross-repo integration: export solder joint reliability / wire bond fatigue to `tpt-electronics`
- [ ] Cross-repo integration: export implant osseointegration / corrosion to `tpt-medical`
- [ ] Validation standards checklist: ASTM E8/E8M, ASTM E112, ASTM E2103, ISO 6892-1, NIST reference microstructures
- [ ] `docs/book/` — mdBook user guide complete
- [ ] `docs/api/` — full API docs published

**Milestone:** End-to-end micro-to-macro pipeline (microstructure → device property)

---

## Ongoing / Cross-Phase

- [ ] Maintain `cargo fmt` / `cargo clippy` / `cargo test` passing on every PR
- [ ] Maintain `cargo deny check licenses` passing (MIT chain enforcement, spec §8)
- [ ] SemVer releases on 6-week cadence
- [ ] RFC discussion required for each new constitutive model
- [ ] 2-approval merge policy maintained
