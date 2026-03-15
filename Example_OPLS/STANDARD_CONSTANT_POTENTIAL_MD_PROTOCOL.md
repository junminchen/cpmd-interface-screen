# Standard Constant-Potential MD Simulation Protocol

## 1. Purpose

This protocol standardizes constant-potential interface simulations across many electrolyte formulations so that extracted descriptors are comparable across the dataset.

The main goal is not to optimize each formulation separately.

The main goal is to keep the simulation protocol fixed enough that descriptor differences are driven by formulation chemistry rather than by changing simulation details.

## 2. Principle

When screening many formulations, the following should stay fixed whenever possible:

- electrode model
- cell geometry
- temperature
- timestep and thermostat settings
- electrostatic treatment
- constant-potential solver settings
- trajectory reporting interval
- analysis window definitions

Only formulation chemistry and the planned voltage schedule should vary.

## 3. Recommended standard protocol

### 3.1 Electrode model

- electrode type: fixed multilayer electrode slabs
- layers per electrode: `5`
- interlayer spacing: `2.5 A`
- in-plane spacing: `3.5 A`
- electrode atoms fixed during MD
- use charge constraint: `true`
- charge constraint target: `0 e`

Reason:

- fixed geometry removes electrode drift
- charge-neutral total electrode state improves consistency across formulations

### 3.2 Cell geometry

Use one common interface box for the first screening phase.

Recommended default:

- `box_x_angstrom = 55.0`
- `box_y_angstrom = 55.0`
- `box_z_angstrom = 90.0`
- `z_cathode_surface_angstrom = 16.0`
- `z_anode_surface_angstrom = 74.0`
- `z_liq_min_angstrom = 20.0`
- `z_liq_max_angstrom = 70.0`

Reason:

- keeps interfacial area identical
- simplifies direct comparison of surface densities and capacitances

If formulation density differs strongly, change composition counts rather than changing box geometry in the first pass.

### 3.3 OpenMM settings

Recommended defaults:

- nonbonded method: `PME`
- nonbonded cutoff: `1.0 nm`
- constraints: `HBonds`
- rigid water: `false`
- ewald error tolerance: `5e-4`
- remove center-of-mass motion: `true` if supported in the workflow

### 3.4 Thermostat / ensemble

Recommended defaults:

- temperature: `298.15 K`
- pressure: `1.0 bar` for bulk pre-equilibration only
- interface production ensemble: `NVT`
- timestep: `2.0 fs`
- friction coefficient: `1.0 ps^-1`

Important:

- do **not** use isotropic NPT during interface production for slab systems
- use barostat only in bulk reference simulations

### 3.5 Constant-potential solver settings

Recommended defaults:

- Gaussian width: `0.20 nm`
- Thomas-Fermi scale: `5.0 nm^-1`
- conjugate-gradient error tolerance: `1e-4`

These should remain fixed across the whole screening campaign.

### 3.6 Voltage schedule

For descriptor screening, use one common voltage grid for all formulations.

Recommended first-pass grid:

- `0.5 V`
- `1.0 V`
- `1.5 V`
- `2.0 V`

Optional extension:

- `3.0 V`

Reason:

- enough points to capture descriptor trends versus voltage
- still affordable for a medium-size screening campaign

## 4. Recommended trajectory workflow

### 4.1 Bulk reference

Purpose:

- obtain density-consistent composition packing
- provide reference environment for gap setup if needed

Recommended:

- bulk NPT pre-equilibration before interface assembly

### 4.2 Interface stages

Recommended three-stage interface workflow:

1. neutral relaxation
2. constant-potential activation
3. constant-potential production

Recommended default lengths:

- neutral relaxation: `200000 steps`
- cp activation: `100000 steps`
- production: `500000 steps`
- report interval: `1000 steps`

At `2 fs` timestep this corresponds to:

- neutral relaxation: `0.4 ns`
- activation: `0.2 ns`
- production: `1.0 ns`

This is a reasonable first-pass screening length.

For higher-confidence finalists, extend production to `>= 2-5 ns`.

## 5. Standard analysis definitions

These must also remain fixed across formulations.

### 5.1 Interface region

- electrode margin: `1.0 A`
- interface width: `5.0 A`

### 5.2 Bulk region

- use the region between the two interface slabs
- for default geometry this is the middle electrolyte slab after excluding both interfacial windows

### 5.3 Li+ coordination cutoff

- `2.8 A`

### 5.4 Reporting units

Use the same output units everywhere:

- number density: `1 / nm^3`
- surface adsorption: `1 / nm^2`
- capacitance: `uF / cm^2`
- residence time: `ps`

## 6. Standard metadata to record for every run

Every descriptor row should preserve the following metadata:

- `formulation_id`
- `voltage_v`
- `temperature_k`
- `box_x_angstrom`
- `box_y_angstrom`
- `box_z_angstrom`
- `timestep_fs`
- `report_interval`
- `production_steps`
- `gaussian_width_nm`
- `thomas_fermi_scale_invnm`
- `charge_constraint_target_e`
- `electrode_margin_angstrom`
- `interface_width_angstrom`
- `li_coordination_cutoff_angstrom`

Without this metadata, cross-run comparison becomes fragile.

## 7. What should be fixed vs variable

### Fixed across the campaign

- electrode geometry
- OpenMM nonbonded settings
- thermostat settings
- constant-potential solver settings
- analysis window definitions
- voltage grid

### Variable across formulations

- electrolyte composition
- species identities
- molecule counts needed to realize the formulation

### Optional second-phase variables

Only after the first campaign is complete:

- extended production length
- finer voltage grid
- different electrode materials

## 8. Recommended first-pass descriptor campaign

For each formulation:

1. build one standard interface cell
2. run the same multi-voltage schedule
3. extract the same descriptor family at every voltage
4. aggregate to one wide-format row

Recommended first-pass descriptor family:

- additive enrichment factor
- additive surface contact probability
- additive first-layer residence time
- Li-shell replacement fraction
- interfacial capacitance
- interfacial anion ratio
- bulk anion ratio
- delta anion ratio
- anion enrichment factor
- voltage-response features derived from the above

## 9. Common failure modes to avoid

- changing box size between formulations without documenting it
- changing production length unevenly across the dataset
- changing voltage grid between formulations
- using NPT for slab production runs
- redefining interface windows from one system to another
- mixing descriptors extracted with different cutoffs or units

## 10. Minimum recommendation for publication-quality comparability

If the campaign is large, the minimum acceptable standardization is:

- one fixed electrode model
- one fixed geometry family
- one fixed voltage schedule
- one fixed production length
- one fixed analysis protocol
- one fixed descriptor schema

This is the minimum needed for a defensible descriptor-learning dataset.
