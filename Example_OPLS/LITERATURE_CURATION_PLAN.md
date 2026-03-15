# Literature Curation Plan for CE + MD Descriptor Validation

## 1. Purpose

Before committing to a fully standardized high-throughput campaign, use literature data to test whether the planned constant-potential MD descriptors are likely to carry predictive value.

This stage is for feasibility validation, not final model benchmarking.

## 2. Main objective

Build a literature-derived electrochemical label table that is sufficiently clean to support a pilot comparison of:

1. composition-only baselines
2. bulk-solvation or formula-based proxies
3. voltage-resolved interfacial MD descriptors

## 3. Scope of the first literature dataset

Start from a narrow domain:

- lithium metal related electrolytes
- prioritize `Li|Cu` Coulombic efficiency datasets
- room-temperature or near-room-temperature studies
- protocols close to each other in current density and areal capacity

Only after this pilot works should the dataset be expanded.

## 4. Recommended inclusion criteria

### Include first

- `Li|Cu` half-cell CE data
- clearly reported electrolyte composition
- clearly reported salt concentration
- clearly reported solvent/additive identity
- clearly reported current density
- clearly reported areal capacity
- clearly reported temperature

### Prioritize

- datasets with repeated measurements or error bars
- datasets with EIS or impedance information
- datasets with SEI characterization
- datasets with compositions that can be reproduced in MD

### Exclude for the first pass

- missing or ambiguous composition
- missing salt concentration
- missing test current density or areal capacity
- highly unusual protocols that are not comparable to the main set
- studies where the CE definition is unclear

## 5. Recommended pilot homogeneity window

To reduce protocol-induced label noise, the first pilot subset should be filtered to something like:

- `Li|Cu`
- `0.5-1.0 mA cm^-2`
- `0.5-1.0 mAh cm^-2`
- near room temperature
- similar Cu substrate conditions if reported

This is not a permanent restriction; it is for the first test of descriptor value.

## 6. Data schema

The first master table should be a CSV with one row per experimental condition.

Required columns:

- `literature_id`
- `reference_short`
- `doi`
- `year`
- `formulation_id`
- `cell_type`
- `salt_name`
- `salt_concentration_m`
- `solvent_1`
- `solvent_1_fraction`
- `solvent_2`
- `solvent_2_fraction`
- `solvent_3`
- `solvent_3_fraction`
- `additive_1`
- `additive_1_fraction`
- `additive_2`
- `additive_2_fraction`
- `anion_name`
- `temperature_c`
- `current_density_ma_cm2`
- `areal_capacity_mah_cm2`
- `ce_mean`
- `ce_std`
- `ce_metric_definition`
- `cycle_count`

Recommended additional columns:

- `interfacial_resistance_ohm_cm2`
- `conductivity_ms_cm`
- `viscosity_cp`
- `transference_number`
- `substrate_type`
- `separator_type`
- `notes`
- `raw_source_location`

## 7. `formulation_id` rule

`formulation_id` must be stable and reproducible.

Recommended rule:

- encode salt identity
- encode salt concentration
- encode solvent mixture
- encode additive identity and fraction

Example style:

- `LiPF6_1p0M_EC50_DMC50_FEC10`
- `LiFSI_1p2M_DME50_DOL50_LiNO3_2wt`

The exact formatting can be standardized later, but the same formulation must never appear under multiple IDs.

## 8. Curation workflow

### Step 1

Seed the table with well-known structured datasets:

- Cui 2023 PNAS CE dataset
- closely related Li|Cu literature with explicit composition tables

### Step 2

Tag each row with:

- `usable_for_md = yes/no`
- `pilot_priority = high/medium/low`

`usable_for_md = yes` means:

- composition is sufficiently explicit
- force-field coverage is plausible
- protocol is not too exotic

### Step 3

Construct a first pilot subset of around `20-40` formulations with:

- diverse anions
- diverse additives
- a reasonable CE range
- minimal protocol spread

## 9. Pilot MD subset selection rules

Choose formulations that satisfy most of the following:

1. composition can be unambiguously encoded into simulation input
2. additive and solvent species exist in the current force-field library or can be added with manageable effort
3. protocol is near the pilot homogeneity window
4. CE span is broad enough to test ranking power
5. multiple chemical families are represented

The first pilot should intentionally include:

- several high-CE formulations
- several mid-CE formulations
- several poor formulations

Otherwise the descriptor test will be underpowered.

## 10. What to validate in the literature phase

The literature phase should answer:

1. do the planned interfacial descriptors vary meaningfully across formulations
2. do they correlate with CE in the expected direction
3. do they outperform simple composition-only baselines
4. are voltage-resolved interfacial anion proxies better than bulk-only proxies

## 11. Success criteria for moving to HTE

Move to standardized high-throughput experiments only if the literature pilot shows at least one of the following:

- meaningful rank correlation between MD descriptors and CE
- clear improvement over composition-only models
- chemically interpretable trends that align with known interface physics

If none of these hold, revise the descriptor family before scaling up experimentally.

## 12. Immediate next tasks

1. create `literature_labels_master.csv`
2. create a `README` describing each column and its allowed units
3. ingest the Cui dataset first
4. manually add a small number of high-quality Li|Cu studies
5. tag `usable_for_md` and `pilot_priority`
6. freeze the first pilot subset for MD
