# Descriptor-to-Target Development Log

## Current status

The current codebase implements descriptor construction and dataset aggregation.

Campaign-standard simulation details are now documented in:

- `Example_OPLS/STANDARD_CONSTANT_POTENTIAL_MD_PROTOCOL.md`
- `Example_OPLS/standard_cpmd_campaign_template.json`

Implemented components:

- `analyze_additive_descriptors.py`
  - converts one constant-potential MD trajectory into one descriptor row
- `aggregate_descriptor_dataset.py`
  - merges many descriptor rows into one table for later supervised learning
- `descriptor_config_template.json`
  - defines additive species, reference species, and coordination-shell settings

Current descriptor families:

- interfacial enrichment
- first-layer contact probability and residence time
- Li+ first-shell replacement fractions
- interfacial capacitance from electrode charge logs

Current interpretation:

- these descriptors are mechanistic intermediate variables
- they are intended to bridge formulation chemistry and electrochemical targets
- they are not themselves a predictive model

## What is not implemented yet

No predictive model is implemented at the moment.

The repository does not currently contain:

- a trained regression model for CE or other electrochemical targets
- a neural network
- an attention-based architecture
- a graph neural network
- an end-to-end model from raw trajectory to target

In the current workflow, MD is used only to produce physically interpretable intermediate variables.

## Intended learning pipeline

The intended training path is:

1. run one MD simulation per formulation
2. extract one descriptor row per formulation
3. aggregate all rows into one descriptor table
4. join with experimental targets by `formulation_id`
5. train supervised models on the descriptor table

The learned mapping will therefore be:

`formulation -> MD trajectory -> descriptors -> supervised model -> target`

not:

`formulation -> raw MD frames -> neural network -> target`

## New design direction: voltage-resolved interfacial anion proxies

The next descriptor expansion should not only measure additive behavior.

It should explicitly measure voltage-dependent interfacial anion behavior, because this is closer to CE-relevant interface physics than a bulk-only solvation proxy.

Target concept:

- build a descriptor family that plays the role of a higher-order proxy relative to the `anion ratio` idea used in Bamboo-Mixer
- but compute it under constant-potential interface MD, not bulk mixture MD

Recommended new descriptor family:

- `interfacial_anion_ratio_at_V`
  - fraction of interfacial Li+ first-shell coordination contributed by anions
- `bulk_anion_ratio_at_V`
  - bulk reference value
- `delta_anion_ratio_at_V`
  - interfacial minus bulk
- `anion_enrichment_factor_at_V`
  - interfacial anion density divided by bulk anion density
- `li_anion_coordination_number_interface_at_V`
- `li_anion_coordination_number_bulk_at_V`
- `anion_residence_time_near_interfacial_li_at_V`
- `interfacial_capacitance_at_V`
- voltage-response descriptors such as:
  - `d_interfacial_anion_ratio_dV`
  - `d_anion_enrichment_dV`
  - `d_capacitance_dV`

This would create a voltage-resolved descriptor tensor that is still small and interpretable, but richer than a single bulk solvation proxy.

## Next development steps

### Phase 1: data plumbing

1. finalize a stable `formulation_id` mapping between MD cases and the experimental label table
2. prepare a clean label table with fields such as:
   - `formulation_id`
   - `coulombic_efficiency_percent`
   - `capacity_retention_percent`
   - `interfacial_resistance_ohm_cm2`
3. run all formulation MD jobs and extract descriptors
4. aggregate all descriptor rows into one training dataset

### Immediate next tasks

1. extend descriptor extraction to support explicit anion species definitions
2. add interface-specific anion-ratio outputs for one voltage
3. extend the workflow to multiple voltages per formulation
4. define a canonical wide-format output schema:
   - one row per formulation
   - voltage-indexed columns for each interfacial anion proxy
5. derive voltage-response summary features from the multi-voltage runs
6. only after that, join with experimental CE labels

### Phase 2: baseline supervised models

Start with simple tabular models:

- ridge / lasso
- random forest
- xgboost or lightgbm if available later

Reason:

- sample count will likely be limited
- descriptor dimension is modest
- interpretability matters

### Phase 3: only if data volume is large enough

Consider neural models only after enough labeled formulations exist.

Possible later options:

- MLP on descriptor vectors
- attention-based tabular models such as TabTransformer / FT-Transformer

These are not justified yet unless the labeled dataset grows beyond the small-data regime.

## Practical recommendation right now

The correct immediate focus is not architecture design.

The correct immediate focus is:

- descriptor reproducibility across all formulations
- voltage consistency across the same formulation
- robust interfacial anion proxy definitions
- protocol consistency across the screening campaign
- label alignment quality
- train/validation split design
- baseline model sanity checks

Only after those are stable does it make sense to test neural-network or attention-based models.
