# Additive Descriptor Workflow

This workflow turns constant-potential MD outputs into sortable scalar descriptors for additive screening.

## Target descriptors

- `*_interface_enrichment_factor`
  - Definition: `rho_interface / rho_bulk`
  - Meaning: whether the additive enriches in the electrode-adjacent slab.
- `*_surface_contact_probability`
  - Definition: fraction of additive molecule-frames located in the first adsorption layer.
  - Meaning: how often the additive contacts the electrode.
- `*_first_layer_residence_time_ps`
  - Definition: mean continuous dwell time in the first adsorption layer.
  - Meaning: kinetic stability of additive adsorption.
- `interfacial_capacitance_uF_per_cm2`
  - Definition: `|Q| / (A * DeltaV)` from the constant-potential charge log.
  - Meaning: electrical response of the interface.
- `*_li_shell_replacement_fraction_*`
  - Definition: `CN_additive / (CN_additive + CN_reference)`
  - Meaning: how strongly the additive replaces reference solvent donors in the Li+ first shell.
- planned next extension: `anion-ratio-like interfacial descriptors across voltage`
  - Definition: voltage-resolved anion participation in interfacial Li+ first-shell coordination
  - Meaning: higher-order proxy for CE-relevant interface reconstruction under bias

`CN` is counted from donor-atom contacts within `li_coordination_cutoff_angstrom`.

## Run

Example for the local inert-electrode demo:

```bash
python Example_OPLS/analyze_additive_descriptors.py \
  --top Example_OPLS/LiPF6_EC_DMC_inert_electrode/start_with_electrodes_mc.pdb \
  --traj Example_OPLS/LiPF6_EC_DMC_inert_electrode/traj_inert.dcd \
  --config Example_OPLS/LiPF6_EC_DMC_inert_electrode/config.json \
  --descriptor-config Example_OPLS/descriptor_config_demo_ec.json \
  --charge-log Example_OPLS/LiPF6_EC_DMC_inert_electrode/electrode_charges.log \
  --out-csv Example_OPLS/LiPF6_EC_DMC_inert_electrode/results/additive_descriptors_demo.csv \
  --out-json Example_OPLS/LiPF6_EC_DMC_inert_electrode/results/additive_descriptors_demo.json
```

## How to adapt for screening

1. Copy `Example_OPLS/descriptor_config_template.json`.
2. Replace `additive_species` with the residue name(s) used by your additive system.
3. Keep `reference_species` as the solvents or other competitors you want the additive to replace around Li+.
4. Run one simulation per formulation and write one descriptor row per run.
5. Aggregate all descriptor rows into one training table.

If anions should count as shell competitors, include them in `reference_species`.

## Aggregate many formulations

```bash
python Example_OPLS/aggregate_descriptor_dataset.py \
  --root Example_OPLS \
  --pattern "additive_descriptors*.json" \
  --out-csv Example_OPLS/results/all_md_descriptors.csv \
  --out-summary Example_OPLS/results/all_md_descriptors_summary.txt
```

If you already have experimental labels, add `--labels your_labels.csv --join-key formulation_id`.

Recommended table structure for later supervised learning:

- `formulation_id`
- MD descriptor columns from `analyze_additive_descriptors.py`
- future voltage-resolved interfacial anion proxy columns
- experimental targets such as `coulombic_efficiency_percent`, `capacity_retention_percent`, `interfacial_resistance_ohm_cm2`

## Next implementation target

The next descriptor expansion should add:

1. explicit anion species tracking
2. interfacial `anion_ratio` and `delta_anion_ratio`
3. the same quantities at multiple voltages for the same formulation
4. derived voltage-response features suitable for supervised learning
