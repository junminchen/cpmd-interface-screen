# cpmd-interface-screen

Constant-potential MD interface screening workflow for electrolyte formulations.

This repository combines:

- constant-potential interface MD workflows
- mechanistic descriptor extraction
- literature CE dataset ingestion and curation
- pilot formulation selection for MD validation

## Main components

- `Example_OPLS/analyze_additive_descriptors.py`
  - extract interfacial descriptor rows from one MD trajectory
- `Example_OPLS/aggregate_descriptor_dataset.py`
  - combine many descriptor rows into one training table
- `Example_OPLS/STANDARD_CONSTANT_POTENTIAL_MD_PROTOCOL.md`
  - standardized simulation protocol for screening campaigns
- `Example_OPLS/standard_cpmd_campaign_template.json`
  - campaign-level parameter template
- `Example_OPLS/run_standard_cpmd_test.sh`
  - one-command smoke test / standard test entrypoint

## Standard test system

The current standard test system is:

- `Example_OPLS/LiPF6_EC_DMC_thick_Li_electrode/three_stage_cp_workflow`

Quick smoke test:

```bash
bash Example_OPLS/run_standard_cpmd_test.sh smoke CPU
```

Standard protocol-length test:

```bash
bash Example_OPLS/run_standard_cpmd_test.sh standard CPU
```

If CUDA is available:

```bash
bash Example_OPLS/run_standard_cpmd_test.sh smoke CUDA
```

## Literature data

Cleaned Cui 2023 PNAS tables are stored under:

- `data/processed/cui_pnas_2023_ce_labels.csv`
- `data/processed/cui_all_unique_formulations.csv`
- `data/processed/literature_labels_master.csv`

## Current project direction

The current modeling direction is:

`formulation -> constant-potential MD -> interfacial descriptors -> supervised learning target`

The next planned descriptor extension is:

- voltage-resolved interfacial anion-ratio-like proxies

## Notes

- This repo currently focuses on descriptor construction and dataset preparation.
- It does not yet contain a final trained CE prediction model.
