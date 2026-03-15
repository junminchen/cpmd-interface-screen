# Cui PNAS 2023 CE Dataset

Files in this folder:

- `cui_pnas_2023_ce_labels.csv`
- `cui_pnas_2023_molecular_database.csv`
- `cui_pnas_2023_references.csv`
- `cui_pnas_2023_import_summary.txt`

## Source

The files were cleaned from:

- `data/raw/cui_pnas_2023_dataset_s01.xlsx`
- article DOI: `10.1073/pnas.2214357120`

## Important usage note

Use `literature_id` as the unique row identifier.

Do not assume `formulation_id` is unique in this dataset.

Reason:

- the same formulation may appear under multiple experimental conditions
- some normalized formulation strings can collapse distinct textual entries

If later you want one row per formulation, you should aggregate explicitly, for example by:

- restricting to a homogeneous protocol subset first
- then averaging `ce_percent` within matched conditions

## Recommended first-pass modeling subset

For the first pilot, filter to a more homogeneous window such as:

- `current_ma_cm2` between `0.5` and `1.0`
- `capacity_mah_cm2` between `0.5` and `1.0`
- non-missing `cycle`
- clearly reproducible electrolyte compositions

This is the safest starting point for comparing:

- composition-only baselines
- bulk proxy baselines
- MD descriptor models
