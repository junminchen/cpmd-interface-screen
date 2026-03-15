#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


def slugify(text: str) -> str:
    text = str(text).strip().lower()
    text = re.sub(r"[^0-9a-zA-Z]+", "_", text).strip("_")
    return text or "row"


def clean_dataset_sheet(path: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    raw = pd.read_excel(path, sheet_name="Dataset")
    units_row = raw.iloc[0].fillna("")
    data = raw.iloc[1:].copy().reset_index(drop=True)

    data.columns = [
        "row_id",
        "electrolyte",
        "solvent_1_volume_fraction",
        "solvent_2_volume_fraction",
        "solvent_3_volume_fraction",
        "solvent_1_molarity",
        "solvent_2_molarity",
        "solvent_3_molarity",
        "salt_1_molarity",
        "salt_2_molarity",
        "salt_3_molarity",
        "fc",
        "oc",
        "fo",
        "inor",
        "f_total",
        "f_solvent",
        "f_anion",
        "o_total",
        "o_solvent",
        "o_anion",
        "c_total",
        "c_solvent",
        "c_anion",
        "ce_percent",
        "lce",
        "method",
        "current_ma_cm2",
        "capacity_mah_cm2",
        "cycle",
        "reference_id",
    ]

    unit_map = {col: str(units_row.iloc[idx]).strip() for idx, col in enumerate(data.columns)}

    data["literature_id"] = [f"cui2023_{idx:03d}" for idx in range(1, len(data) + 1)]
    data["reference_short"] = "Kim et al. PNAS 2023"
    data["doi"] = "10.1073/pnas.2214357120"
    data["year"] = 2023
    data["source_dataset"] = "cui_pnas_2023_dataset_s01"
    data["formulation_id"] = data["electrolyte"].map(slugify)

    numeric_cols = [
        "solvent_1_volume_fraction",
        "solvent_2_volume_fraction",
        "solvent_3_volume_fraction",
        "solvent_1_molarity",
        "solvent_2_molarity",
        "solvent_3_molarity",
        "salt_1_molarity",
        "salt_2_molarity",
        "salt_3_molarity",
        "fc",
        "oc",
        "fo",
        "inor",
        "f_total",
        "f_solvent",
        "f_anion",
        "o_total",
        "o_solvent",
        "o_anion",
        "c_total",
        "c_solvent",
        "c_anion",
        "ce_percent",
        "lce",
        "current_ma_cm2",
        "capacity_mah_cm2",
        "cycle",
        "reference_id",
    ]
    for col in numeric_cols:
        data[col] = pd.to_numeric(data[col], errors="coerce")

    desired_order = [
        "literature_id",
        "reference_short",
        "doi",
        "year",
        "source_dataset",
        "formulation_id",
        "electrolyte",
        "method",
        "current_ma_cm2",
        "capacity_mah_cm2",
        "cycle",
        "ce_percent",
        "lce",
        "reference_id",
        "solvent_1_volume_fraction",
        "solvent_2_volume_fraction",
        "solvent_3_volume_fraction",
        "solvent_1_molarity",
        "solvent_2_molarity",
        "solvent_3_molarity",
        "salt_1_molarity",
        "salt_2_molarity",
        "salt_3_molarity",
        "fc",
        "oc",
        "fo",
        "inor",
        "f_total",
        "f_solvent",
        "f_anion",
        "o_total",
        "o_solvent",
        "o_anion",
        "c_total",
        "c_solvent",
        "c_anion",
    ]
    data = data[desired_order]
    return data, unit_map


def clean_molecular_sheet(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name="Molecular Database")
    data = raw.iloc[2:].copy().reset_index(drop=True)
    data.columns = [
        "row_id",
        "molecule",
        "count_c",
        "count_h",
        "count_o",
        "count_f",
        "count_n",
        "count_s",
        "count_p",
        "unused",
        "other_elements",
        "molecular_weight",
        "density_g_ml",
        "molarity_pure",
        "smiles",
    ]
    data = data.drop(columns=["row_id", "unused"])
    for col in [
        "count_c",
        "count_h",
        "count_o",
        "count_f",
        "count_n",
        "count_s",
        "count_p",
        "molecular_weight",
        "density_g_ml",
        "molarity_pure",
    ]:
        data[col] = pd.to_numeric(data[col], errors="coerce")
    return data


def clean_reference_sheet(path: Path) -> pd.DataFrame:
    refs = pd.read_excel(path, sheet_name="References")
    refs = refs.rename(columns={"References": "reference_id", "Unnamed: 2": "reference_text"})
    refs = refs[["reference_id", "reference_text"]].copy()
    refs["reference_id"] = pd.to_numeric(refs["reference_id"], errors="coerce")
    refs = refs.dropna(subset=["reference_id"]).reset_index(drop=True)
    refs["reference_id"] = refs["reference_id"].astype(int)
    return refs


def write_summary(path: Path, labels: pd.DataFrame, references: pd.DataFrame, unit_map: dict[str, str]) -> None:
    lines = []
    lines.append("Cui PNAS 2023 CE dataset import summary")
    lines.append(f"rows: {len(labels)}")
    lines.append(f"unique_formulations: {labels['formulation_id'].nunique()}")
    lines.append(f"references: {len(references)}")
    lines.append(f"ce_min: {labels['ce_percent'].min():.3f}")
    lines.append(f"ce_max: {labels['ce_percent'].max():.3f}")
    lines.append(f"current_min_ma_cm2: {labels['current_ma_cm2'].min():.3f}")
    lines.append(f"current_max_ma_cm2: {labels['current_ma_cm2'].max():.3f}")
    lines.append("units_from_source:")
    for key, value in unit_map.items():
        if value:
            lines.append(f"  {key}: {value}")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Ingest Cui PNAS 2023 Dataset S01 into local CSV tables.")
    parser.add_argument("--input-xlsx", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    input_path = Path(args.input_xlsx).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    labels, unit_map = clean_dataset_sheet(input_path)
    molecules = clean_molecular_sheet(input_path)
    references = clean_reference_sheet(input_path)
    labels = labels.merge(references, on="reference_id", how="left")

    labels_path = out_dir / "cui_pnas_2023_ce_labels.csv"
    molecules_path = out_dir / "cui_pnas_2023_molecular_database.csv"
    references_path = out_dir / "cui_pnas_2023_references.csv"
    summary_path = out_dir / "cui_pnas_2023_import_summary.txt"

    labels.to_csv(labels_path, index=False)
    molecules.to_csv(molecules_path, index=False)
    references.to_csv(references_path, index=False)
    write_summary(summary_path, labels, references, unit_map)

    print(f"Wrote {labels_path}")
    print(f"Wrote {molecules_path}")
    print(f"Wrote {references_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
