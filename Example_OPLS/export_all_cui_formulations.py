#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Export all Cui 2023 electrolyte formulations without MD feasibility filtering.")
    parser.add_argument("--master", required=True, help="literature_labels_master.csv")
    parser.add_argument("--out-unique", required=True, help="Unique formulation table")
    parser.add_argument("--out-all-rows", required=True, help="All rows table")
    parser.add_argument("--out-summary", required=True, help="Summary text")
    args = parser.parse_args()

    master = pd.read_csv(Path(args.master).expanduser().resolve())

    all_rows = master.copy()

    agg = (
        master.groupby("formulation_id", dropna=False)
        .agg(
            raw_electrolyte_text=("raw_electrolyte_text", "first"),
            cell_type=("cell_type", "first"),
            salt_name=("salt_name", "first"),
            salt_name_2=("salt_name_2", "first"),
            salt_name_3=("salt_name_3", "first"),
            salt_concentration_m=("salt_concentration_m", "first"),
            salt_concentration_2_m=("salt_concentration_2_m", "first"),
            salt_concentration_3_m=("salt_concentration_3_m", "first"),
            anion_name=("anion_name", "first"),
            solvent_1=("solvent_1", "first"),
            solvent_1_fraction=("solvent_1_fraction", "first"),
            solvent_2=("solvent_2", "first"),
            solvent_2_fraction=("solvent_2_fraction", "first"),
            solvent_3=("solvent_3", "first"),
            solvent_3_fraction=("solvent_3_fraction", "first"),
            additive_1=("additive_1", "first"),
            additive_1_fraction=("additive_1_fraction", "first"),
            additive_2=("additive_2", "first"),
            additive_2_fraction=("additive_2_fraction", "first"),
            n_conditions=("literature_id", "count"),
            ce_mean_min=("ce_mean", "min"),
            ce_mean_max=("ce_mean", "max"),
            current_density_min_ma_cm2=("current_density_ma_cm2", "min"),
            current_density_max_ma_cm2=("current_density_ma_cm2", "max"),
            areal_capacity_min_mah_cm2=("areal_capacity_mah_cm2", "min"),
            areal_capacity_max_mah_cm2=("areal_capacity_mah_cm2", "max"),
            notes=("notes", lambda s: " | ".join(sorted({str(v) for v in s if pd.notna(v) and str(v).strip()}))),
            references=("raw_source_location", lambda s: " || ".join(sorted({str(v) for v in s if pd.notna(v) and str(v).strip()}))),
        )
        .reset_index()
        .sort_values(["anion_name", "ce_mean_min", "formulation_id"], na_position="last")
    )

    out_unique = Path(args.out_unique).expanduser().resolve()
    out_all = Path(args.out_all_rows).expanduser().resolve()
    out_summary = Path(args.out_summary).expanduser().resolve()
    out_unique.parent.mkdir(parents=True, exist_ok=True)
    out_all.parent.mkdir(parents=True, exist_ok=True)
    out_summary.parent.mkdir(parents=True, exist_ok=True)

    agg.to_csv(out_unique, index=False)
    all_rows.to_csv(out_all, index=False)

    lines = []
    lines.append("All Cui formulations export summary")
    lines.append(f"all_rows: {len(all_rows)}")
    lines.append(f"unique_formulations: {agg['formulation_id'].nunique()}")
    lines.append("anion_counts:")
    for name, count in agg["anion_name"].fillna("NA").value_counts().items():
        lines.append(f"  {name}: {int(count)}")
    lines.append("top_duplicate_formulations:")
    dup = agg.sort_values("n_conditions", ascending=False).head(15)
    for row in dup.itertuples():
        lines.append(f"  {row.formulation_id}: n_conditions={row.n_conditions}, ce_range=[{row.ce_mean_min:.2f}, {row.ce_mean_max:.2f}]")
    out_summary.write_text("\n".join(lines) + "\n")

    print(f"Wrote {out_unique}")
    print(f"Wrote {out_all}")
    print(f"Wrote {out_summary}")


if __name__ == "__main__":
    main()
