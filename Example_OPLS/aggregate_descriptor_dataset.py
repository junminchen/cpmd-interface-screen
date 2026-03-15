#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def discover_files(root: Path, pattern: str) -> list[Path]:
    return sorted(path.resolve() for path in root.rglob(pattern) if path.is_file())


def load_descriptor_row(path: Path) -> dict:
    if path.suffix.lower() == ".json":
        record = json.loads(path.read_text())
    elif path.suffix.lower() == ".csv":
        frame = pd.read_csv(path)
        if frame.empty:
            raise RuntimeError(f"Descriptor CSV is empty: {path}")
        record = frame.iloc[0].to_dict()
    else:
        raise RuntimeError(f"Unsupported descriptor file type: {path}")
    record["descriptor_source_path"] = str(path)
    record["descriptor_source_dir"] = str(path.parent)
    return record


def normalize_record(record: dict, fallback_id: str) -> dict:
    out = dict(record)
    if "system_name" not in out or not str(out["system_name"]).strip():
        out["system_name"] = fallback_id
    if "formulation_id" not in out or not str(out["formulation_id"]).strip():
        out["formulation_id"] = out["system_name"]
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Aggregate many MD descriptor rows into one training table."
    )
    parser.add_argument("--root", required=True, help="Root directory to scan")
    parser.add_argument(
        "--pattern",
        default="additive_descriptors*.json",
        help="Glob used with rglob under --root",
    )
    parser.add_argument(
        "--labels",
        default=None,
        help="Optional CSV of experimental targets to join",
    )
    parser.add_argument(
        "--join-key",
        default="formulation_id",
        help="Join key used for --labels, e.g. formulation_id or system_name",
    )
    parser.add_argument("--out-csv", required=True, help="Aggregated dataset CSV")
    parser.add_argument("--out-summary", required=True, help="Text summary")
    args = parser.parse_args()

    root = Path(args.root).expanduser().resolve()
    files = discover_files(root, args.pattern)
    if not files:
        raise RuntimeError(f"No descriptor files found under {root} with pattern {args.pattern}")

    rows = []
    for idx, path in enumerate(files, start=1):
        fallback_id = f"row_{idx:04d}"
        rows.append(normalize_record(load_descriptor_row(path), fallback_id))

    dataset = pd.DataFrame(rows)

    duplicated = dataset[args.join_key].duplicated().sum() if args.join_key in dataset.columns else 0
    if duplicated:
        raise RuntimeError(f"Found {duplicated} duplicated values in join key '{args.join_key}'")

    label_columns = []
    if args.labels:
        labels_path = Path(args.labels).expanduser().resolve()
        labels = pd.read_csv(labels_path)
        if args.join_key not in labels.columns:
            raise RuntimeError(f"Join key '{args.join_key}' not found in labels file {labels_path}")
        label_columns = [col for col in labels.columns if col != args.join_key]
        dataset = dataset.merge(labels, on=args.join_key, how="left", validate="one_to_one")

    out_csv = Path(args.out_csv).expanduser().resolve()
    out_summary = Path(args.out_summary).expanduser().resolve()
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_summary.parent.mkdir(parents=True, exist_ok=True)
    dataset.to_csv(out_csv, index=False)

    numeric_cols = dataset.select_dtypes(include="number").columns.tolist()
    lines = []
    lines.append("Aggregated descriptor dataset")
    lines.append(f"root: {root}")
    lines.append(f"pattern: {args.pattern}")
    lines.append(f"rows: {len(dataset)}")
    lines.append(f"columns: {len(dataset.columns)}")
    lines.append(f"join_key: {args.join_key}")
    lines.append(f"numeric_columns: {len(numeric_cols)}")
    if args.labels:
        lines.append(f"labels_file: {Path(args.labels).expanduser().resolve()}")
        lines.append(f"label_columns_joined: {', '.join(label_columns) if label_columns else '(none)'}")
        missing_counts = dataset[label_columns].isna().sum() if label_columns else pd.Series(dtype=int)
        for name, count in missing_counts.items():
            lines.append(f"missing_{name}: {int(count)}")
    lines.append("descriptor_files:")
    for path in files:
        lines.append(f"  {path}")
    out_summary.write_text("\n".join(lines) + "\n")

    print(f"Rows aggregated: {len(dataset)}")
    print(f"Wrote dataset: {out_csv}")
    print(f"Wrote summary: {out_summary}")


if __name__ == "__main__":
    main()
