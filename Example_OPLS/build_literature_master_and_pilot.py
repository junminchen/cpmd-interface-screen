#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import pandas as pd


SALT_ALIASES = {
    "LIPF6": ("LiPF6", "PF6"),
    "LIBF4": ("LiBF4", "BF4"),
    "LITFSI": ("LiTFSI", "TFSI"),
    "LIFSI": ("LiFSI", "FSI"),
    "LINO3": ("LiNO3", "NO3"),
    "LIBOB": ("LiBOB", "BOB"),
    "LIDFOB": ("LiDFOB", "DFOB"),
    "LICLO4": ("LiClO4", "ClO4"),
    "LIBETI": ("LiBETI", "BETI"),
    "LIASI6": ("LiAsF6", "AsF6"),
    "LIASF6": ("LiAsF6", "AsF6"),
    "LITFPFB": ("LiTFPFB", "TFPFB"),
    "LIPO2F2": ("LiPO2F2", "PO2F2"),
}

SPECIES_ALIASES = {
    "EC": "EC",
    "DMC": "DMC",
    "DEC": "DEC",
    "EMC": "EMC",
    "PC": "PC",
    "FEC": "FEC",
    "DOL": "DOL",
    "DME": "DME",
    "GBL": "GBL",
    "SL": "SL",
    "PS": "PS",
    "PP": "PP",
    "EP": "EP",
    "DFEC": "DFEC",
    "DFEA": "DFEA",
    "FEMC": "FEMC",
    "HFDEC": "HFDEC",
    "VC": "VC",
    "BTFE": "BTFE",
    "FDMB": "FDMB",
    "D2": "D2",
    "TTE": "TTE",
    "FM": "FM",
    "DMB": "DMB",
    "TMS": "TMS",
    "TFEO": "TFEO",
    "TFEC": "TFEC",
    "TEP": "TEP",
    "DEE": "DEE",
}

SUPPORTED_SPECIES = {
    "BF4",
    "BOB",
    "DEC",
    "DFEA",
    "DFEC",
    "DFOB",
    "DFP",
    "DMC",
    "DME",
    "DOL",
    "EC",
    "EMC",
    "EP",
    "FAN",
    "FEC",
    "FEMC",
    "FSI",
    "G2",
    "G3",
    "GBL",
    "HFDEC",
    "NO3",
    "PC",
    "PEO",
    "PF6",
    "PP",
    "PS",
    "SL",
    "TFSI",
}


def tokenize_formula(text: str) -> list[str]:
    upper = str(text).upper()
    upper = upper.replace("(", " ").replace(")", " ").replace("/", " ").replace("-", " ")
    return re.findall(r"[A-Z0-9]+", upper)


def parse_electrolyte(text: str) -> dict:
    tokens = tokenize_formula(text)

    salts = []
    anions = []
    for token, (salt_name, anion_name) in SALT_ALIASES.items():
        if token in tokens or token in str(text).upper():
            if salt_name not in salts:
                salts.append(salt_name)
            if anion_name not in anions:
                anions.append(anion_name)

    species = []
    for token, canonical in SPECIES_ALIASES.items():
        if token in tokens:
            species.append(canonical)

    # Keep order of first appearance in the original string.
    ordered_species = []
    upper = str(text).upper()
    for token, canonical in sorted(SPECIES_ALIASES.items(), key=lambda kv: upper.find(kv[0]) if kv[0] in upper else 10**9):
        if token in upper and canonical not in ordered_species:
            ordered_species.append(canonical)

    solvents = ordered_species[:3]
    additives = ordered_species[3:5]

    unknown_tokens = []
    ignored = {
        "M",
        "WT",
        "VOL",
        "V",
        "W",
    }
    known_tokens = set(SALT_ALIASES) | set(SPECIES_ALIASES) | ignored
    for token in tokens:
        if token in known_tokens:
            continue
        if re.fullmatch(r"\d+(\.\d+)?", token):
            continue
        if re.fullmatch(r"\d+[A-Z]*", token):
            continue
        unknown_tokens.append(token)

    return {
        "salt_names": salts,
        "anion_names": anions,
        "solvents": solvents,
        "additives": additives,
        "unknown_tokens": sorted(set(unknown_tokens)),
    }


def normalize_fraction(value):
    if pd.isna(value):
        return pd.NA
    return float(value)


def build_master(cui: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in cui.itertuples(index=False):
        parsed = parse_electrolyte(row.electrolyte)
        supported_species = set()
        unsupported_species = set()
        for salt in parsed["salt_names"]:
            anion = SALT_ALIASES[salt.upper().replace("LI", "LI")][1] if False else None
        for item in parsed["solvents"] + parsed["additives"] + parsed["anion_names"]:
            if item in SUPPORTED_SPECIES:
                supported_species.add(item)
            else:
                unsupported_species.add(item)

        usable_for_md = (
            len(parsed["salt_names"]) >= 1
            and len(parsed["unknown_tokens"]) == 0
            and len(unsupported_species) == 0
            and len(parsed["solvents"]) >= 1
        )

        md_risk_notes = []
        if len(parsed["salt_names"]) == 0:
            md_risk_notes.append("salt_not_parsed")
        if len(parsed["unknown_tokens"]) > 0:
            md_risk_notes.append("unknown_tokens:" + "|".join(parsed["unknown_tokens"]))
        if len(unsupported_species) > 0:
            md_risk_notes.append("unsupported_species:" + "|".join(sorted(unsupported_species)))

        solvents = parsed["solvents"] + [pd.NA] * max(0, 3 - len(parsed["solvents"]))
        additives = parsed["additives"] + [pd.NA] * max(0, 2 - len(parsed["additives"]))
        salts = parsed["salt_names"] + [pd.NA] * max(0, 3 - len(parsed["salt_names"]))

        rows.append(
            {
                "literature_id": row.literature_id,
                "reference_short": row.reference_short,
                "doi": row.doi,
                "year": row.year,
                "source_dataset": row.source_dataset,
                "formulation_id": row.formulation_id,
                "cell_type": "Li|Cu",
                "salt_name": salts[0],
                "salt_name_2": salts[1],
                "salt_name_3": salts[2],
                "salt_concentration_m": row.salt_1_molarity,
                "salt_concentration_2_m": row.salt_2_molarity,
                "salt_concentration_3_m": row.salt_3_molarity,
                "anion_name": parsed["anion_names"][0] if parsed["anion_names"] else pd.NA,
                "solvent_1": solvents[0],
                "solvent_1_fraction": normalize_fraction(row.solvent_1_volume_fraction),
                "solvent_2": solvents[1],
                "solvent_2_fraction": normalize_fraction(row.solvent_2_volume_fraction),
                "solvent_3": solvents[2],
                "solvent_3_fraction": normalize_fraction(row.solvent_3_volume_fraction),
                "additive_1": additives[0],
                "additive_1_fraction": pd.NA,
                "additive_2": additives[1],
                "additive_2_fraction": pd.NA,
                "temperature_c": pd.NA,
                "current_density_ma_cm2": row.current_ma_cm2,
                "areal_capacity_mah_cm2": row.capacity_mah_cm2,
                "ce_mean": row.ce_percent,
                "ce_std": pd.NA,
                "ce_metric_definition": row.method,
                "cycle_count": row.cycle,
                "interfacial_resistance_ohm_cm2": pd.NA,
                "conductivity_ms_cm": pd.NA,
                "viscosity_cp": pd.NA,
                "transference_number": pd.NA,
                "substrate_type": pd.NA,
                "separator_type": pd.NA,
                "usable_for_md": "yes" if usable_for_md else "no",
                "pilot_priority": "pending",
                "raw_source_location": row.reference_text,
                "notes": ";".join(md_risk_notes) if md_risk_notes else "",
                "raw_electrolyte_text": row.electrolyte,
                "raw_reference_id": row.reference_id,
                "fc": row.fc,
                "oc": row.oc,
                "fo": row.fo,
                "inor": row.inor,
                "f_total": row.f_total,
                "f_solvent": row.f_solvent,
                "f_anion": row.f_anion,
                "o_total": row.o_total,
                "o_solvent": row.o_solvent,
                "o_anion": row.o_anion,
                "c_total": row.c_total,
                "c_solvent": row.c_solvent,
                "c_anion": row.c_anion,
            }
        )
    return pd.DataFrame(rows)


def assign_pilot_priority(df: pd.DataFrame) -> pd.DataFrame:
    pilot = df.copy()
    in_window = (
        pilot["current_density_ma_cm2"].between(0.5, 1.0, inclusive="both")
        & pilot["areal_capacity_mah_cm2"].between(0.5, 1.0, inclusive="both")
    )
    pilot["in_pilot_window"] = in_window

    unique_formulations = (
        pilot[pilot["usable_for_md"] == "yes"]
        .sort_values(["formulation_id", "ce_mean"], ascending=[True, True])
        .drop_duplicates("formulation_id", keep="first")
        .copy()
    )
    unique_formulations = unique_formulations[unique_formulations["in_pilot_window"]].copy()

    low_cut = unique_formulations["ce_mean"].quantile(1 / 3) if len(unique_formulations) else math.nan
    high_cut = unique_formulations["ce_mean"].quantile(2 / 3) if len(unique_formulations) else math.nan

    def ce_bin(value: float) -> str:
        if pd.isna(value):
            return "unknown"
        if value <= low_cut:
            return "low"
        if value >= high_cut:
            return "high"
        return "mid"

    unique_formulations["ce_bin"] = unique_formulations["ce_mean"].map(ce_bin)
    unique_formulations["pilot_priority"] = "medium"
    unique_formulations.loc[unique_formulations["ce_bin"].isin(["low", "high"]), "pilot_priority"] = "high"
    unique_formulations.loc[
        unique_formulations["notes"].str.contains("unsupported_species|unknown_tokens", na=False),
        "pilot_priority",
    ] = "low"

    pilot_map = unique_formulations.set_index("literature_id")["pilot_priority"].to_dict()
    ce_bin_map = unique_formulations.set_index("literature_id")["ce_bin"].to_dict()

    pilot["pilot_priority"] = pilot["literature_id"].map(pilot_map).fillna("low")
    pilot["ce_bin"] = pilot["literature_id"].map(ce_bin_map).fillna("out_of_window")
    return pilot, unique_formulations.sort_values(["pilot_priority", "ce_mean"], ascending=[True, True])


def pick_top30(shortlist: pd.DataFrame) -> pd.DataFrame:
    frames = []
    targets = [("low", 10), ("mid", 10), ("high", 10)]
    for ce_bin, n_pick in targets:
        block = shortlist[shortlist["ce_bin"] == ce_bin].copy()
        block = block.sort_values(
            ["pilot_priority", "anion_name", "ce_mean"],
            ascending=[True, True, True],
        )
        frames.append(block.head(n_pick))
    out = pd.concat(frames, ignore_index=True)
    return out.drop_duplicates("formulation_id").reset_index(drop=True)


def write_summary(path: Path, master: pd.DataFrame, shortlist: pd.DataFrame, top30: pd.DataFrame) -> None:
    lines = []
    lines.append("Literature master + pilot build summary")
    lines.append(f"master_rows: {len(master)}")
    lines.append(f"master_unique_formulations: {master['formulation_id'].nunique()}")
    lines.append(f"usable_for_md_yes: {(master['usable_for_md'] == 'yes').sum()}")
    lines.append(f"pilot_shortlist_rows: {len(shortlist)}")
    lines.append(f"pilot_top30_rows: {len(top30)}")
    lines.append(f"pilot_priority_high: {(shortlist['pilot_priority'] == 'high').sum()}")
    lines.append(f"pilot_priority_medium: {(shortlist['pilot_priority'] == 'medium').sum()}")
    lines.append("pilot_ce_bin_counts:")
    for name, count in shortlist["ce_bin"].value_counts().items():
        lines.append(f"  {name}: {int(count)}")
    lines.append("pilot_anion_counts:")
    for name, count in shortlist["anion_name"].fillna("NA").value_counts().items():
        lines.append(f"  {name}: {int(count)}")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build unified literature labels master table and pilot shortlist.")
    parser.add_argument("--cui-labels", required=True)
    parser.add_argument("--out-master", required=True)
    parser.add_argument("--out-pilot", required=True)
    parser.add_argument("--out-pilot-top30", required=True)
    parser.add_argument("--out-summary", required=True)
    args = parser.parse_args()

    cui = pd.read_csv(Path(args.cui_labels).expanduser().resolve())
    master = build_master(cui)
    pilot_labeled, shortlist = assign_pilot_priority(master)
    top30 = pick_top30(shortlist)

    out_master = Path(args.out_master).expanduser().resolve()
    out_pilot = Path(args.out_pilot).expanduser().resolve()
    out_pilot_top30 = Path(args.out_pilot_top30).expanduser().resolve()
    out_summary = Path(args.out_summary).expanduser().resolve()
    out_master.parent.mkdir(parents=True, exist_ok=True)
    out_pilot.parent.mkdir(parents=True, exist_ok=True)
    out_pilot_top30.parent.mkdir(parents=True, exist_ok=True)
    out_summary.parent.mkdir(parents=True, exist_ok=True)

    master.to_csv(out_master, index=False)
    shortlist.to_csv(out_pilot, index=False)
    top30.to_csv(out_pilot_top30, index=False)
    write_summary(out_summary, master, shortlist, top30)

    print(f"Wrote {out_master}")
    print(f"Wrote {out_pilot}")
    print(f"Wrote {out_pilot_top30}")
    print(f"Wrote {out_summary}")


if __name__ == "__main__":
    main()
