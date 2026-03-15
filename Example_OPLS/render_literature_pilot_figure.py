#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd


COLORS = {
    "all": "#b8c4d6",
    "usable": "#4c78a8",
    "top30": "#d95f02",
    "PF6": "#1b9e77",
    "FSI": "#d95f02",
    "TFSI": "#7570b3",
}


def load_tables(master_path: Path, top30_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    master = pd.read_csv(master_path)
    top30 = pd.read_csv(top30_path)
    return master, top30


def render(master: pd.DataFrame, top30: pd.DataFrame, out_path: Path) -> None:
    usable = master[master["usable_for_md"] == "yes"].copy()

    fig = plt.figure(figsize=(14, 10), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.1])

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])

    bins = 18
    ax1.hist(master["ce_mean"], bins=bins, color=COLORS["all"], alpha=0.7, label="All literature")
    ax1.hist(usable["ce_mean"], bins=bins, color=COLORS["usable"], alpha=0.75, label="Usable for MD")
    ax1.hist(top30["ce_mean"], bins=bins, color=COLORS["top30"], alpha=0.85, label="Top30 shortlist")
    ax1.set_xlabel("Coulombic efficiency (%)")
    ax1.set_ylabel("Count")
    ax1.set_title("CE Distribution")
    ax1.legend(frameon=False)
    ax1.grid(alpha=0.2)

    ax2.scatter(
        master["current_density_ma_cm2"],
        master["areal_capacity_mah_cm2"],
        s=34,
        alpha=0.35,
        color=COLORS["all"],
        label="All literature",
    )
    ax2.scatter(
        usable["current_density_ma_cm2"],
        usable["areal_capacity_mah_cm2"],
        s=42,
        alpha=0.75,
        color=COLORS["usable"],
        label="Usable for MD",
    )
    ax2.scatter(
        top30["current_density_ma_cm2"],
        top30["areal_capacity_mah_cm2"],
        s=52,
        alpha=0.95,
        color=COLORS["top30"],
        edgecolor="black",
        linewidth=0.4,
        label="Top30 shortlist",
    )
    ax2.set_xlabel("Current density (mA cm$^{-2}$)")
    ax2.set_ylabel("Areal capacity (mAh cm$^{-2}$)")
    ax2.set_title("Protocol Window")
    ax2.legend(frameon=False, loc="lower right")
    ax2.grid(alpha=0.2)

    plot_df = top30.sort_values("ce_mean").reset_index(drop=True).copy()
    plot_df["rank"] = range(1, len(plot_df) + 1)
    for anion, group in plot_df.groupby("anion_name"):
        color = COLORS.get(anion, "#666666")
        ax3.scatter(
            group["rank"],
            group["ce_mean"],
            s=78,
            color=color,
            alpha=0.9,
            label=anion,
        )
        for row in group.itertuples():
            short_name = str(row.formulation_id)[:28]
            ax3.text(
                row.rank + 0.2,
                row.ce_mean,
                short_name,
                fontsize=7,
                alpha=0.8,
                va="center",
            )

    ax3.set_xlabel("Shortlist rank (sorted by CE)")
    ax3.set_ylabel("Coulombic efficiency (%)")
    ax3.set_title("Top30 Pilot MD Shortlist")
    ax3.grid(alpha=0.2)
    ax3.legend(title="Anion", frameon=False, ncol=3, loc="upper left")

    fig.suptitle(
        "Cui 2023 Literature Dataset: All Rows vs MD-Usable Subset vs Pilot Shortlist",
        fontsize=15,
    )
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Render literature/pilot dataset summary figure.")
    parser.add_argument("--master", required=True)
    parser.add_argument("--top30", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    master, top30 = load_tables(Path(args.master), Path(args.top30))
    out_path = Path(args.out).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    render(master, top30, out_path)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
