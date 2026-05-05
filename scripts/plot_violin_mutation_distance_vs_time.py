#!/usr/bin/env python3
"""
2D jitter plot: mutation DNA distance distribution vs time for selected groups.

Uses the same selection criteria as the 3D plot:
- replicate 1, starting WT/E
- replicate 3, starting WT/E
- replicate 2, starting D3
- replicate 3, starting D3

Data source:
- all_summary.xlsx sheet plotted_reads_gt_1000 for metadata + results_folder
- each sample's sample_results.xlsx sheet dna_distance_distribution for counts

Dedup:
- if there are two rows with (timepoint=7, replicate=3, Starting point=D3), keep only the first.

Output:
- PDF jitter (strip) plot where x=time (continuous mapped), y=mutation distance, hue=group.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

TIMEPOINT_MAP = {
    0: 0.0,
    1: 5.6,
    2: 26.1,
    3: 41.5,
    4: 49.5,
    5: 64.7,
    6: 73.6,
    7: 90.0,
}


def _prefer_arial() -> None:
    try:
        import matplotlib.font_manager as _fm

        preferred = ["Arial", "Liberation Sans", "DejaVu Sans"]
        available = {f.name for f in _fm.fontManager.ttflist}
        for name in preferred:
            if name in available:
                plt.rcParams["font.family"] = name
                return
    except Exception:
        return


def load_distance_counts(sample_results_xlsx: Path) -> pd.DataFrame:
    df = pd.read_excel(sample_results_xlsx, sheet_name="dna_distance_distribution", engine="openpyxl")
    if df.empty or "dna_distance" not in df.columns:
        return pd.DataFrame(columns=["dna_distance", "read_count"])
    out = df.copy()
    out["dna_distance"] = pd.to_numeric(out["dna_distance"], errors="coerce")
    if "read_count" in out.columns:
        out["read_count"] = pd.to_numeric(out["read_count"], errors="coerce")
    else:
        # Fall back to fraction if counts are absent.
        if "fraction" in out.columns:
            out["fraction"] = pd.to_numeric(out["fraction"], errors="coerce")
        out["read_count"] = pd.NA
    out = out.dropna(subset=["dna_distance"])
    out = out[out["dna_distance"] >= 0]
    return out[["dna_distance", "read_count"]]


def main() -> None:
    _prefer_arial()

    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--base_dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent,
        help="Project root containing all_summary.xlsx and aa_analysis_results_* folders",
    )
    ap.add_argument(
        "--summary_xlsx",
        type=Path,
        default=None,
        help="Path to all_summary.xlsx (default: <base_dir>/all_summary.xlsx)",
    )
    ap.add_argument("--sheet", default="plotted_reads_gt_1000")
    ap.add_argument(
        "--out_pdf",
        type=Path,
        default=None,
        help="Output PDF path (default: <base_dir>/plots_mutation_distance_violin/mutation_distance_violin.pdf)",
    )
    args = ap.parse_args()

    base = args.base_dir.resolve()
    summary = (args.summary_xlsx or (base / "all_summary.xlsx")).resolve()
    out_pdf = args.out_pdf or (base / "plots_mutation_distance_violin" / "mutation_distance_violin.pdf")
    out_pdf.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_excel(summary, sheet_name=args.sheet, engine="openpyxl")
    required = {"sample_id", "results_folder", "timepoint", "replicate", "Starting point"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {args.sheet}: {sorted(missing)}")

    df = df.copy()
    df["timepoint"] = pd.to_numeric(df["timepoint"], errors="coerce")
    df["replicate"] = pd.to_numeric(df["replicate"], errors="coerce")
    df = df.dropna(subset=["timepoint", "replicate", "Starting point"])
    df["Starting point"] = df["Starting point"].astype(str)

    keep_groups = {
        (1, "WT/E"),
        (3, "WT/E"),
        (2, "D3"),
        (3, "D3"),
    }
    df = df[df.apply(lambda r: (int(r["replicate"]), r["Starting point"]) in keep_groups, axis=1)].copy()

    mask_dup = (df["timepoint"] == 7) & (df["replicate"] == 3) & (df["Starting point"] == "D3")
    dup_idx = list(df.index[mask_dup])
    if len(dup_idx) >= 2:
        df = df.drop(index=dup_idx[1:])

    df["time"] = df["timepoint"].astype(int).map(lambda t: TIMEPOINT_MAP.get(t, float(t)))
    df["group"] = df.apply(lambda r: f"rep{int(r['replicate'])}_{r['Starting point']}", axis=1)

    rows = []
    for _, r in df.iterrows():
        sample_xlsx = base / str(r["results_folder"]) / "sample_results.xlsx"
        if not sample_xlsx.exists():
            continue
        dist = load_distance_counts(sample_xlsx)
        if dist.empty:
            continue
        dist = dist.dropna(subset=["read_count"])
        dist["read_count"] = pd.to_numeric(dist["read_count"], errors="coerce").fillna(0).astype(int)
        dist = dist[dist["read_count"] > 0]
        if dist.empty:
            continue

        for d, c in zip(dist["dna_distance"].astype(int), dist["read_count"].astype(int)):
            rows.append(
                {
                    "time": float(r["time"]),
                    "group": str(r["group"]),
                    "dna_distance": int(d),
                    "read_count": int(c),
                }
            )

    if not rows:
        raise ValueError("No read_count data found to build violins.")

    long = pd.DataFrame(rows)

    # Expand by read_count for true violin density.
    expanded = long.loc[long.index.repeat(long["read_count"])].copy()
    expanded = expanded.drop(columns=["read_count"])

    # Jitter plot with seaborn.
    times = sorted(expanded["time"].unique())
    groups = sorted(expanded["group"].unique())
    group_colors = {
        "rep1_WT/E": "#1f77b4",
        "rep3_WT/E": "#ff7f0e",
        "rep2_D3": "#2ca02c",
        "rep3_D3": "#d62728",
    }

    # Use categorical x for stable ordering and jitter.
    expanded["time_str"] = expanded["time"].map(lambda t: f"{t:g}")
    time_order = [f"{t:g}" for t in times]

    fig, ax = plt.subplots(figsize=(16, 7))
    sns.stripplot(
        data=expanded,
        x="time_str",
        y="dna_distance",
        hue="group",
        order=time_order,
        dodge=True,
        jitter=0.25,
        alpha=0.25,
        size=1.5,
        palette=group_colors,
        ax=ax,
    )
    ax.set_xlabel("Time")
    ax.set_ylabel("Mutation distance")
    ax.grid(False)
    ax.legend(title="", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

    fig.tight_layout()
    fig.savefig(out_pdf, dpi=600)
    plt.close(fig)

    print(f"Wrote {out_pdf}")


if __name__ == "__main__":
    main()

