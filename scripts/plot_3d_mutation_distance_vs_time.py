#!/usr/bin/env python3
"""
3D plot: mutation DNA distance distribution vs time for selected groups.

Input:
- all_summary.xlsx sheet plotted_reads_gt_1000
- per-sample: <results_folder>/sample_results.xlsx sheet dna_distance_distribution

Plot:
- x: dna_distance (integer bins)
- y: mapped time value (continuous)
- z: fraction of plotted reads at that distance (0..1)

Color:
- constant per (replicate, Starting point) group across timepoints

Dedup:
- if there are two rows with (timepoint=7, replicate=3, Starting point=D3), keep only the first
  as it appears in the summary sheet.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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


def load_distance_distribution(sample_results_xlsx: Path) -> pd.DataFrame:
    df = pd.read_excel(sample_results_xlsx, sheet_name="dna_distance_distribution", engine="openpyxl")
    if df.empty or "dna_distance" not in df.columns:
        return pd.DataFrame(columns=["dna_distance", "fraction"])

    out = df.copy()
    out["dna_distance"] = pd.to_numeric(out["dna_distance"], errors="coerce")
    if "fraction" in out.columns:
        out["fraction"] = pd.to_numeric(out["fraction"], errors="coerce")
    elif "read_count" in out.columns:
        rc = pd.to_numeric(out["read_count"], errors="coerce").fillna(0.0)
        total = float(rc.sum())
        out["fraction"] = rc / total if total > 0 else 0.0
    else:
        out["fraction"] = np.nan

    out = out.dropna(subset=["dna_distance", "fraction"])
    out = out[(out["dna_distance"] >= 0) & (out["fraction"] >= 0)]
    out = out.sort_values("dna_distance", kind="stable")
    return out[["dna_distance", "fraction"]]


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
    ap.add_argument(
        "--sheet",
        default="plotted_reads_gt_1000",
        help="Sheet containing sample metadata",
    )
    ap.add_argument(
        "--out_pdf",
        type=Path,
        default=None,
        help="Output PDF path (default: <base_dir>/plots_3d_mutation_distance/mutation_distance_vs_time_3d.pdf)",
    )
    args = ap.parse_args()

    base = args.base_dir.resolve()
    summary = (args.summary_xlsx or (base / "all_summary.xlsx")).resolve()
    out_pdf = args.out_pdf or (base / "plots_3d_mutation_distance" / "mutation_distance_vs_time_3d.pdf")
    out_pdf.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_excel(summary, sheet_name=args.sheet, engine="openpyxl")
    required = {"sample_id", "results_folder", "timepoint", "replicate", "Starting point"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {args.sheet}: {sorted(missing)}")

    # Filter to requested groups.
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

    # Dedup: (timepoint=7, replicate=3, Starting point=D3) keep first row in sheet order.
    mask_dup = (df["timepoint"] == 7) & (df["replicate"] == 3) & (df["Starting point"] == "D3")
    dup_idx = list(df.index[mask_dup])
    if len(dup_idx) >= 2:
        df = df.drop(index=dup_idx[1:])

    if df.empty:
        raise ValueError("No rows matched requested groups after filtering.")

    # Map timepoint to continuous.
    df["time_mapped"] = df["timepoint"].astype(int).map(lambda t: TIMEPOINT_MAP.get(t, float(t)))

    # Assign colors per group.
    group_colors = {
        (1, "WT/E"): "#1f77b4",
        (3, "WT/E"): "#ff7f0e",
        (2, "D3"): "#2ca02c",
        (3, "D3"): "#d62728",
    }

    # Load all distributions and determine global x support.
    rows = []
    all_distances = set()
    for _, r in df.iterrows():
        results_folder = str(r["results_folder"])
        sample_xlsx = base / results_folder / "sample_results.xlsx"
        if not sample_xlsx.exists():
            continue
        dist = load_distance_distribution(sample_xlsx)
        if dist.empty:
            continue
        rep = int(r["replicate"])
        sp = str(r["Starting point"])
        tm = float(r["time_mapped"])
        for d in dist["dna_distance"].astype(int).tolist():
            all_distances.add(int(d))
        rows.append((rep, sp, tm, dist))

    if not rows:
        raise ValueError("No dna_distance_distribution data found for selected samples.")

    xs = np.array(sorted(all_distances), dtype=float)
    if xs.size == 0:
        raise ValueError("No dna_distance bins found.")

    fig = plt.figure(figsize=(18, 10))
    ax = fig.add_subplot(111, projection="3d")

    # Plot each timepoint curve for each group with constant color.
    for rep, sp, tm, dist in rows:
        g = (rep, sp)
        color = group_colors.get(g, "#444444")
        frac_by_d = dict(zip(dist["dna_distance"].astype(int), dist["fraction"].astype(float)))
        zs = np.array([frac_by_d.get(int(d), 0.0) for d in xs], dtype=float)
        ys = np.full_like(xs, tm, dtype=float)
        ax.plot(xs, ys, zs, color=color, linewidth=2.0, alpha=0.9)

    ax.set_xlabel("Mutation distance")
    ax.set_ylabel("Time")
    ax.set_zlabel("Fraction of plotted reads")
    ax.set_zlim(0.0, 1.0)
    ax.set_zticks(np.arange(0.0, 1.01, 0.2))
    ax.grid(False)
    ax.view_init(elev=22, azim=35)

    fig.tight_layout()
    fig.savefig(out_pdf, dpi=600)
    plt.close(fig)

    print(f"Wrote {out_pdf}")


if __name__ == "__main__":
    main()

