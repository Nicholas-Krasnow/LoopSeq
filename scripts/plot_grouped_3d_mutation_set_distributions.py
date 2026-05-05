#!/usr/bin/env python3
"""
Plot 3D mutation-set distributions grouped by (replicate, Starting point).

For each sample listed in all_summary.xlsx / plotted_reads_gt_1000:
- Load <results_folder>/sample_results.xlsx sheet "plotted_mutation_sets"
- Build a mutation-set -> fraction distribution from read_count values

For each unique (replicate, Starting point) combination:
- Sort samples by increasing timepoint (then sample_id)
- Plot all sample distributions on one 3D axis with a shared mutation-set axis
  x = shared mutation-set index (ordered by pooled abundance in group)
  y = timepoint
  z = fraction of plotted reads
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

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

try:
    import matplotlib.font_manager as _fm
    _preferred = ["Arial", "Liberation Sans", "DejaVu Sans"]
    _available = {f.name for f in _fm.fontManager.ttflist}
    for _name in _preferred:
        if _name in _available:
            plt.rcParams["font.family"] = _name
            break
except Exception:
    pass


def slugify(value: str) -> str:
    clean = "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in str(value))
    while "__" in clean:
        clean = clean.replace("__", "_")
    return clean.strip("_")


def load_sample_distribution(sample_results_xlsx: Path) -> dict[str, float]:
    df = pd.read_excel(sample_results_xlsx, sheet_name="plotted_mutation_sets", engine="openpyxl")
    required = {"mutation_set_label", "read_count"}
    if not required.issubset(df.columns):
        raise ValueError("Missing required columns in plotted_mutation_sets")
    df = df.copy()
    df["read_count"] = pd.to_numeric(df["read_count"], errors="coerce")
    df = df.dropna(subset=["mutation_set_label", "read_count"])
    df = df[df["read_count"] > 0]
    if df.empty:
        return {}
    grouped = df.groupby("mutation_set_label", as_index=False)["read_count"].sum()
    total = float(grouped["read_count"].sum())
    if total <= 0:
        return {}
    grouped["fraction"] = grouped["read_count"] / total
    return dict(zip(grouped["mutation_set_label"].astype(str), grouped["fraction"].astype(float)))


def load_sample_categories(sample_results_xlsx: Path) -> dict[str, str]:
    df = pd.read_excel(sample_results_xlsx, sheet_name="plotted_mutation_sets", engine="openpyxl")
    if not {"mutation_set_label", "category"}.issubset(df.columns):
        return {}
    out = {}
    for row in df.itertuples(index=False):
        label = str(getattr(row, "mutation_set_label"))
        cat = str(getattr(row, "category"))
        out[label] = cat
    return out


def make_group_plot(group_df: pd.DataFrame, base_dir: Path, out_dir: Path) -> Path | None:
    group_df = group_df.sort_values(["timepoint", "sample_id"], kind="stable").reset_index(drop=True)
    replicate = str(group_df["replicate"].iloc[0])
    starting_point = str(group_df["Starting point"].iloc[0])

    series = []
    pooled = {}
    category_votes = {}
    for row in group_df.itertuples(index=False):
        sample_id = str(row.sample_id)
        results_folder = str(row.results_folder)
        tp = float(row.timepoint)
        tp_mapped = TIMEPOINT_MAP.get(int(tp), tp)
        sample_xlsx = base_dir / results_folder / "sample_results.xlsx"
        if not sample_xlsx.is_file():
            continue
        try:
            dist = load_sample_distribution(sample_xlsx)
            cats = load_sample_categories(sample_xlsx)
        except Exception:
            continue
        if not dist:
            continue
        for label, frac in dist.items():
            pooled[label] = pooled.get(label, 0.0) + float(frac)
        for label, cat in cats.items():
            if label not in category_votes:
                category_votes[label] = {}
            category_votes[label][cat] = category_votes[label].get(cat, 0) + 1
        series.append((sample_id, float(tp_mapped), dist))

    if not series:
        return None

    # Sort mutation sets by abundance at final timepoint (ascending),
    # so the most abundant at final timepoint is farthest from origin.
    final_tp = max(tp for _, tp, _ in series)
    final_dists = [dist for _, tp, dist in series if tp == final_tp]
    final_abundance = {}
    for lbl in pooled.keys():
        final_abundance[lbl] = float(sum(d.get(lbl, 0.0) for d in final_dists))
    shared_labels = [k for k, _ in sorted(final_abundance.items(), key=lambda kv: kv[1], reverse=True)]
    labels_by_abundance = list(shared_labels)

    # Stretch and evenly distribute mutation-set positions across the full x-axis.
    x_stretch = 6.5
    x_max = max(1.0, len(shared_labels) * x_stretch)
    if len(shared_labels) == 1:
        xvals = np.array([x_max / 2.0], dtype=float)
    else:
        xvals = np.linspace(1.0, x_max, len(shared_labels), dtype=float)

    fig = plt.figure(figsize=(60, 18))
    ax = fig.add_subplot(111, projection="3d")

    # Color by mutation set category from individual sample plots.
    category_colors = {
        "Shared": "#BFC2BB",
        "D3-specific": "#118040",
        "WT-specific": "#464747",
        "starting_sequence": "#eeeeee",
    }

    color_by_label = {}
    for lbl in labels_by_abundance:
        votes = category_votes.get(lbl, {})
        if votes:
            best_cat = sorted(votes.items(), key=lambda kv: kv[1], reverse=True)[0][0]
        else:
            best_cat = "Shared"
        color_by_label[lbl] = category_colors.get(best_cat, "#BFC2BB")

    # 3D lines by mutation set across time.
    # Each line keeps x fixed (mutation set) and connects its abundance across time.
    series_sorted = sorted(series, key=lambda t: t[1])
    time_values = np.array([tp for _, tp, _ in series_sorted], dtype=float)
    for xv, lbl in zip(xvals, shared_labels):
        zvals = np.array([float(dist.get(lbl, 0.0)) for _, _, dist in series_sorted], dtype=float)
        if np.all(zvals == 0):
            continue
        xs = np.full_like(time_values, xv, dtype=float)
        line_color = color_by_label[lbl]
        ax.plot(
            xs,
            time_values,
            zvals,
            color=line_color,
            linewidth=1.8,
            alpha=0.95,
        )
        # Filled area under each line (to z=0) with same color at 50% transparency.
        verts = [(xv, float(t), float(z)) for t, z in zip(time_values, zvals)]
        verts += [(xv, float(t), 0.0) for t in time_values[::-1]]
        poly = Poly3DCollection([verts], facecolor=line_color, edgecolor="none", alpha=0.5)
        ax.add_collection3d(poly)

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_zlabel("")
    ax.set_title(f"Replicate {replicate}, Starting point {starting_point}")
    # Stretch x-axis span to reduce overlap among dense mutation-set labels.
    ax.set_xlim(1.0, x_max)
    ax.grid(False)
    mapped_times = pd.to_numeric(group_df["timepoint"], errors="coerce").dropna().astype(float).map(
        lambda t: TIMEPOINT_MAP.get(int(t), float(t))
    )
    if not mapped_times.empty:
        ax.set_ylim(mapped_times.min() - 2.0, mapped_times.max() + 2.0)
    # Use a fixed z-axis scale and ticks across all plots.
    ax.set_zlim(0.0, 1.0)
    zticks = np.arange(0.0, 1.01, 0.2)
    ax.set_zticks(zticks)
    ax.set_zticklabels([f"{t:.1f}" for t in zticks])
    ax.zaxis.labelpad = 18
    ax.tick_params(axis="z", pad=8, labelsize=21.12)
    # Make z-axis label read in the same direction as increasing axis values.
    ax.zaxis.set_rotate_label(False)
    ax.zaxis.label.set_rotation(90)
    # View angle: bring timepoint axis to the front and mutation-set axis to the side.
    ax.view_init(elev=22, azim=28)

    # Remove gray background panes.
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        try:
            axis.pane.fill = False
            axis.pane.set_edgecolor((1, 1, 1, 0))
        except Exception:
            pass

    # Label mutation sets directly on x-axis (instead of numeric index).
    xlabels = shared_labels
    ax.set_xticks(xvals)
    ax.set_xticklabels(xlabels, fontsize=6, rotation=90, ha="center", va="top")
    ax.tick_params(axis="x", pad=16)
    ax.xaxis.labelpad = 44
    ax.tick_params(axis="y", labelsize=21.12)

    # Thicken 3D axis lines.
    try:
        ax.xaxis.line.set_linewidth(2.0)
        ax.yaxis.line.set_linewidth(2.0)
        ax.zaxis.line.set_linewidth(2.0)
    except Exception:
        pass

    # Maximize usable plot area for dense x labels while keeping within PDF canvas.
    fig.subplots_adjust(left=0.03, right=0.992, bottom=0.26, top=0.96)
    out_name = f"mutation_set_distribution_3d_rep{slugify(replicate)}_start{slugify(starting_point)}.pdf"
    out_path = out_dir / out_name
    fig.savefig(out_path, dpi=600)
    plt.close(fig)
    return out_path


def main() -> None:
    ap = argparse.ArgumentParser(description="Grouped 3D mutation-set distribution plots")
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
        help="Sheet in summary workbook with sample metadata",
    )
    ap.add_argument(
        "--out_dir",
        type=Path,
        default=None,
        help="Output directory for plots (default: <base_dir>/plots_3d_mutation_set_distributions)",
    )
    args = ap.parse_args()

    base_dir = args.base_dir.resolve()
    summary_xlsx = (args.summary_xlsx or (base_dir / "all_summary.xlsx")).resolve()
    out_dir = (args.out_dir or (base_dir / "plots_3d_mutation_set_distributions")).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_excel(summary_xlsx, sheet_name=args.sheet, engine="openpyxl")
    required = {"sample_id", "results_folder", "timepoint", "replicate", "Starting point"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {args.sheet}: {sorted(missing)}")

    # Keep rows with grouping/timepoint metadata available.
    df = df.dropna(subset=["replicate", "Starting point", "timepoint"]).copy()
    df["timepoint"] = pd.to_numeric(df["timepoint"], errors="coerce")
    df = df.dropna(subset=["timepoint"])
    # Special-case de-duplication requested by user:
    # remove only the second occurrence (workbook order) of
    # timepoint == 7, replicate == 3, Starting point == D3.
    mask_specific = (
        (df["Starting point"].astype(str) == "D3")
        & (df["timepoint"] == 7)
        & (pd.to_numeric(df["replicate"], errors="coerce") == 3)
    )
    match_indices = list(df.index[mask_specific])
    if len(match_indices) >= 2:
        second_idx = match_indices[1]
        df = df.loc[df.index != second_idx].copy()
    if df.empty:
        raise ValueError("No usable rows after filtering for replicate/Starting point/timepoint")

    created = []
    for (_, _), g in df.groupby(["replicate", "Starting point"], sort=True):
        out_path = make_group_plot(g, base_dir=base_dir, out_dir=out_dir)
        if out_path is not None:
            created.append(out_path)

    print(f"Wrote {len(created)} plot(s) to {out_dir}")
    for p in created:
        print(p)


if __name__ == "__main__":
    main()

