#!/usr/bin/env python3
"""
Build all_summary.xlsx: one row per aa_analysis_results_* folder that has sample_results.xlsx.

Sample id = substring between 'results_' and '_cluster' in the folder basename
(e.g. aa_analysis_results_1-48_F2_n14_cluster_1based -> 1-48_F2_n14).
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

RESULT_DIR_RE = re.compile(r"results_(.+)_cluster")


def sample_id_from_dirname(name: str) -> str | None:
    m = RESULT_DIR_RE.search(name)
    return m.group(1) if m else None


def read_sample_summary(xlsx_path: Path) -> dict:
    """Pull totals from mutation_set_categories sheet."""
    out = {
        "total_plotted_reads": None,
        "shannon_entropy_bits_plotted_sets": None,
        "fraction_shared_reads": None,
        "shared_read_count": None,
        "avg_dna_distance_plotted_sets": None,
        "error": None,
    }
    try:
        df = pd.read_excel(xlsx_path, sheet_name="mutation_set_categories", engine="openpyxl")
    except Exception as e:
        out["error"] = str(e)
        return out

    if df.empty:
        out["error"] = "empty sheet"
        return out

    if "total_plotted_reads" in df.columns:
        out["total_plotted_reads"] = int(df["total_plotted_reads"].iloc[0])

    if "shannon_entropy_bits_plotted_sets" in df.columns:
        out["shannon_entropy_bits_plotted_sets"] = float(df["shannon_entropy_bits_plotted_sets"].iloc[0])

    if "category" in df.columns and "fraction_of_plotted_reads" in df.columns:
        shared = df[df["category"].astype(str) == "Shared"]
        if not shared.empty:
            out["fraction_shared_reads"] = float(shared["fraction_of_plotted_reads"].iloc[0])
        else:
            out["fraction_shared_reads"] = 0.0
        if "read_count" in df.columns:
            shared_c = df[df["category"].astype(str) == "Shared"]
            if not shared_c.empty:
                out["shared_read_count"] = int(shared_c["read_count"].iloc[0])
            else:
                out["shared_read_count"] = 0

    # Average DNA distance over plotted mutation sets (weighted by read_count)
    try:
        dist = pd.read_excel(xlsx_path, sheet_name="dna_distance_distribution", engine="openpyxl")
        if not dist.empty and {"dna_distance", "read_count"}.issubset(dist.columns):
            d = pd.to_numeric(dist["dna_distance"], errors="coerce")
            c = pd.to_numeric(dist["read_count"], errors="coerce")
            m = (d.notna()) & (c.notna()) & (c > 0)
            if m.any():
                total = float(c[m].sum())
                if total > 0:
                    out["avg_dna_distance_plotted_sets"] = float((d[m] * c[m]).sum() / total)
    except Exception:
        # Missing sheet is OK (older results)
        pass

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Aggregate sample_results.xlsx into all_summary.xlsx")
    ap.add_argument(
        "base_dir",
        type=Path,
        nargs="?",
        default=None,
        help="Project root (parent of aa_analysis_results_*). Default: parent of scripts/",
    )
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output path (default: <base_dir>/all_summary.xlsx)",
    )
    args = ap.parse_args()

    script_dir = Path(__file__).resolve().parent
    base = (args.base_dir or script_dir.parent).resolve()
    out_path = args.output or (base / "all_summary.xlsx")

    rows = []
    for d in sorted(base.glob("aa_analysis_results_*")):
        if not d.is_dir():
            continue
        sid = sample_id_from_dirname(d.name)
        if not sid:
            continue
        sr = d / "sample_results.xlsx"
        if not sr.is_file():
            rows.append(
                {
                    "sample_id": sid,
                    "results_folder": d.name,
                    "total_plotted_reads": None,
                    "shannon_entropy_bits_plotted_sets": None,
                    "fraction_shared_reads": None,
                    "shared_read_count": None,
                    "avg_dna_distance_plotted_sets": None,
                    "note": "missing sample_results.xlsx",
                }
            )
            continue

        data = read_sample_summary(sr)
        note = data.pop("error", None)
        rows.append(
            {
                "sample_id": sid,
                "results_folder": d.name,
                "total_plotted_reads": data["total_plotted_reads"],
                "shannon_entropy_bits_plotted_sets": data["shannon_entropy_bits_plotted_sets"],
                "fraction_shared_reads": data["fraction_shared_reads"],
                "shared_read_count": data["shared_read_count"],
                "avg_dna_distance_plotted_sets": data["avg_dna_distance_plotted_sets"],
                "note": note or "",
            }
        )

    df = pd.DataFrame(rows)
    if df.empty:
        print("No aa_analysis_results_* directories found under {0}".format(base), file=sys.stderr)
        sys.exit(1)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    reads = pd.to_numeric(df["total_plotted_reads"], errors="coerce")
    df_high = df[reads > 1000].copy()

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="all_samples", index=False)
        df_high.to_excel(writer, sheet_name="plotted_reads_gt_1000", index=False)

    print(
        "Wrote {0} (all_samples: {1} rows; plotted_reads_gt_1000: {2} rows)".format(
            out_path, len(df), len(df_high)
        )
    )


if __name__ == "__main__":
    main()
