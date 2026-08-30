#!/usr/bin/env python3
"""Reproducible controlled-benchmark plots from summary TSVs only."""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SUM = ROOT / "benchmark/summaries/block12"
OUT = ROOT / "benchmark/plots/block12"
OUT.mkdir(parents=True, exist_ok=True)
DEPTHS = [100, 250, 500, 1000, 2000, 5000]
TARGETS = [0.005, 0.015, 0.035, 0.075, 0.25, 0.75]
COLORS = {"gatk": "#2166ac", "mutserve2": "#b2182b", "A": "#666666", "B": "#1b9e77", "C": "#d95f02"}


def read(name):
    with (SUM / name).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def mean(rows, field):
    values = [float(r[field]) for r in rows if r.get(field) not in {None, "", "NA"}]
    return np.mean(values) if values else np.nan


def save(fig, name):
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    caller = read("block12_caller_summary.tsv")
    variants = read("block12_variant_metrics.tsv")
    numt = read("block12_numt_summary.tsv")
    audit = read("block12_numt_observability_audit.tsv")
    circular = read("block12_circular_summary.tsv")
    depth = read("block12_depth_summary.tsv")
    runs = read("block12_run_metrics.tsv")

    # A. Recall versus nominal target VAF, compact panels by depth.
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharex=True, sharey=True)
    for ax, d in zip(axes.flat, DEPTHS):
        for c in ("gatk", "mutserve2"):
            y = [mean([r for r in caller if int(r["depth"]) == d and r["caller"] == c and float(r["target_vaf"]) == v], "recall") for v in TARGETS]
            ax.plot(TARGETS, y, marker="o", label=c, color=COLORS[c])
        ax.set_title(f"{d}×"); ax.set_xscale("log"); ax.set_ylim(-.03, 1.03); ax.grid(alpha=.25)
    axes[0, 0].legend(); fig.supxlabel("Nominal VAF target (log scale)"); fig.supylabel("Recall")
    save(fig, "block12_recall_vs_vaf")

    # B. VAF error against achieved VAF for baseline canonical callsets.
    fig, ax = plt.subplots(figsize=(8, 5))
    for c, source in (("gatk", "gatk_canonical_A"), ("mutserve2", "mutserve2_canonical_A")):
        subset = [r for r in variants if r["callset_source"] == source and r["absolute_vaf_error"] != "NA"]
        points = defaultdict(list)
        for r in subset: points[float(r["achieved_vaf"])].append(float(r["absolute_vaf_error"]))
        xs = sorted(points); ys = [np.mean(points[x]) for x in xs]
        ax.plot(xs, ys, marker="o", label=c, color=COLORS[c])
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("Achieved input VAF"); ax.set_ylabel("Mean absolute VAF error"); ax.grid(alpha=.25); ax.legend()
    save(fig, "block12_absolute_vaf_error")

    # C. Observable NUMT-associated FPs only.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, c in zip(axes, ("gatk", "mutserve2")):
        x = np.arange(len(DEPTHS)); width = .24
        for i, strategy in enumerate(("A", "B", "C")):
            vals = [sum(int(r["observable_numt_fp"]) for r in numt if int(r["depth"]) == d and r["caller"] == c and r["numt_strategy"] == strategy) for d in DEPTHS]
            ax.bar(x + (i-1)*width, vals, width, label=strategy, color=COLORS[strategy])
        ax.set_title(c); ax.set_xticks(x, DEPTHS, rotation=35); ax.set_xlabel("Depth"); ax.grid(axis="y", alpha=.25)
    axes[0].set_ylabel("Observable challenge FP count (3 replicates)"); axes[0].legend(title="Strategy")
    save(fig, "block12_numt_fp_by_strategy")

    # D. Observability (deduplicate caller-expanded audit rows).
    unique = {(r["benchmark_id"], r["truth_id"]): r for r in audit}.values()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for context, marker in (("NUMT_LIKE", "o"), ("AMBIGUOUS_CONTEXT", "s")):
        vals = []
        for d in DEPTHS:
            q = [r for r in unique if int(r["depth"]) == d and r["numt_context"] == context]
            vals.append(100 * sum(r["observable_numt_challenge"] == "true" for r in q) / len(q))
        ax.plot(DEPTHS, vals, marker=marker, label=context)
    ax.set_xscale("log"); ax.set_ylim(0, 105); ax.set_xlabel("Depth"); ax.set_ylabel("Observable challenge (%)"); ax.grid(alpha=.25); ax.legend()
    save(fig, "block12_numt_observability_vs_depth")

    # E. Boundary/non-boundary recall by standard versus shifted-aware.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, boundary in zip(axes, ("BOUNDARY", "NON_BOUNDARY")):
        for mode, color in (("standard_only", "#762a83"), ("standard_plus_shifted", "#1b7837")):
            vals = [mean([r for r in circular if int(r["depth"]) == d and r["boundary_class"] == boundary and r["circular_mode"] == mode], "recall") for d in DEPTHS]
            ax.plot(DEPTHS, vals, marker="o", label=mode, color=color)
        ax.set_xscale("log"); ax.set_title(boundary); ax.set_xlabel("Depth"); ax.grid(alpha=.25)
    axes[0].set_ylabel("Recall"); axes[0].legend(fontsize=8)
    save(fig, "block12_boundary_recall")

    # F. Recall and F1 versus depth for baseline Strategy A.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, metric_name in zip(axes, ("recall", "f1")):
        for c in ("gatk", "mutserve2"):
            vals = [mean([r for r in depth if int(r["depth"]) == d and r["caller"] == c and r["numt_strategy"] == "A"], metric_name) for d in DEPTHS]
            ax.plot(DEPTHS, vals, marker="o", label=c, color=COLORS[c])
        ax.set_xscale("log"); ax.set_title(metric_name.upper()); ax.set_xlabel("Depth"); ax.grid(alpha=.25)
    axes[0].set_ylabel("Metric value"); axes[0].legend()
    save(fig, "block12_accuracy_vs_depth")

    # G/H/I resource scaling from one physical run per replicate/depth.
    for field, ylabel, name in (("wall_clock_seconds", "Wall time (seconds)", "block12_runtime_vs_depth"),
                                ("cpu_hours", "CPU-hours", "block12_cpu_hours_vs_depth"),
                                ("disk_work_bytes_peak", "Peak work bytes", "block12_work_storage_vs_depth")):
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for rep in ("1", "2", "3"):
            vals = [float(next(r[field] for r in runs if r["benchmark_id"] == f"block12_R{rep}_d{d}")) for d in DEPTHS]
            ax.plot(DEPTHS, vals, marker="o", alpha=.55, label=f"R{rep}")
        means = [mean([r for r in runs if r["benchmark_id"].endswith(f"_d{d}")], field) for d in DEPTHS]
        ax.plot(DEPTHS, means, color="black", linewidth=2.5, marker="o", label="mean")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("Depth"); ax.set_ylabel(ylabel); ax.grid(alpha=.25); ax.legend()
        save(fig, name)

    # J. Caller recall heatmaps across VAF/depth.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    for ax, c in zip(axes, ("gatk", "mutserve2")):
        grid = np.array([[mean([r for r in caller if r["caller"] == c and int(r["depth"]) == d and float(r["target_vaf"]) == v], "recall") for v in TARGETS] for d in DEPTHS])
        im = ax.imshow(grid, vmin=0, vmax=1, cmap="viridis", aspect="auto")
        ax.set_title(c); ax.set_xticks(range(6), [str(v) for v in TARGETS], rotation=40); ax.set_yticks(range(6), DEPTHS); ax.set_xlabel("Nominal VAF")
    axes[0].set_ylabel("Depth"); fig.colorbar(im, ax=axes, label="Recall", shrink=.85)
    save(fig, "block12_caller_recall_heatmap")
    print(f"PASS generated={len(list(OUT.glob('*.png')))}_png+{len(list(OUT.glob('*.svg')))}_svg {OUT}")


if __name__ == "__main__":
    main()
