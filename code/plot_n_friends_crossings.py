"""Histograms of crossing numbers for n-friends found (n = -3..-1, 1..5)."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
CSV_PATH = REPO_ROOT / "Data" / "n_friends(master version).csv"
OUT_DIR = REPO_ROOT / "Data" / "results"

N_VALUES = [-3, -2, -1, 1, 2, 3, 4, 5]
BIN_WIDTH = 50

plt.rcParams.update({
    "font.size": 15,
    "axes.titlesize": 18,
    "axes.labelsize": 17,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
    "figure.titlesize": 22,
})


def shared_bins(all_values: np.ndarray, width: int = BIN_WIDTH) -> np.ndarray:
    lo = int(np.floor(all_values.min() / width) * width)
    hi = int(np.ceil(all_values.max() / width) * width)
    return np.arange(lo, hi + width, width)


def plot_histogram(values: np.ndarray, n: int, bins: np.ndarray, out_path: Path) -> None:
    q25, q50, q75 = np.percentile(values, [25, 50, 75])
    mean = values.mean()

    fig, ax = plt.subplots(figsize=(10.5, 6.5))
    ax.hist(values, bins=bins, color="#4C72B0", edgecolor="white")

    marks = [
        ("25th pct", q25, "#DD8452", "--"),
        ("median",   q50, "#55A868", "--"),
        ("75th pct", q75, "#DD8452", "--"),
        ("mean",     mean, "#C44E52", "-"),
    ]
    for label, x, color, style in marks:
        ax.axvline(x, color=color, linestyle=style, linewidth=1.8,
                   label=f"{label} = {x:.2f}")

    ax.set_xlabel("Number of crossings")
    ax.set_ylabel("Count")
    ax.set_title(f"Number of crossings of {n}-friend")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_tiled(groups: dict[int, np.ndarray], bins: np.ndarray, out_path: Path) -> None:
    ns_sorted = sorted(groups)
    ncols = 2
    nrows = int(np.ceil(len(ns_sorted) / ncols))

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(7.0 * ncols, 3.8 * nrows),
        sharex=True, sharey=True,
    )
    axes_flat = np.atleast_1d(axes).flatten()

    mark_specs = [
        ("25th pct", "#DD8452", "--"),
        ("median",   "#55A868", "--"),
        ("75th pct", "#DD8452", ":"),
        ("mean",     "#C44E52", "-"),
    ]

    for ax, n in zip(axes_flat, ns_sorted):
        vals = groups[n]
        q25, q50, q75 = np.percentile(vals, [25, 50, 75])
        mean = vals.mean()
        stats = [q25, q50, q75, mean]

        ax.hist(vals, bins=bins, color="#4C72B0", edgecolor="white")
        for (_, color, style), x in zip(mark_specs, stats):
            ax.axvline(x, color=color, linestyle=style, linewidth=1.4)

        ax.set_title(f"n = {n}  (N = {len(vals)})", fontsize=17)
        stats_text = (
            f"25th = {q25:.1f}\n"
            f"med  = {q50:.1f}\n"
            f"75th = {q75:.1f}\n"
            f"mean = {mean:.1f}"
        )
        ax.text(
            0.98, 0.97, stats_text,
            transform=ax.transAxes,
            ha="right", va="top",
            fontsize=13, family="monospace",
            bbox=dict(boxstyle="round,pad=0.3",
                      facecolor="white", edgecolor="0.7", alpha=0.85),
        )

    # Hide any unused axes
    for ax in axes_flat[len(ns_sorted):]:
        ax.set_visible(False)

    # Shared axis labels
    fig.supxlabel("Number of crossings", y=0.035)
    fig.supylabel("Count")

    # One legend for the marker line styles, below the shared x-axis label
    legend_handles = [
        plt.Line2D([0], [0], color=color, linestyle=style, linewidth=2.2, label=label)
        for label, color, style in mark_specs
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=len(mark_specs),
        frameon=False,
        bbox_to_anchor=(0.5, -0.005),
    )

    fig.suptitle("Number of crossings of n-friends", fontsize=22)
    fig.tight_layout(rect=(0, 0.055, 1, 0.975))
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_overlay(groups: dict[int, np.ndarray], bins: np.ndarray, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 7))
    cmap = plt.get_cmap("viridis")
    ns_sorted = sorted(groups)
    for i, n in enumerate(ns_sorted):
        color = cmap(i / max(1, len(ns_sorted) - 1))
        ax.hist(
            groups[n],
            bins=bins,
            histtype="step",
            linewidth=1.8,
            color=color,
            label=f"n = {n:>2}  (N = {len(groups[n])})",
        )
    ax.set_xlabel("Number of crossings")
    ax.set_ylabel("Count")
    ax.set_title("Crossing-number distributions across n")
    ax.legend(ncol=2, fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(CSV_PATH)

    groups: dict[int, np.ndarray] = {}
    for n in N_VALUES:
        sub = df.loc[df["n"] == n, "num_crossings"].to_numpy()
        if sub.size:
            groups[n] = sub

    all_vals = np.concatenate(list(groups.values()))
    bins = shared_bins(all_vals)

    for n, vals in groups.items():
        tag = f"neg{abs(n)}" if n < 0 else str(n)
        out_path = OUT_DIR / f"crossings_hist_n{tag}.png"
        plot_histogram(vals, n, bins, out_path)
        print(f"[wrote] {out_path.relative_to(REPO_ROOT)}  (N={len(vals)})")

    overlay_path = OUT_DIR / "crossings_hist_overlay.png"
    plot_overlay(groups, bins, overlay_path)
    print(f"[wrote] {overlay_path.relative_to(REPO_ROOT)}")

    tiled_path = OUT_DIR / "crossings_hist_tiled.png"
    plot_tiled(groups, bins, tiled_path)
    print(f"[wrote] {tiled_path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
