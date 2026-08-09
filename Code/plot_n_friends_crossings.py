"""Histograms of crossing numbers for n-friends found (n = -3..-1, 1..5)."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
CSV_PATH = REPO_ROOT / "data" / "n_friends(master version).csv"
OUT_DIR = REPO_ROOT / "preliminary_results"

N_VALUES = [-3, -2, -1, 1, 2, 3, 4, 5]
BIN_WIDTH = 50


def shared_bins(all_values: np.ndarray, width: int = BIN_WIDTH) -> np.ndarray:
    lo = int(np.floor(all_values.min() / width) * width)
    hi = int(np.ceil(all_values.max() / width) * width)
    return np.arange(lo, hi + width, width)


def plot_histogram(values: np.ndarray, n: int, bins: np.ndarray, out_path: Path) -> None:
    q25, q50, q75 = np.percentile(values, [25, 50, 75])
    mean = values.mean()

    fig, ax = plt.subplots(figsize=(9, 5.5))
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


def plot_overlay(groups: dict[int, np.ndarray], bins: np.ndarray, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
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
    ax.legend(ncol=2, fontsize=9)
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


if __name__ == "__main__":
    main()
