"""Plot crossing-number and n-friend distributions from n_friends.csv."""

import argparse
import csv
import os
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

CACHE_ROOT = Path(tempfile.gettempdir()) / "plausibly_slice_plot_cache"
os.environ.setdefault("MPLCONFIGDIR", str(CACHE_ROOT / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(CACHE_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import PercentFormatter


HERE = Path(__file__).resolve().parent
DEFAULT_INPUT = HERE / "results" / "n_friends.csv"
DEFAULT_OUTPUT = HERE / "results" / "plots"
KNOT_CROSSINGS = re.compile(r"^(?:K)?(\d+)")
N_VALUES = (1, 2, 3)


def load_data(csv_path):
    """Return friend crossing counts grouped by searched-knot crossings."""
    friend_crossings = defaultdict(list)
    n_counts = defaultdict(Counter)

    with csv_path.open(newline="") as file:
        for line_number, row in enumerate(csv.DictReader(file), start=2):
            match = KNOT_CROSSINGS.match(row["n_friend_name"])
            if match is None:
                raise ValueError(
                    f"Could not read crossing number on line {line_number}: "
                    f"{row['n_friend_name']!r}"
                )

            searched_crossings = int(match.group(1))
            n = int(row["n"])
            if n not in N_VALUES:
                raise ValueError(f"Unexpected n={n} on line {line_number}")

            friend_crossings[searched_crossings].append(int(row["num_crossings"]))
            n_counts[searched_crossings][n] += 1

    return dict(friend_crossings), dict(n_counts)


def histogram_edges(values):
    """Choose readable integer-aligned bins, capped at 35 per panel."""
    edges = np.histogram_bin_edges(values, bins="fd")
    bin_count = min(max(len(edges) - 1, 5), 35)
    low, high = min(values), max(values)
    if low == high:
        return np.array([low - 0.5, high + 0.5])
    return np.linspace(low - 0.5, high + 0.5, bin_count + 1)


def plot_crossing_histograms(friend_crossings, output_path):
    crossings = sorted(friend_crossings)
    columns = 2
    rows = (len(crossings) + columns - 1) // columns
    figure, axes = plt.subplots(rows, columns, figsize=(13, 4.2 * rows))
    axes = np.atleast_1d(axes).ravel()

    for axis, crossing in zip(axes, crossings):
        values = friend_crossings[crossing]
        axis.hist(
            values,
            bins=histogram_edges(values),
            color="#4472C4",
            edgecolor="white",
            linewidth=0.7,
        )
        median = float(np.median(values))
        axis.axvline(median, color="#C44E52", linestyle="--", linewidth=1.5)
        axis.set_title(
            f"Friends of {crossing}-crossing knots\n"
            f"{len(values):,} friends; median = {median:g}"
        )
        axis.set_xlabel("Crossings in friend diagram")
        axis.set_ylabel("Number of friends")
        axis.grid(axis="y", alpha=0.25)

    for axis in axes[len(crossings) :]:
        axis.remove()

    figure.suptitle(
        "Distribution of friend crossing numbers (n = 1, 2, and 3 combined)",
        fontsize=15,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.97))
    figure.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(figure)


def plot_n_friend_proportions(n_counts, output_path):
    crossings = sorted(n_counts)
    x_positions = np.arange(len(crossings))
    width = 0.24
    colors = ("#4C78A8", "#F58518", "#54A24B")

    figure, axis = plt.subplots(figsize=(11, 6.5))
    for index, (n, color) in enumerate(zip(N_VALUES, colors)):
        counts = np.array([n_counts[crossing][n] for crossing in crossings])
        totals = np.array([sum(n_counts[crossing].values()) for crossing in crossings])
        proportions = counts / totals
        bars = axis.bar(
            x_positions + (index - 1) * width,
            proportions,
            width,
            label=f"{n}-friends",
            color=color,
        )
        axis.bar_label(
            bars,
            labels=[f"{count:,}" for count in counts],
            padding=3,
            fontsize=9,
        )

    sample_sizes = [sum(n_counts[crossing].values()) for crossing in crossings]
    axis.set_xticks(
        x_positions,
        [f"{crossing}\n(total {total:,})" for crossing, total in zip(crossings, sample_sizes)],
    )
    axis.set_xlabel("Crossings in knot being searched")
    axis.set_ylabel("Proportion of friends found")
    axis.set_title("Proportion of 1-, 2-, and 3-friends by knot crossing number")
    axis.yaxis.set_major_formatter(PercentFormatter(1.0))
    axis.set_ylim(0, max(0.75, axis.get_ylim()[1] * 1.08))
    axis.legend(title="Friend type", ncols=3, loc="upper center")
    axis.grid(axis="y", alpha=0.25)
    axis.set_axisbelow(True)
    axis.text(
        0.99,
        0.98,
        "Bar labels are raw counts",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        color="dimgray",
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    friend_crossings, n_counts = load_data(args.input)
    if not friend_crossings:
        raise ValueError(f"No data found in {args.input}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    histogram_path = args.output_dir / "friend_crossing_histograms.png"
    proportions_path = args.output_dir / "n_friend_proportions.png"
    plot_crossing_histograms(friend_crossings, histogram_path)
    plot_n_friend_proportions(n_counts, proportions_path)

    print(f"Read {sum(map(len, friend_crossings.values())):,} friends.")
    print(f"Wrote {histogram_path}")
    print(f"Wrote {proportions_path}")


if __name__ == "__main__":
    main()
