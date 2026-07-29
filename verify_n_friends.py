"""Verify the n-friend pairs stored in a CSV file."""

import argparse
import ast
import csv
from pathlib import Path
import tempfile

import snappy

from find_n_friends import _closed_isometric, n_surgery_slope


HERE = Path(__file__).resolve().parent
DEFAULT_CSV = HERE / "results" / "n_friends.csv"
KNOT_DATA = HERE / "data" / "plausibly_unknown.csv"
REQUIRED_FIELDS = {"name", "n", "PD_code"}


def read_csv(path):
    with path.open(newline="") as file:
        lines = file.readlines()

    for header_index, line in enumerate(lines):
        if set(next(csv.reader([line]))) >= REQUIRED_FIELDS:
            break
    else:
        raise ValueError(
            f"{path} must have the columns name, n, and PD_code."
        )

    reader = csv.DictReader(lines[header_index:])
    return lines[:header_index], reader.fieldnames, list(reader)


def n_surgery(exterior, n):
    filled = snappy.ManifoldHP(exterior)
    filled.dehn_fill(n_surgery_slope(filled, n))
    return filled


def load_base_pd_codes():
    with KNOT_DATA.open(newline="") as file:
        return {
            row["name"]: ast.literal_eval(row["PD_codes"])
            for row in csv.DictReader(file)
        }


def base_knot_exterior(name, pd_codes):
    try:
        return snappy.Manifold(name)
    except (OSError, ValueError):
        if name not in pd_codes:
            raise ValueError(
                f"Could not load {name} by name or find its PD code in "
                f"{KNOT_DATA}."
            )
        return snappy.Link(pd_codes[name]).exterior()


def write_csv(path, preamble, fieldnames, rows):
    with tempfile.NamedTemporaryFile(
        "w", newline="", dir=path.parent, delete=False
    ) as temporary_file:
        temporary_path = Path(temporary_file.name)
        temporary_file.writelines(preamble)
        writer = csv.DictWriter(temporary_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    temporary_path.chmod(path.stat().st_mode)
    temporary_path.replace(path)


def verify_n_friends(csv_path):
    """Verify every pair and replace mirrored friend diagrams in the CSV."""
    path = Path(csv_path)
    preamble, fieldnames, rows = read_csv(path)
    targets = {}
    base_exteriors = {}
    base_pd_codes = load_base_pd_codes()
    corrections = 0

    for row_number, row in enumerate(rows, start=1):
        name = row["name"]
        n = int(row["n"])
        if n < 0:
            raise ValueError(f"Row {row_number} has a negative n: {n}.")

        key = (name, n)
        if name not in base_exteriors:
            base_exteriors[name] = base_knot_exterior(name, base_pd_codes)
        if key not in targets:
            targets[key] = n_surgery(base_exteriors[name], n)
        target = targets[key]

        friend = snappy.Link(ast.literal_eval(row["PD_code"]))
        if _closed_isometric(target, n_surgery(friend.exterior(), n)):
            print(f"Verified row {row_number}: {name}, n={n}.")
            continue

        mirrored_friend = friend.mirror()
        if _closed_isometric(
            target, n_surgery(mirrored_friend.exterior(), n)
        ):
            row["PD_code"] = str(
                mirrored_friend.PD_code(min_strand_index=1)
            )
            corrections += 1
            print(f"Corrected mirror in row {row_number}: {name}, n={n}.")
            continue

        raise ValueError(
            f"Row {row_number} is not an {n}-friend of {name}, "
            "even after mirroring its PD code."
        )

    if corrections:
        write_csv(path, preamble, fieldnames, rows)
    print(
        f"Verified {len(rows)} n-friends and corrected "
        f"{corrections} mirrored PD codes."
    )
    return corrections


def main():
    parser = argparse.ArgumentParser(
        description="Verify n-friends and correct mirrored PD codes."
    )
    parser.add_argument(
        "csv_path",
        nargs="?",
        type=Path,
        default=DEFAULT_CSV,
        help=f"CSV to verify (default: {DEFAULT_CSV})",
    )
    args = parser.parse_args()
    verify_n_friends(args.csv_path)


if __name__ == "__main__":
    main()
