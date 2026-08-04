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
REQUIRED_FIELDS = {
    "n",
    "n_friend_name",
    "n_friend_index",
    "knot_PD_code",
    "n_friend_PD_code",
}


def read_csv(path):
    with path.open(newline="") as file:
        lines = file.readlines()

    for header_index, line in enumerate(lines):
        if set(next(csv.reader([line]))) >= REQUIRED_FIELDS:
            break
    else:
        raise ValueError(
            f"{path} does not use the n-friend CSV schema."
        )

    reader = csv.DictReader(lines[header_index:])
    return lines[:header_index], reader.fieldnames, list(reader)


def n_surgery(exterior, n):
    filled = snappy.ManifoldHP(exterior)
    filled.dehn_fill(n_surgery_slope(filled, n))
    return filled


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
    corrections = 0

    for row_number, row in enumerate(rows, start=1):
        name = row["n_friend_name"]
        n = int(row["n"])
        if n < 0:
            raise ValueError(f"Row {row_number} has a negative n: {n}.")

        key = (row["n_friend_index"], n)
        if key not in targets:
            original = snappy.Link(
                ast.literal_eval(row["n_friend_PD_code"])
            )
            targets[key] = n_surgery(original.exterior(), n)
        target = targets[key]

        friend = snappy.Link(ast.literal_eval(row["knot_PD_code"]))
        if _closed_isometric(target, n_surgery(friend.exterior(), n)):
            print(f"Verified row {row_number}: {name}, n={n}.")
            continue

        mirrored_friend = friend.mirror()
        if _closed_isometric(
            target, n_surgery(mirrored_friend.exterior(), n)
        ):
            row["knot_PD_code"] = str(mirrored_friend.PD_code())
            row["num_crossings"] = str(len(mirrored_friend.crossings))
            row["volume"] = str(float(mirrored_friend.exterior().volume()))
            row["verification"] = "True"
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
