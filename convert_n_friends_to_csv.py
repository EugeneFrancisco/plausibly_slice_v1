"""Convert saved n-friend link diagrams to a CSV file."""

import csv
from pathlib import Path

from plink import LinkManager


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "results" / "n_friends"
OUTPUT = HERE / "results" / "n_friends.csv"
FIELDNAMES = ["name", "n", "PD_code"]


def pd_code(link_file):
    manager = LinkManager()
    manager._from_string(link_file.read_text())
    return manager.PD_code()


def main():
    rows_written = 0

    with OUTPUT.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=FIELDNAMES)
        writer.writeheader()

        for knot_directory in sorted(SOURCE.iterdir()):
            if not knot_directory.is_dir():
                continue

            friend_directories = sorted(
                knot_directory.glob("*_friends"),
                key=lambda path: int(path.name.removesuffix("_friends")),
            )
            for friend_directory in friend_directories:
                n = int(friend_directory.name.removesuffix("_friends"))
                for link_file in sorted(friend_directory.glob("*.lnk")):
                    writer.writerow(
                        {
                            "name": knot_directory.name,
                            "n": n,
                            "PD_code": pd_code(link_file),
                        }
                    )
                    rows_written += 1

    print(f"Wrote {rows_written} friends to {OUTPUT}.")


if __name__ == "__main__":
    main()
