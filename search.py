"""Search for n-friends of knots with unknown smooth slice status."""

import ast
import csv
from pathlib import Path

import snappy

from find_n_friends import find_common_n_surgery_via_words

# Only search the first SEARCH_LENGTH of the unknown knots.
SEARCH_START = 101
# Exclusive
SEARCH_END = 200

MIN_N = 1
MAX_N = 5

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "results" / "n_friends.csv"
FIELDNAMES = ["name", "n", "PD_code"]


def save_friends(name, n, friends):
    rows = []
    if OUTPUT.exists():
        with OUTPUT.open(newline="") as file:
            rows = list(csv.DictReader(file))

    rows = [
        row
        for row in rows
        if not (row["name"] == name and int(row["n"]) == n)
    ]
    rows.extend(
        {
            "name": name,
            "n": n,
            "PD_code": friend.PD_code(min_strand_index=1),
        }
        for friend in friends
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def main():
    with open(HERE / "data" / "plausibly_unknown.csv", newline="") as file:
        unknown_knots = [
            knot for knot in csv.DictReader(file) if knot["slice"] == "0"
        ]

    for knot in unknown_knots[SEARCH_START:SEARCH_END]:
        name = knot["name"]
        exterior = snappy.Link(ast.literal_eval(knot["PD_codes"])).exterior()

        for n in range(MIN_N, MAX_N + 1):
            print(f"Searching for {n}-friends of {name}.")
            friends = find_common_n_surgery_via_words(exterior, n) or []

            friend_links = []
            for friend in friends:
                friend_exterior = snappy.Manifold(friend[3])
                friend_links.append(friend_exterior.exterior_to_link())
            save_friends(name, n, friend_links)
            print(f"Saved {len(friends)} friends.")


if __name__ == "__main__":
    main()
