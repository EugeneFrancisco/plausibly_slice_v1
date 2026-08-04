"""Search for n-friends or n-RBG links."""

import argparse
import ast
import csv
import importlib.util
import os
from pathlib import Path

import snappy

from find_n_friends import (
    _closed_isometric,
    find_common_n_surgery_via_words,
    n_surgery_slope,
)

MIN_N = 4
MAX_N = 5

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "results" / "n_friends.csv"
FRIEND_SEARCH_PROGRESS = HERE / "results" / "n_friends_4_5_progress.csv"
RBG_OUTPUT = HERE / "results" / "rbg_links"
S_INVARIANT_FIELDS = ("s", "s_2", "s_3")
FIELDNAMES = [
    "id_num",
    "num_crossings",
    "volume",
    "n",
    "verification",
    "n_friend_name",
    "n_friend_index",
    "knot_PD_code",
    "n_friend_PD_code",
]


def n_surgery(exterior, n):
    filled = snappy.ManifoldHP(exterior)
    filled.dehn_fill(n_surgery_slope(filled, n))
    return filled


def orient_and_verify_friend(friend, target, n):
    """Orient a friend so its positive n-surgery matches the target."""
    if _closed_isometric(target, n_surgery(friend.exterior(), n)):
        return friend, True

    mirrored_friend = friend.mirror()
    if _closed_isometric(
        target, n_surgery(mirrored_friend.exterior(), n)
    ):
        return mirrored_friend, True
    return friend, False


def read_output():
    if not OUTPUT.exists():
        return list(FIELDNAMES), []

    with OUTPUT.open(newline="") as file:
        reader = csv.DictReader(file)
        fieldnames = reader.fieldnames or []
        rows = list(reader)
    missing = [field for field in FIELDNAMES if field not in fieldnames]
    if missing:
        raise ValueError(
            f"{OUTPUT} does not use the n-friend schema; missing {missing}."
        )
    return fieldnames, rows


def save_friends(name, source_index, original_pd_code, exterior, n, friends):
    """Append new friends to the output CSV and return the number added."""
    fieldnames, rows = read_output()
    if not friends:
        return 0

    target = n_surgery(exterior, n)
    next_id = max((int(row["id_num"]) for row in rows), default=0) + 1
    existing_pd_codes = {
        row["knot_PD_code"]
        for row in rows
        if row["n_friend_index"] == str(source_index)
        and row["n"] == str(n)
    }
    new_rows = []
    for friend in friends:
        friend, verified = orient_and_verify_friend(friend, target, n)
        friend_exterior = friend.exterior()
        friend_pd_code = str(friend.PD_code())
        if friend_pd_code in existing_pd_codes:
            continue
        row = {
            "id_num": next_id + len(new_rows),
            "num_crossings": len(friend.crossings),
            "volume": float(friend_exterior.volume()),
            "n": n,
            "verification": verified,
            "n_friend_name": name,
            "n_friend_index": source_index,
            "knot_PD_code": friend_pd_code,
            "n_friend_PD_code": str(original_pd_code),
        }
        for field in fieldnames:
            row.setdefault(field, "na")
        new_rows.append(row)
        existing_pd_codes.add(friend_pd_code)

    if not new_rows:
        return 0

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("a", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writerows(new_rows)
        file.flush()
        os.fsync(file.fileno())
    return len(new_rows)


def friend_links(search_results):
    rows = []
    for result in search_results:
        friend_exterior = snappy.Manifold(result[3])
        rows.append(friend_exterior.exterior_to_link())
    return rows


def read_friend_search_progress(rows):
    """Return completed (base-knot index, n) searches."""
    completed = {
        (row["n_friend_index"], int(row["n"]))
        for row in rows
        if row["n"] in {str(MIN_N), str(MAX_N)}
    }
    if FRIEND_SEARCH_PROGRESS.exists():
        with FRIEND_SEARCH_PROGRESS.open(newline="") as file:
            for row in csv.DictReader(file):
                completed.add((row["n_friend_index"], int(row["n"])))
    return completed


def mark_friend_search_complete(source_index, name, n):
    """Record and sync one completed search, including zero-result searches."""
    exists = FRIEND_SEARCH_PROGRESS.exists()
    FRIEND_SEARCH_PROGRESS.parent.mkdir(parents=True, exist_ok=True)
    with FRIEND_SEARCH_PROGRESS.open("a", newline="") as file:
        writer = csv.DictWriter(
            file, fieldnames=("n_friend_index", "n_friend_name", "n")
        )
        if not exists or FRIEND_SEARCH_PROGRESS.stat().st_size == 0:
            writer.writeheader()
        writer.writerow(
            {
                "n_friend_index": source_index,
                "n_friend_name": name,
                "n": n,
            }
        )
        file.flush()
        os.fsync(file.fileno())


def base_knots(rows):
    """Return each base knot represented in the n-friends CSV once."""
    knots = {}
    for row in rows:
        source_index = row["n_friend_index"]
        knots.setdefault(
            source_index,
            (row["n_friend_name"], row["n_friend_PD_code"]),
        )
    return sorted(knots.items(), key=lambda item: int(item[0]))


def search_n_friends():
    """Find missing 4- and 5-friends of base knots in the output CSV."""
    _, rows = read_output()
    knots = base_knots(rows)
    completed = read_friend_search_progress(rows)
    searches = [
        (source_index, name, pd_code, n)
        for source_index, (name, pd_code) in knots
        for n in range(MIN_N, MAX_N + 1)
        if (source_index, n) not in completed
    ]

    print(
        f"Loaded {len(knots)} base knots; "
        f"{len(searches)} of {len(knots) * 2} searches remain.",
        flush=True,
    )
    for position, (source_index, name, pd_code, n) in enumerate(
        searches, start=1
    ):
        original_pd_code = ast.literal_eval(pd_code)
        exterior = snappy.Link(original_pd_code).exterior()
        print(
            f"[{position}/{len(searches)}] Searching for {n}-friends "
            f"of {name} (base index {source_index}).",
            flush=True,
        )
        results = find_common_n_surgery_via_words(exterior, n) or []
        friends = friend_links(results)
        added = save_friends(
            name,
            source_index,
            original_pd_code,
            exterior,
            n,
            friends,
        )
        mark_friend_search_complete(source_index, name, n)
        print(
            f"Found {len(friends)} friends; appended {added}. "
            "Checkpoint saved.",
            flush=True,
        )


def has_slice_obstructing_s_invariant(row):
    """Return whether a recorded s-invariant is at least 2."""
    for field in S_INVARIANT_FIELDS:
        try:
            if float(row.get(field, "")) >= 2:
                return True
        except (TypeError, ValueError):
            continue
    return False


def load_n_rbg_search():
    """Load the n-RBG search function from n-rbg.py."""
    path = HERE / "n-rbg.py"
    spec = importlib.util.spec_from_file_location("n_rbg", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load {path}.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.forms_special_NRBG_link


def save_rbg_link(rbg_link, row, output_dir=RBG_OUTPUT):
    """Save an n-RBG link as a SnapPy Link Editor projection."""
    from plink import LinkManager

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    filename = (
        f"{row['n_friend_name']}_{row['n']}_friend_"
        f"row_{row['id_num']}_RBG.lnk"
    )
    path = output_dir / filename
    manager = LinkManager()
    rbg_link.link.view(manager)
    path.write_text(manager.SnapPea_projection_file())
    save_rbg_manifest_entry(rbg_link, row, path)
    return path


def save_rbg_manifest_entry(rbg_link, row, link_path):
    """Record the framing data needed to verify a saved RBG link."""
    fieldnames = [
        "original_knot",
        "n",
        "friend_row",
        "framings",
        "link_file",
    ]
    manifest = link_path.parent / "certificates.csv"
    records = []
    if manifest.exists():
        with manifest.open(newline="") as file:
            records = list(csv.DictReader(file))

    record = {
        "original_knot": row["n_friend_name"],
        "n": row["n"],
        "friend_row": row["id_num"],
        "framings": str(rbg_link.framings),
        "link_file": link_path.name,
    }
    records = [
        old for old in records if old.get("link_file") != link_path.name
    ]
    records.append(record)
    with manifest.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def search_n_rbg_links(csv_path=OUTPUT, output_dir=RBG_OUTPUT):
    """Search obstructing friend pairs for n-special RBG links."""
    with Path(csv_path).open(newline="") as file:
        reader = csv.DictReader(file)
        fieldnames = set(reader.fieldnames or [])
        required = {
            "n",
            "n_friend_name",
            "knot_PD_code",
            "n_friend_PD_code",
        }
        missing = required - fieldnames
        if missing:
            raise ValueError(
                f"{csv_path} does not use the n-friend schema; "
                f"missing {sorted(missing)}."
            )
        rows = [
            row for row in reader if has_slice_obstructing_s_invariant(row)
        ]

    print(f"Found {len(rows)} friend pairs with s-invariant at least 2.")
    forms_special_n_rbg_link = load_n_rbg_search()
    for row in rows:
        name = row["n_friend_name"]
        n = int(row["n"])
        friend = snappy.Link(ast.literal_eval(row["knot_PD_code"]))
        original = snappy.Link(ast.literal_eval(row["n_friend_PD_code"]))

        print(f"Searching for a special {n}-RBG link for {name}.")
        rbg_link = forms_special_n_rbg_link(
            n, friend.exterior(), original.exterior()
        )
        if rbg_link is not None:
            print(f"{name} was found to not be slice.")
            path = save_rbg_link(rbg_link, row, output_dir)
            print(f"Saved RBG link to {path}.")


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=("friends", "rbg"),
        default="friends",
        help="search for n-friends or search saved pairs for n-RBG links",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.mode == "friends":
        search_n_friends()
    else:
        search_n_rbg_links()


if __name__ == "__main__":
    main()
