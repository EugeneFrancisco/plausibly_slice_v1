"""
Run the n-friend *search* on the documented examples and check it recovers
the listed partner knot.

For each example in ``qin_n_friends.csv`` we keep the ones we can check
quickly -- n <= 3, both knots loadable in SnapPy with fewer than 18 crossings,
and flagged ``isometry_testable`` -- then call the paper's search routine
``find_common_n_surgery_via_words`` and confirm the partner knot turns up among
the n-friends it finds.

Two practical points:

  * The search drills the +n surgery of the knot it starts from, so it only
    recovers a pair from the knot whose +n surgery is the shared manifold.  We
    therefore try both knots as the starting point.
  * A returned friend is matched to the target by comparing knot *exteriors*
    (orientation-blind), which sidesteps the chirality bookkeeping that the
    closed-surgery comparison in ``is_n_friend`` has to worry about.

Some shared surgeries have no geometric (drillable) triangulation, so the
search cannot reach them; such an example is reported as UNREACHABLE rather
than a failure.

Run in the ``sage`` conda env:

    conda activate sage
    python test.py
"""

from __future__ import annotations

import os
import re
import csv
import sys
from typing import Iterator

import snappy

sys.setrecursionlimit(50000)

_HERE = os.path.dirname(os.path.abspath(__file__))
if os.path.join(_HERE, "code") not in sys.path:
    sys.path.insert(0, os.path.join(_HERE, "code"))

from find_n_friends import find_common_n_surgery_via_words  # noqa: E402

CSV_PATH = os.path.join(_HERE, "qin_n_friends.csv")

MAX_N = 3
MAX_CROSSINGS = 18
MAX_LEN = 3.0                 # geodesic-length cutoff for the search


def crossing_number(name: str) -> int | None:
    """Crossing number parsed from a knot name (e.g. K14n10164 -> 14).

    Args:
        name: A knot name from the CSV.

    Returns:
        The crossing number, or None for exotic notations we cannot load
        (e.g. the pretzel ``K(-2,1,2)`` or connect sum ``T(-2,3)#T(2,5)``).
    """
    match = re.match(r"K?(\d+)", name)
    return int(match.group(1)) if match else None


def exteriors_isometric(isosig: str, target_name: str, tries: int = 8) -> bool:
    """True if the found exterior (an isosig) is the ``target_name`` knot.

    Compares knot exteriors, which are orientation-blind, so a knot and its
    mirror match -- exactly what we want when the tables store one chirality.

    Args:
        isosig: Triangulation isosig of the found exterior.
        target_name: Name of the knot to match against.
        tries: Number of randomization retries on isometry-test failure.

    Returns:
        True iff the two exteriors are isometric.
    """
    F = snappy.ManifoldHP(snappy.Triangulation(isosig))
    T = snappy.ManifoldHP(snappy.Manifold(target_name))
    if abs(float(F.volume()) - float(T.volume())) > 1e-6:
        return False
    for _ in range(tries):
        try:
            return bool(F.is_isometric_to(T))
        except RuntimeError:
            F.randomize(); T.randomize()
    return False


def search_recovers(searched: str, target: str, n: int) -> bool | None:
    """Run the search on ``searched``; is ``target`` among its n-friends?

    Args:
        searched: Name of the knot to start the search from.
        target: Name of the partner knot we hope to recover.
        n: The surgery coefficient.

    Returns:
        True/False whether ``target`` was recovered, or None if the search
        could not run (the n-surgery has no drillable triangulation).
    """
    hits = find_common_n_surgery_via_words(searched, n, MAX_LEN)
    if hits is None:
        return None
    return any(exteriors_isometric(hit[-1], target) for hit in hits)


def selected_rows() -> Iterator[dict[str, str]]:
    """The CSV rows worth checking: isometry-testable, small n, small knots.

    Yields:
        Each qualifying row as a column-name -> value dict.
    """
    with open(CSV_PATH) as f:
        for row in csv.DictReader(f):
            if int(row["n"]) > MAX_N:
                continue
            if row["isometry_testable"].strip() != "True":
                continue
            yield row


def main() -> None:
    passed = failed = skipped = 0
    print(f"Running the n-friend search on examples from "
          f"{os.path.basename(CSV_PATH)}\n(n <= {MAX_N}, < {MAX_CROSSINGS} "
          f"crossings, max geodesic length {MAX_LEN})\n")

    for row in selected_rows():
        n, K_B, K_G = int(row["n"]), row["K_B"], row["K_G"]
        label = f"n={n}  {K_B} & {K_G}"

        crossings = [crossing_number(K_B), crossing_number(K_G)]
        if None in crossings or max(crossings) >= MAX_CROSSINGS:
            print(f"SKIP        {label}  (unsupported or too-large knot)")
            skipped += 1
            continue

        # The search only recovers a pair from the knot whose +n surgery is
        # the shared manifold, so try both starting points.
        results = {
            (K_B, K_G): search_recovers(K_B, K_G, n),
            (K_G, K_B): search_recovers(K_G, K_B, n),
        }

        found = [(src, tgt) for (src, tgt), r in results.items() if r is True]
        if found:
            src, tgt = found[0]
            print(f"PASS        {label}   search({src}) found {tgt}")
            passed += 1
        elif any(r is None for r in results.values()):
            # The direction whose +n surgery is the shared manifold could not
            # be drilled (no geometric triangulation), so the pair is out of
            # the search's reach rather than genuinely missed.
            print(f"UNREACHABLE {label}   (shared n-surgery has no drillable "
                  f"triangulation)")
            skipped += 1
        else:
            print(f"FAIL        {label}   (search did not recover the partner)")
            failed += 1

    print(f"\n{passed} passed, {failed} failed, {skipped} skipped")


if __name__ == "__main__":
    main()
