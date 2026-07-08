"""
Verify the documented n-friend examples from ``qin_n_friends.csv``.

Two knots are n-friends when their integer n-surgeries are the same closed
3-manifold (the surgery generalization of the 0-friends the paper uses).
This script reads the CSV, keeps the examples we can check quickly -- n <= 3,
both knots loadable in SnapPy with fewer than 18 crossings, and flagged
``isometry_testable`` -- and confirms each pair's n-surgeries are isometric.

Surgery signs: SnapPy stores each knot with a fixed chirality that need not
match the paper's, so a pair can relate the +n surgery of one knot to the -n
surgery of the other.  We therefore accept a match at either sign and report
the slopes that matched.

Run in the ``sage`` conda env:

    conda activate sage
    python test.py
"""

from __future__ import annotations

import os
import re
import csv
import sys

import snappy

sys.setrecursionlimit(50000)

_HERE = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(_HERE, "qin_n_friends.csv")

MAX_N = 3
MAX_CROSSINGS = 18


def crossing_number(name: str) -> int | None:
    """Crossing number parsed from a knot name (e.g. K14n10164 -> 14).

    Returns None for exotic notations we cannot load, such as the pretzel
    ``K(-2,1,2)`` or the connect sum ``T(-2,3)#T(2,5)``.
    """
    match = re.match(r"K?(\d+)", name)
    return int(match.group(1)) if match else None


def load_exterior(name: str) -> snappy.Manifold | None:
    """The knot exterior as a high-precision Manifold, or None if unknown."""
    try:
        return snappy.ManifoldHP(snappy.Manifold(name))
    except Exception:
        return None


def surgery_slopes(exterior: snappy.Manifold, n: int) -> list[tuple[int, int]]:
    """The +n and -n integer surgery slopes (meridian is (1, 0))."""
    a, b = exterior.homological_longitude()
    return [(n + a, b), (a - n, b)]


def _isometric(A: snappy.Manifold, B: snappy.Manifold, tries: int = 8) -> bool:
    """True if the two closed surgeries are isometric (randomizing on error)."""
    if abs(float(A.volume()) - float(B.volume())) > 1e-6:
        return False
    for _ in range(tries):
        try:
            return bool(A.is_isometric_to(B))
        except RuntimeError:
            A.randomize(); B.randomize()
    return False


def matching_slopes(A: snappy.Manifold, B: snappy.Manifold, n: int):
    """Return an isometric (slope_A, slope_B) pair of n-surgeries, or None."""
    for slope_A in surgery_slopes(A, n):
        for slope_B in surgery_slopes(B, n):
            Af, Bf = A.copy(), B.copy()
            Af.dehn_fill(slope_A)
            Bf.dehn_fill(slope_B)
            if _isometric(Af, Bf):
                return slope_A, slope_B
    return None


def selected_rows():
    """The CSV rows worth checking: isometry-testable, small n, small knots."""
    with open(CSV_PATH) as f:
        for row in csv.DictReader(f):
            if int(row["n"]) > MAX_N:
                continue
            if row["isometry_testable"].strip() != "True":
                continue
            yield row


def main() -> None:
    passed = failed = skipped = 0
    print(f"Checking n-friend examples from {os.path.basename(CSV_PATH)} "
          f"(n <= {MAX_N}, < {MAX_CROSSINGS} crossings)\n")

    for row in selected_rows():
        n, K_B, K_G = int(row["n"]), row["K_B"], row["K_G"]
        label = f"n={n}  {K_B} & {K_G}"

        crossings = [crossing_number(K_B), crossing_number(K_G)]
        if None in crossings:
            print(f"SKIP  {label}  (unsupported knot notation)")
            skipped += 1
            continue
        if max(crossings) >= MAX_CROSSINGS:
            print(f"SKIP  {label}  (>= {MAX_CROSSINGS} crossings)")
            skipped += 1
            continue

        A, B = load_exterior(K_B), load_exterior(K_G)
        if A is None or B is None:
            print(f"SKIP  {label}  (could not load in SnapPy)")
            skipped += 1
            continue

        match = matching_slopes(A, B, n)
        if match:
            slope_A, slope_B = match
            print(f"PASS  {label}   {K_B}^{slope_A} == {K_G}^{slope_B}")
            passed += 1
        else:
            print(f"FAIL  {label}   (n-surgeries not isometric)")
            failed += 1

    print(f"\n{passed} passed, {failed} failed, {skipped} skipped")


if __name__ == "__main__":
    main()
