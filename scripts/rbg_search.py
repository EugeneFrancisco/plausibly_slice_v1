"""Demonstrate the RBG-link workflow for a pair of n-friends.

Run in the ``sage`` conda environment. The s-invariant calculation also
requires KnotJob and Java 23 or newer.

Examples:

    python code/rbg_search.py 6_2 K13n3596 1
    python code/rbg_search.py K11n34 <friend-isosig> 0
"""

from __future__ import annotations

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "code"))

from knot import Knot  # noqa: E402


def find_rbg_link(knot_b: Knot, knot_g: Knot, n: int):
    """Find an RBG link encoding the given n-friends."""
    if n == 0:
        return Knot.search_for_rbg_link(knot_b, knot_g)
    return Knot.search_for_n_rbg_link(knot_b, knot_g, n)


def is_special(rbg, n: int) -> bool:
    """Check the relevant special-link condition."""
    if n == 0:
        return rbg.is_super_special()
    return rbg.is_n_special()


def traces_are_diffeomorphic(rbg, n: int) -> bool:
    """Check whether the special-link results certify diffeomorphic traces."""
    return is_special(rbg, n) and rbg.framings[0] == (0, 1)


def encoded_knots(rbg) -> tuple[Knot, Knot]:
    """Build the knots encoded by the blue and green components."""
    knot_b = Knot(exterior=rbg.blue_exterior, name="K_B")
    knot_g = Knot(exterior=rbg.green_exterior, name="K_G")
    return knot_b, knot_g


def analyze(knot_b_name: str, knot_g_name: str, n: int) -> None:
    """Run the RBG search and its trace and s-invariant checks."""
    knot_b = Knot(knot_b_name)
    knot_g = Knot(knot_g_name)

    if not knot_b.is_n_friend(knot_g, n):
        print(f"{knot_b_name} and {knot_g_name} are not verified {n}-friends.")
        return

    print(f"Verified that {knot_b_name} and {knot_g_name} are {n}-friends.")
    print("\n[1] Searching for an RBG link...")
    rbg = find_rbg_link(knot_b, knot_g, n)
    if rbg is None:
        print("No RBG link was found within the search bounds.")
        return

    print(f"Found {rbg}")
    print("\n[2] Checking the link...")
    special = is_special(rbg, n)
    red_framing = rbg.framings[0][0]
    label = "super-special" if n == 0 else f"{n}-special"
    print(f"{label}: {special}")
    print(f"red framing r: {red_framing}")

    if not traces_are_diffeomorphic(rbg, n):
        print("The special-link results used here do not certify diffeomorphic traces.")
        return

    print(f"The {n}-traces of K_B and K_G are diffeomorphic.")
    print("\n[3] Computing s-invariants of K_B...")
    encoded_b, _ = encoded_knots(rbg)
    s_invariants = encoded_b.rasmussen_s()
    print(f"s(K_B) = {s_invariants}")

    witness = next(
        ((prime, value) for prime, value in s_invariants.items() if value != 0),
        None,
    )
    if witness is None:
        print("The computed s-invariants vanish, so they give no slice conclusion.")
        return

    prime, value = witness
    print(f"s over F_{prime} is {value}, so K_B is not smoothly slice.")
    if n == 0:
        print("Because the 0-traces are diffeomorphic, K_G is also not smoothly slice.")
    else:
        print(
            "For n > 0, diffeomorphic n-traces do not by themselves transfer "
            "ordinary smooth sliceness to K_G. Qin's n-sliceness bounds are "
            "needed for that positive-n analysis."
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Search for an RBG link associated to two n-friends."
    )
    parser.add_argument("knot_b", help="SnapPy name or isosig for K_B")
    parser.add_argument("knot_g", help="SnapPy name or isosig for K_G")
    parser.add_argument("n", type=int, help="nonnegative surgery coefficient")
    args = parser.parse_args()
    if args.n < 0:
        parser.error("n must be nonnegative")
    return args


if __name__ == "__main__":
    arguments = parse_args()
    analyze(arguments.knot_b, arguments.knot_g, arguments.n)
