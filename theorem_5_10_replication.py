"""
Replicate Theorem 5.10 of Dunfield-Gong (25 knots that are not smoothly
slice) end to end for one base knot, in three parts mirroring the proof.

All of the machinery now lives on the ``Knot`` class in ``knot.py``; this
file is just the demonstration driver.  The three parts are:

    1. Find a pair of 0-friends  K, K'  (knots with the same 0-surgery).
    2. Turn the pair into a super-special RBG link.  Either load the
       precomputed link from ``data/unknown_with_0-friend_final.csv``
       (fast) or search for it live following Section 5.11 (slow).
    3. Recover the two knots K_B, K_G the link encodes and compute their
       Rasmussen s-invariants.  If either is nonzero then, by Theorem 5.9,
       both K_B and K_G are not smoothly slice.

Computing the Rasmussen s-invariant uses KnotJob, which needs Java >= 23.

To run:

    conda activate sage
    python theorem_5_10_replication.py [BASE_KNOT]            # fast, from CSV
    python theorem_5_10_replication.py BASE_KNOT PARTNER      # live search

The path to KnotJob's Java is auto-detected (see knot._find_java); override
it with the KNOTJOB_JAVA environment variable if needed.
"""

from __future__ import annotations

import sys

from knot import Knot


def demonstrate(base_knot):
    """
    Run Parts 1-3 end to end and print the Theorem 5.10 conclusion, using
    the precomputed RBG link from the paper's CSV.

    Args:
        base_knot: The base knot from Table 10 to demonstrate.

    Returns:
        True if the knot is shown to be not smoothly slice.
    """
    print("=" * 68)
    print(f"Theorem 5.10 for base knot {base_knot}")
    print("=" * 68)

    # ---- Part 2: rebuild and verify the super-special RBG link --------
    print("\n[Part 2] Rebuilding the RBG link from the recorded DT code...")
    rbg = Knot.recorded_rbg_link(base_knot)
    print(f"    {rbg}")
    print(f"    super-special? {rbg.is_super_special()}")
    assert rbg.is_super_special(), "recorded RBG link is not super-special!"

    # ---- Part 1: the two knots it encodes really are 0-friends --------
    print("\n[Part 1] The blue/green fillings share a 0-surgery (0-friends):")
    k_b, k_g = Knot.blue_green_knots(rbg)
    friends = k_b.is_zero_friend(k_g)
    print(f"    K_B and K_G are 0-friends? {friends}")
    assert friends, "blue/green fillings are not certified 0-friends!"

    # ---- Part 3: extract K_B, K_G and compute their s-invariants ------
    print("\n[Part 3] Extracting K_B, K_G and computing Rasmussen s-invariants...")
    print(f"    K_B: {len(k_b.diagram().PD_code())} crossings   "
          f"K_G: {len(k_g.diagram().PD_code())} crossings")
    print(f"    s(K_B) = {k_b.rasmussen_s()}")
    print(f"    s(K_G) = {k_g.rasmussen_s()}")

    # ---- Conclusion via Theorem 5.9 -----------------------------------
    not_slice, witness = Knot.conclude_via_theorem_5_9(k_b, k_g)
    print("\n[Conclusion]")
    if not_slice:
        print(f"    {witness}")
        print("    => by Theorem 5.9, BOTH K_B and K_G are NOT smoothly slice.")
        print(f"    => {base_knot} is not smoothly slice.  [Theorem 5.10]")
    else:
        print("    All computed s-invariants vanish; this field does not")
        print("    obstruct sliceness (try more primes / the Sq^1 refinements).")
    return not_slice


def demonstrate_search(base_knot, partner):
    """
    Run the full Theorem 5.10 pipeline live for a 0-friend pair.

    Unlike ``demonstrate``, which loads a precomputed RBG link from the
    CSV, this drives ``Knot.conclude_slice_status_via_rbg`` on the pair
    (base_knot, partner): it searches from scratch (Section 5.11) for a
    super-special RBG link, recovers the two knots K_B, K_G it encodes,
    computes their Rasmussen s-invariants, and applies Theorem 5.9.

    This is the expensive path -- it drills geodesics and does isometry
    checks in SnapPy -- so it is kept separate from the fast CSV-based
    ``demonstrate``.  Use it to see the search working end to end on a
    pair you already know are 0-friends.

    Args:
        base_knot: One knot of the pair (name/isosig/Manifold).
        partner: The other knot of the pair (name/isosig/Manifold).

    Returns:
        True if the pair is shown to be not smoothly slice.
    """
    K = Knot(base_knot)
    K_prime = Knot(partner)
    print("=" * 68)
    print(f"Live RBG pipeline for the 0-friend pair "
          f"({K.label}, {K_prime.label})")
    print("=" * 68)
    return K.conclude_slice_status_via_rbg(K_prime)


if __name__ == "__main__":
    # Usage:
    #   python theorem_5_10_replication.py [BASE_KNOT]
    #       Fast path: rebuild the recorded RBG link from the CSV and run
    #       the full Theorem 5.10 pipeline (Parts 1-3).
    #   python theorem_5_10_replication.py BASE_KNOT PARTNER
    #       Slow path: live Section 5.11 search showing how the RBG link
    #       is found from scratch for the 0-friend pair (BASE_KNOT, PARTNER).
    #
    # Example:
    # python theorem_5_10_replication.py K11n34 'sLLvLvzPzQQQQcdgnmpkqkpqrmnqoorrhsgpxhpxudhbnwsssko_bBba'

    if len(sys.argv) > 2:
        demonstrate_search(sys.argv[1], sys.argv[2])
    else:
        demonstrate(sys.argv[1] if len(sys.argv) > 1 else "K11n34")
