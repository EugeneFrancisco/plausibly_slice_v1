"""
Reproduce Theorem 5.10 of Dunfield-Gong for the 25 knots of Table 10.

For each base knot K in Table 10 we obtain a super-special RBG link relating
K to a 0-friend K', recover the two simpler knots K_B, K_G the link encodes,
compute their Rasmussen s-invariants, and apply Theorem 5.9 to conclude that
K (and its 0-friend) are *not* smoothly slice.

There are two ways to get the RBG link, controlled by ``--search``:

  * recorded (default):  rebuild the link from the DT code recorded against
    the base knot in ``data/unknown_with_0-friend_final.csv``.  Fast
    (seconds/knot) and works for all 25 knots, including the ``*nh_*``
    census knots whose names SnapPy cannot resolve on its own.

  * live (``--search``):  actually run the Section 5.11 search from the two
    0-friend exteriors -- the base knot's exterior (from its name) and its
    0-friend's exterior (from the ``tri`` isosig in the CSV).  This is the
    faithful "search for the RBG link" path but is expensive (minutes to
    hours per knot) and, without the PlausibleKnots census installed, only
    works for base knots SnapPy can name (K11n34, K13n866, K15n25044,
    16n180537, 16n74539).

Both paths then run the identical downstream pipeline (verify super-special
-> extract K_B, K_G -> confirm 0-friends -> s-invariants -> Theorem 5.9).

Environment: the ``sage`` conda env (SnapPy + Sage); the s-invariant step
shells out to KnotJob and needs Java >= 23 (auto-detected, see knot._find_java).

Usage:

    conda activate sage
    python main.py                    # all 25, recorded links (fast)
    python main.py K11n34             # one knot, recorded link
    python main.py --search K11n34    # one knot, live Section 5.11 search
    python main.py --search           # all 25 via live search (very slow)
"""

from __future__ import annotations

import csv
import sys

import snappy

from knot import Knot, DATA_CSV


# The 25 knots of Table 10, in the paper's order.  Each is a base knot
# K in PS_19 for which a super-special RBG link to a 0-friend K' was found.
TABLE_10 = [
    "K11n34", "K13n866", "K15n25044", "16n180537", "16n74539",
    "17nh_0001844", "17nh_0002715", "17nh_0212094", "17nh_0212095",
    "18nh_00010270", "18nh_00166702", "18nh_00610378", "18nh_00610381",
    "19nh_000002588", "19nh_000003154", "19nh_000003570", "19nh_000032808",
    "19nh_000076489", "19nh_000018991", "19nh_000066839", "19nh_000066841",
    "19nh_000177115", "19nh_000177116", "19nh_001336127", "19nh_002457201",
]

# The lone Theorem 5.8 (r=0 diffeo-trace) case.  Its RBG link has r=0, so K
# and K' have diffeomorphic traces and the conclusion runs through the Sq^1
# refinement s_{Sq^1_o(K)} != 0 rather than the plain s_2/s_3 that Theorem 5.9
# and KnotJob's -s2/-s3 give us.  We flag it so a vanishing s_2/s_3 here is
# reported as "expected, needs Sq^1" instead of a spurious failure.
NEEDS_SQ1 = {"18nh_00010270"}


def partner_exterior(base_knot: str) -> snappy.Manifold:
    """
    The 0-friend K' of ``base_knot`` as a one-cusped SnapPy exterior.

    Read from the ``tri`` triangulation isosig recorded against the base
    knot in the paper's CSV (this is the exterior that Stage A found sharing
    a 0-surgery with the base knot).

    Args:
        base_knot: A Table 10 base knot name.

    Returns:
        The 0-friend exterior as a Manifold.
    """
    with open(DATA_CSV) as f:
        for row in csv.DictReader(f):
            if row["base_knot"] == base_knot:
                return snappy.Manifold(row["tri"])
    raise KeyError(f"{base_knot} not in {DATA_CSV}")


def get_rbg_link(base_knot: str, *, live: bool, verbose: bool = True):
    """
    Obtain a super-special RBG link for ``base_knot``.

    Args:
        base_knot: A Table 10 base knot name.
        live: If True, run the Section 5.11 search from the base knot's
            exterior and its 0-friend's exterior (slow).  If False, rebuild
            the recorded link from the CSV DT code (fast).
        verbose: If True, print progress during a live search.

    Returns:
        A super-special RedBlueGreenLink, or None if a live search found none.
    """
    if not live:
        return Knot.recorded_rbg_link(base_knot)

    # Live search: the two independent 0-friend exteriors.  E_K comes from
    # the base knot's name (so *nh_* census knots need PlausibleKnots), and
    # E_K' from the recorded 0-friend triangulation.
    try:
        E_K = snappy.Manifold(base_knot)
    except Exception as e:
        raise RuntimeError(
            f"SnapPy cannot resolve base knot {base_knot!r} by name "
            f"({type(e).__name__}); the live search needs the PlausibleKnots "
            f"census for *nh_* knots. Use the recorded path for this knot."
        ) from e
    E_Kprime = partner_exterior(base_knot)
    return Knot.search_for_rbg_link(E_K, E_Kprime, verbose=verbose)


def test_one_knot(base_knot: str, *, live: bool, verbose: bool = True) -> dict:
    """
    Run the full Theorem 5.10 pipeline for one Table 10 base knot.

    Parts, mirroring the proof:
        2. Obtain and verify a super-special RBG link (recorded or searched).
        1. Confirm the blue/green fillings K_B, K_G are 0-friends.
        3. Compute Rasmussen s-invariants and apply Theorem 5.9.

    Args:
        base_knot: A Table 10 base knot name.
        live: Whether to search for the link live (see ``get_rbg_link``).
        verbose: If True, print the step-by-step narrative.

    Returns:
        A result dict with keys: base_knot, not_slice (bool), s_b, s_g
        (s-invariant dicts), witness (str or None), note (str).
    """
    result = {"base_knot": base_knot, "not_slice": False,
              "s_b": {}, "s_g": {}, "witness": None, "note": ""}

    def say(*a):
        if verbose:
            print(*a, flush=True)

    say("=" * 70)
    say(f"Theorem 5.10 for base knot {base_knot}"
        f"   ({'live search' if live else 'recorded link'})")
    say("=" * 70)

    # ---- Part 2: obtain and verify the super-special RBG link ----------
    say("\n[Part 2] Obtaining the super-special RBG link...")
    rbg = get_rbg_link(base_knot, live=live, verbose=verbose)
    if rbg is None:
        result["note"] = "no super-special RBG link found by the live search"
        say(f"    {result['note']}")
        return result
    say(f"    {rbg}")
    ss = rbg.is_super_special()
    say(f"    super-special? {ss}")
    if not ss:
        result["note"] = "RBG link is not super-special"
        say(f"    {result['note']}")
        return result

    # ---- Part 1: the two encoded knots are 0-friends ------------------
    say("\n[Part 1] Checking the blue/green fillings share a 0-surgery...")
    k_b, k_g = Knot.blue_green_knots(rbg)
    friends = k_b.is_zero_friend(k_g)
    say(f"    K_B ({len(k_b.diagram().PD_code())} crossings) and "
        f"K_G ({len(k_g.diagram().PD_code())} crossings) are 0-friends? {friends}")

    # ---- Part 3: s-invariants and Theorem 5.9 -------------------------
    say("\n[Part 3] Computing Rasmussen s-invariants (KnotJob)...")
    result["s_b"] = k_b.rasmussen_s()
    result["s_g"] = k_g.rasmussen_s()
    say(f"    s(K_B) = {result['s_b']}")
    say(f"    s(K_G) = {result['s_g']}")

    not_slice, witness = Knot.conclude_via_theorem_5_9(k_b, k_g)
    result["not_slice"] = not_slice
    result["witness"] = witness

    say("\n[Conclusion]")
    if not_slice:
        say(f"    {witness}")
        say(f"    => by Theorem 5.9, K_B and K_G are NOT smoothly slice.")
        say(f"    => {base_knot} is not smoothly slice.  [Theorem 5.10]")
    elif base_knot in NEEDS_SQ1:
        result["note"] = (
            "expected: this is the r=0 diffeo-trace case (Theorem 5.8); "
            "the conclusion needs the Sq^1-refined s-invariant, not s_2/s_3"
        )
        say(f"    s_2/s_3 vanish here as expected -- {base_knot} is the r=0")
        say(f"    diffeo-trace knot, concluded via Theorem 5.8 and the Sq^1")
        say(f"    refinement (outside the s_2/s_3 pipeline).")
    else:
        result["note"] = "all computed s-invariants vanish; try more primes"
        say(f"    {result['note']}")
    return result


def main(argv: list[str]) -> None:
    """
    Run the pipeline over one or all Table 10 knots and print a summary.

    Args:
        argv: Command-line arguments after the program name.  A ``--search``
            flag selects the live search; a remaining argument names a single
            base knot (default: all of Table 10).
    """
    live = "--search" in argv
    args = [a for a in argv if a != "--search"]
    knots = [args[0]] if args else TABLE_10

    results = []
    for base_knot in knots:
        try:
            results.append(test_one_knot(base_knot, live=live))
        except Exception as e:
            print(f"\n[error] {base_knot}: {type(e).__name__}: {e}", flush=True)
            results.append({"base_knot": base_knot, "not_slice": False,
                            "s_b": {}, "s_g": {}, "witness": None,
                            "note": f"{type(e).__name__}: {e}"})

    # ---- summary ------------------------------------------------------
    print("\n" + "=" * 70)
    print(f"SUMMARY  ({'live search' if live else 'recorded links'})")
    print("=" * 70)
    shown = 0
    for r in results:
        if r["not_slice"]:
            shown += 1
            print(f"  NOT SLICE   {r['base_knot']:<16} {r['witness']}")
        else:
            print(f"  --          {r['base_knot']:<16} {r['note']}")
    print("-" * 70)
    print(f"  Shown not smoothly slice via Theorem 5.9: {shown}/{len(results)}")
    if any(r['base_knot'] in NEEDS_SQ1 for r in results):
        print(f"  (18nh_00010270 is the separate Theorem 5.8 / Sq^1 case.)")


if __name__ == "__main__":
    main(sys.argv[1:])
