"""
Reproduce Theorem 5.10 of Dunfield-Gong for the 25 knots of Table 10.

For each base knot K in Table 10 we obtain a super-special RBG link relating
K to a 0-friend K', recover the two simpler knots K_B, K_G the link encodes,
compute their Rasmussen s-invariants, and apply Theorem 5.9 to conclude that
K (and its 0-friend) are *not* smoothly slice.

Two ways to get the RBG link, controlled by ``--search``:

  * recorded (default): rebuild it from the DT code in the CSV. Instead of searching for
  RBG links, we just used the ones given to us in data/unknown_with_0-friend_final.csv
  * live (``--search``): run the Section 5.11 search from the base knot's
    exterior and its 0-friend's exterior (CSV ``tri`` isosig).  This will actually
    search for the RBG links that can prove non-sliceness for the base knots. But it is
    much slower.

Regardless of how we find RBG linkes, we always (a) verify it is super-special (b) extract
K_B, K_G (c) confirm 0-friends (d) get the s-invariants (e) hopefully invoke Theorem 5.9.

Run in the ``sage`` conda env; the s-invariant step needs KnotJob + Java >= 23.

Usage:

    conda activate sage
    python main.py                    # all 25, recorded links (fast)
    python main.py K11n34             # one knot, recorded link
    python main.py --search K11n34    # one knot, live Section 5.11 search
    python main.py --search           # all 25 via live search (very slow)
"""

from __future__ import annotations

import csv
import os
import sys
from datetime import datetime

import snappy

from knot import Knot, DATA_CSV

# Verdicts are appended here one knot at a time so progress survives a long
# (or interrupted) run; see ``main``.
RESULTS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "table_10_results.txt")


# The 25 base knots of Table 10, in the paper's order.
TABLE_10 = [
    "K11n34",
    "K13n866",
    "K15n25044",
    "16n180537",
    "16n74539",
    "17nh_0001844",
    "17nh_0002715",
    "17nh_0212094",
    "17nh_0212095",
    "18nh_00010270",
    "18nh_00166702",
    "18nh_00610378",
    "18nh_00610381",
    "19nh_000002588",
    "19nh_000003154",
    "19nh_000003570",
    "19nh_000032808",
    "19nh_000076489",
    "19nh_000018991",
    "19nh_000066839",
    "19nh_000066841",
    "19nh_000177115",
    "19nh_000177116",
    "19nh_001336127",
    "19nh_002457201",
]

# The lone r=0 case (Theorem 5.8): concluded via the Sq^1-refined s-invariant,
# not s_2/s_3, so we expect s_2/s_3 to vanish here.  See Table 10 / Thm 5.8.
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

    # Section 5.11 search from the two 0-friend exteriors: E_K from the base
    # knot's name (*nh_* knots need PlausibleKnots), E_K' from the CSV isosig.
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
    say("\n[Part 2] RBG link:")
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
    k_b, k_g = Knot.blue_green_knots(rbg)
    friends = k_b.is_zero_friend(k_g)
    say(f"\n[Part 1] K_B ({len(k_b.diagram().PD_code())} cr) and "
        f"K_G ({len(k_g.diagram().PD_code())} cr) are 0-friends? {friends}")

    # ---- Part 3: s-invariants and Theorem 5.9 -------------------------
    say("\n[Part 3] Rasmussen s-invariants (KnotJob):")
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
        result["note"] = "r=0 case: needs Sq^1-refined s (Thm 5.8), not s_2/s_3"
        say(f"    s_2/s_3 vanish as expected ({result['note']}).")
    else:
        result["note"] = "all computed s-invariants vanish; try more primes"
        say(f"    {result['note']}")
    return result


def nonzero_s(r: dict) -> str:
    """The first nonzero s-invariant among K_B, K_G as e.g. 's_{F_3} = 2'."""
    for s in (r["s_b"], r["s_g"]):
        for p, v in s.items():
            if v != 0:
                field = "Q" if p == 0 else f"F_{p}"
                return f"s_{{{field}}} = {v}"
    return ""


def result_line(r: dict) -> str:
    """One-line verdict for a knot, used in both the file and the summary."""
    if r["not_slice"]:
        return f"{r['base_knot']:<16} not slice   {nonzero_s(r)}"
    return f"{r['base_knot']:<16} --          {r['note']}"


def main(argv: list[str]) -> None:
    """
    Run the pipeline over one or all Table 10 knots and print a summary.

    Each knot's verdict is appended to ``RESULTS_FILE`` as soon as it is
    computed, so partial progress is visible during a long run.

    Args:
        argv: Command-line arguments after the program name.  A ``--search``
            flag selects the live search; a remaining argument names a single
            base knot (default: all of Table 10).
    """
    live = "--search" in argv
    args = [a for a in argv if a != "--search"]
    knots = [args[0]] if args else TABLE_10
    mode = "live search" if live else "recorded links"

    with open(RESULTS_FILE, "a") as f:
        f.write(f"\n# Theorem 5.10 verdicts ({mode}), "
                f"started {datetime.now():%Y-%m-%d %H:%M:%S}\n")
        f.write(f"# {len(knots)} knot(s); each line is written as it finishes.\n")

    results = []
    for base_knot in knots:
        try:
            r = test_one_knot(base_knot, live=live)
        except Exception as e:
            print(f"\n[error] {base_knot}: {type(e).__name__}: {e}", flush=True)
            r = {"base_knot": base_knot, "not_slice": False, "s_b": {},
                 "s_g": {}, "witness": None, "note": f"{type(e).__name__}: {e}"}
        results.append(r)
        with open(RESULTS_FILE, "a") as f:
            f.write(result_line(r) + "\n")

    # ---- summary ------------------------------------------------------
    shown = sum(r["not_slice"] for r in results)
    summary = [
        "=" * 70,
        f"SUMMARY  ({mode})",
        "=" * 70,
        *(f"  {result_line(r)}" for r in results),
        "-" * 70,
        f"  Shown not smoothly slice via Theorem 5.9: {shown}/{len(results)}",
    ]
    if any(r["base_knot"] in NEEDS_SQ1 for r in results):
        summary.append("  (18nh_00010270 is the separate Theorem 5.8 / Sq^1 case.)")
    print("\n" + "\n".join(summary))
    with open(RESULTS_FILE, "a") as f:
        f.write("\n" + "\n".join(summary) + "\n")
    print(f"\n(Verdicts written to {RESULTS_FILE})")


if __name__ == "__main__":
    main(sys.argv[1:])
