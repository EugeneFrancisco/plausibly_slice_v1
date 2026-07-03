"""
This file defines utility to accomplish three tasks:

    1. Find a pair of 0-friends  K, K'  (knots with the same 0-surgery).
    2. Given a pair of 0-friends, find out if they make a super-special RBG link.
    3. Given a super-special RBG link, find what the two knots K_G and K_B are. Then,
        compute their s-invariants. If s(K_B) != 0 or s(K_G) != 0, then by Theorem 5.9
        both K_G and K_B are not smoothly slice.

This file walks through those three parts on a concrete example, the
Conway knot K11n34. To change this, modify what knot is pointed to in
``data/unknown_with_0-friend_final.csv``.

To compute the Rasmussen s-invariant, we use knotjob. This requires Java23 to be
installed.

To run:

    conda activate sage
    python theorem_5_10_replication.py

Note that the path to KnotJob's Java is auto-detected (see _find_java).
You can override this with the KNOTJOB_JAVA environment variable if needed.
"""

from __future__ import annotations

import os
import re
import sys
import ast
import csv
import shutil
import tempfile
import subprocess

import snappy

# The paper's code lives in ``code/``; make it importable before importing from it.
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "code"))

from rbg import RedBlueGreenLink                       # noqa: E402
from find_0_friends import find_common_zero_surgery_via_words  # noqa: E402

DATA_CSV = os.path.join(_HERE, "data", "unknown_with_0-friend_final.csv")
KNOTJOB_JAR = os.path.join(_HERE, "tools", "knotjob", "KnotJob.jar")


#  ====================== Part 1:  Pairs of 0-friends ======================
# Two knots are "0-friends" when their 0-surgeries are the same.
# When that 0-surgery is hyperbolic we can check a match of the surfaces using
# SnapPy's isometry check.

def zero_surgery(knot: str | snappy.Manifold) -> snappy.Manifold:
    """
    Return the 0-surgery of a knot exterior as a closed manifold.

    Args:
        knot: A knot name SnapPy understands, or a Manifold.

    Returns:
        The 0-surgery as a closed SnapPy Manifold.
    """
    M = snappy.Manifold(knot) if isinstance(knot, str) else knot.copy()
    M.dehn_fill((0, 1))
    return M


def verify_zero_friends(
    knot_a: str | snappy.Manifold,
    knot_b: str | snappy.Manifold,
    tries: int = 5,
) -> bool:
    """
    Certify whether two knots are 0-friends (share a 0-surgery).

    Args:
        knot_a: A knot name SnapPy understands (e.g. 'K12n309') or a
            knot exterior as a Manifold (e.g. an RBG blue/green exterior).
        knot_b: A knot name SnapPy understands, or a knot exterior.
        tries: Number of randomized retries for the isometry check.
            Each randomized retry quadruples the number of triangles in
            the manifold representation.

    Returns:
        True if the two 0-surgeries are certified isometric, else False.
    """
    A, B = zero_surgery(knot_a), zero_surgery(knot_b)
    for _ in range(tries):
        try:
            return bool(A.is_isometric_to(B))
        except RuntimeError:
            A.randomize(); B.randomize()
    return False


def search_for_zero_friends(
    knot: str | snappy.Manifold, max_len: float = 3.0
) -> list[tuple[str, complex, float, str]] | None:
    """
    Search for 0-friends of a knot via short closed geodesics.

    Thin wrapper around the paper's main search routine: it drills the
    short geodesics of the hyperbolic 0-surgery and keeps those whose
    result is again a knot exterior (a fresh 0-friend).

    Args:
        knot: A knot name SnapPy understands, or a Manifold.
        max_len: Maximum geodesic length to search; larger is slower.

    Returns:
        A list of 4-tuples (word, complex_length, volume, isosig), or None.
    """
    return find_common_zero_surgery_via_words(knot, max_len)


# ============ Part 2:  Pair of 0-friends  -->  RBG link,  check super-special ============
#
# An RBG link R u B u G with framings (r, b, g) encodes a pair of
# 0-friends K_B, K_G.  It is super-special (Definition 5.7) when
# R u B and R u G are both Hopf links and b = g = 0 (so r is an integer).
#
# Instead of searching for such an RBG link from scratch, the code takes
# advantages of precomputed work from data/unknown_with_0-friend_final.csv
# to get the DT code and the framing for any base knot. Given a base knot,
# the code looks at data/unknown_with_0-friend_final.csv which gives the
# DT code and framing for an RBG link that uses the base knot and its
# 0 friend. Note that you only need the base knot for this.
#
# Then, the RBG link is verified to be super-special.

def build_rbg_link(
    dt_or_pd: str | list, framing: list[tuple[int, int]]
) -> RedBlueGreenLink:
    """
    Build a RedBlueGreenLink from a link description and framings.

    Args:
        dt_or_pd: A DT string like 'DT[vcde...]' or a PD code list.
        framing: List of three slopes, e.g. [(0,1),(0,1),(0,1)].

    Returns:
        The constructed RedBlueGreenLink.
    """
    link = snappy.Link(dt_or_pd)
    return RedBlueGreenLink(link, framing)


def check_super_special(rbg: RedBlueGreenLink) -> bool:
    """
    Check whether an RBG link is super-special (Def. 5.7).

    Args:
        rbg: The RedBlueGreenLink to test.

    Returns:
        True iff ``rbg`` is super-special.
    """
    return rbg.is_super_special()


# =========== Part 3:  Super-special RBG link  -->  K_B, K_G  and their s-invariants =======
#
# From a super-special L we recover the two knots K_B (blue) and K_G
# (green) as honest diagrams in S^3, then compute their Rasmussen
# s-invariants.  Theorem 5.9: if s_F(K_G) != 0 for some field F, then
# BOTH K_B and K_G are not smoothly slice.

def extract_knots(
    rbg: RedBlueGreenLink, simplify: bool = True
) -> tuple[snappy.Link, snappy.Link]:
    """
    Extract the two knots K_B and K_G encoded by an RBG link.

    Args:
        rbg: The RedBlueGreenLink to extract from.
        simplify: If True, simplify each diagram before returning.

    Returns:
        A tuple (blue_link, green_link) as SnapPy Link diagrams.
    """
    blue, green = rbg.blue_knot(), rbg.green_knot()
    if simplify:
        blue, green = blue.copy(), green.copy()
        blue.simplify("global")
        green.simplify("global")
    return blue, green


def _find_java() -> str:
    """
    Locate a Java >= 23 capable of running KnotJob.jar.

    Returns:
        Path to a suitable java executable.
    """
    if os.environ.get("KNOTJOB_JAVA"):
        return os.environ["KNOTJOB_JAVA"]
    # Prefer a conda env named jdk23 next to the active env.
    for base in [os.path.expanduser("~/miniforge3/envs/jdk23/bin/java"),
                 os.path.expanduser("~/miniconda3/envs/jdk23/bin/java"),
                 os.path.expanduser("~/anaconda3/envs/jdk23/bin/java")]:
        if os.path.exists(base):
            return base
    found = shutil.which("java")
    if found:
        return found
    raise RuntimeError("No java found; set KNOTJOB_JAVA to a Java >= 23.")


def _pd_to_knotjob(link: snappy.Link, name: str) -> str:
    """
    Format a SnapPy Link as a KnotJob planar-diagram line.

    Args:
        link: The SnapPy Link to format.
        name: Name to label the diagram with.

    Returns:
        A KnotJob 'PD' input line.
    """
    crossings = ", ".join(
        "X[" + ",".join(str(x) for x in c) + "]" for c in link.PD_code()
    )
    return f"{name} = PD {crossings}"


def rasmussen_s(
    link: snappy.Link,
    primes: tuple[int, ...] = (2, 3),
    name: str = "knot",
    max_heap: str = "8g",
) -> dict[int, int]:
    """
    Compute Rasmussen s-invariants of a knot via KnotJob.

    Args:
        link: The knot as a SnapPy Link.
        primes: Fields F_p to compute over (p = 0 means the rationals).
        name: Name to label the diagram with.
        max_heap: Max Java heap size, e.g. '8g'.

    Returns:
        A dict {p: s_value} of s-invariants keyed by prime.
    """
    java = _find_java()
    if not os.path.exists(KNOTJOB_JAR):
        raise RuntimeError(f"KnotJob.jar not found at {KNOTJOB_JAR}")

    with tempfile.TemporaryDirectory() as tmp:
        txt = os.path.join(tmp, "knot.txt")
        with open(txt, "w") as f:
            f.write(_pd_to_knotjob(link, name) + "\n")

        cmd = [java, f"-Xmx{max_heap}", "-jar", KNOTJOB_JAR, txt]
        cmd += [f"-s{p}" for p in primes]
        cmd += ["-nf"]  # print to terminal, don't write a results file
        out = subprocess.run(cmd, capture_output=True, text=True).stdout

    # Lines look like:  "S-Invariant mod 3 : 2"
    result = {}
    for p, val in re.findall(r"S-Invariant mod (\d+)\s*:\s*(-?\d+)", out):
        result[int(p)] = int(val)
    return result


def conclude_via_theorem_5_9(
    s_blue: dict[int, int], s_green: dict[int, int]
) -> tuple[bool, str | None]:
    """
    Apply Theorem 5.9: a nonzero s obstructs sliceness of both knots.

    Args:
        s_blue: s-invariants of K_B as a dict {p: s}.
        s_green: s-invariants of K_G as a dict {p: s}.

    Returns:
        A tuple (is_not_slice, witness_string_or_None).
    """
    for label, s in (("K_B", s_blue), ("K_G", s_green)):
        for p, val in s.items():
            if val != 0:
                field = "Q" if p == 0 else f"F_{p}"
                return True, f"s_{{{field}}}({label}) = {val} != 0"
    return False, None


# =========== Putting the three parts together ============

def _load_row(base_knot: str) -> dict[str, str]:
    """
    Fetch the RBG record for a base knot from the paper's CSV.

    Args:
        base_knot: The base knot name to look up.

    Returns:
        The matching CSV row as a dict.
    """
    with open(DATA_CSV) as f:
        for row in csv.DictReader(f):
            if row["base_knot"] == base_knot:
                return row
    raise KeyError(f"{base_knot} not in {DATA_CSV}")


def demonstrate(base_knot="K11n34"):
    """
    Run Parts 1-3 end to end and print the Theorem 5.10 conclusion.

    Args:
        base_knot: The base knot from Table 10 to demonstrate.

    Returns:
        True if the knot is shown to be not smoothly slice.
    """
    row = _load_row(base_knot)
    framing = ast.literal_eval(row["framing"])
    print("=" * 68)
    print(f"Theorem 5.10 for base knot {base_knot}")
    print("=" * 68)

    # ---- Part 2: rebuild and verify the super-special RBG link --------
    print("\n[Part 2] Rebuilding the RBG link from the recorded DT code...")
    rbg = build_rbg_link(row["RBG_DT"], framing)
    print(f"    {rbg}")
    ss = check_super_special(rbg)
    print(f"    super-special? {ss}")
    assert ss, "recorded RBG link is not super-special!"

    # ---- Part 1: the two knots it encodes really are 0-friends --------
    print("\n[Part 1] The blue/green fillings share a 0-surgery (0-friends):")
    friends = verify_zero_friends(rbg.blue_exterior, rbg.green_exterior)
    print(f"    K_B and K_G are 0-friends? {friends}")
    assert friends, "blue/green fillings are not certified 0-friends!"

    # ---- Part 3: extract K_B, K_G and compute their s-invariants ------
    print("\n[Part 3] Extracting K_B, K_G and computing Rasmussen s-invariants...")
    blue, green = extract_knots(rbg)
    print(f"    K_B: {len(blue.PD_code())} crossings   "
          f"K_G: {len(green.PD_code())} crossings")
    s_blue = rasmussen_s(blue, name="K_B")
    s_green = rasmussen_s(green, name="K_G")
    print(f"    s(K_B) = {s_blue}")
    print(f"    s(K_G) = {s_green}")

    # ---- Conclusion via Theorem 5.9 -----------------------------------
    not_slice, witness = conclude_via_theorem_5_9(s_blue, s_green)
    print("\n[Conclusion]")
    if not_slice:
        print(f"    {witness}")
        print("    => by Theorem 5.9, BOTH K_B and K_G are NOT smoothly slice.")
        print(f"    => {base_knot} is not smoothly slice.  [Theorem 5.10]")
    else:
        print("    All computed s-invariants vanish; this field does not")
        print("    obstruct sliceness (try more primes / the Sq^1 refinements).")
    return not_slice


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else "K11n34"
    demonstrate(target)
