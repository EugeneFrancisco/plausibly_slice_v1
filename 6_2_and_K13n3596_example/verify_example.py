"""Verify Qin's m=0, r=-1 example and test the n-RBG search."""

import os
import sys

import snappy
from plink import LinkManager
from sage.all import matrix

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "code"))

from diagram_info import BLUE_KNOT, FRAMINGS, GREEN_KNOT, N  # noqa: E402
from n_rbg import NRedBlueGreenLink, find_n_special_rbg_link  # noqa: E402
from rbg import is_hopf_link, is_unlink  # noqa: E402


def load_link():
    """Load the SnapPy Link Editor projection without opening the GUI."""
    manager = LinkManager()
    with open(os.path.join(HERE, "diagram.lnk")) as stream:
        manager._from_string(stream.read())
    return snappy.Link(manager.PD_code())


def direct_verification():
    """Check the properties asserted in Qin et al., Example 3.2."""
    link = load_link()
    rbg = NRedBlueGreenLink(link, N, FRAMINGS)
    linking = matrix(link.linking_matrix())
    framed_linking = matrix(linking)
    framed_linking[0, 0] = FRAMINGS[0][0]

    checks = {
        "three components": len(link.link_components) == 3,
        "all components are unknotted": all(
            is_unlink(link.sublink([i])) for i in range(3)
        ),
        "red-blue is a Hopf link": is_hopf_link(link.sublink([0, 1])),
        "red-green is a Hopf link": is_hopf_link(link.sublink([0, 2])),
        "blue and green framings are zero": FRAMINGS[1:] == [(0, 1), (0, 1)],
        "red framing is -1": FRAMINGS[0] == (-1, 1),
        "|det(M_L)| = 1": abs(int(framed_linking.det())) == N,
        "NRedBlueGreenLink recognizes 1-special": rbg.is_n_special(),
        "K_B is 6_2": rbg.blue_exterior.is_isometric_to(
            snappy.Manifold(BLUE_KNOT)
        ),
        "K_G is K13n3596": rbg.green_exterior.is_isometric_to(
            snappy.Manifold(GREEN_KNOT)
        ),
    }

    blue_surgery = snappy.Manifold(BLUE_KNOT)
    blue_surgery.dehn_fill((N, 1))
    green_surgery = snappy.Manifold(GREEN_KNOT)
    green_surgery.dehn_fill((N, 1))
    checks["the 1-surgeries are isometric"] = blue_surgery.is_isometric_to(
        green_surgery
    )

    print("Direct verification")
    print("  linking matrix:", linking)
    print("  framed determinant:", framed_linking.det())
    for description, passed in checks.items():
        print(f"  {'PASS' if passed else 'FAIL'}: {description}")
    return rbg, all(checks.values())


def search_verification(known_rbg):
    """Check whether the public search rediscovers Qin's example."""
    print("\nSearch verification")
    found = find_n_special_rbg_link(
        snappy.Manifold(BLUE_KNOT), snappy.Manifold(GREEN_KNOT), N
    )
    if found is None:
        print("  FAIL: find_n_special_rbg_link returned None")
        return False

    checks = {
        "found link is 1-special": found.is_n_special(),
        "found K_B is 6_2": found.blue_exterior.is_isometric_to(
            snappy.Manifold(BLUE_KNOT)
        ),
        "found K_G is K13n3596": found.green_exterior.is_isometric_to(
            snappy.Manifold(GREEN_KNOT)
        ),
        "found exterior matches saved diagram": found.exterior.is_isometric_to(
            known_rbg.exterior
        ),
    }
    for description, passed in checks.items():
        print(f"  {'PASS' if passed else 'FAIL'}: {description}")
    return all(checks.values())


if __name__ == "__main__":
    known, direct_passed = direct_verification()
    search_passed = search_verification(known)
    if not direct_passed or not search_passed:
        raise SystemExit(1)
