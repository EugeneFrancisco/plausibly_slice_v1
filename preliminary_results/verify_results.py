"""Reload and verify the saved preliminary-result diagrams in Sage."""

from pathlib import Path

import snappy
from plink import LinkManager
from sage.all import matrix

from rbg import is_hopf_link, is_unlink


HERE = Path(__file__).resolve().parent
FRAMINGS = [(-1, 1), (0, 1), (0, 1)]

FRIENDS = {
    "6_2_1_friend_K13n3596": ("6_2", 1),
    "6_2_3_friend_K14n10164_mirror": ("6_2", 3),
    "conway_1_friend_A": ("K11n34", 1),
    "conway_1_friend_B": ("K11n34", 1),
}

RBG_LINKS = {
    "6_2_K13n3596_1_RBG": ("6_2", "6_2_1_friend_K13n3596"),
    "conway_friend_A_1_RBG": ("K11n34", "conway_1_friend_A"),
}


def load_link(path):
    manager = LinkManager()
    manager._from_string(path.read_text())
    return snappy.Link(manager.PD_code())


def isometric(first, second):
    first, second = first.copy(), second.copy()
    for _ in range(10):
        try:
            return bool(first.is_isometric_to(second))
        except RuntimeError:
            first.randomize()
            second.randomize()
    return False


def verify_friend(name, base_name, n):
    link = load_link(HERE / "friends" / f"{name}.lnk")
    exterior = link.exterior()
    triangulation = snappy.Manifold(str(HERE / "friends" / f"{name}.tri"))
    assert isometric(exterior, triangulation)

    base_surgery = snappy.Manifold(base_name)
    base_surgery.dehn_fill((n, 1))
    friend_surgery = exterior.copy()
    friend_surgery.dehn_fill((n, 1))
    assert isometric(base_surgery, friend_surgery)
    print(f"PASS: {name} is a {n}-friend of {base_name}")


def verify_rbg(name, blue_name, green_file):
    link = load_link(HERE / "rbg_links" / f"{name}.lnk")
    exterior = link.exterior()
    triangulation = snappy.Manifold(str(HERE / "rbg_links" / f"{name}.tri"))
    assert isometric(exterior, triangulation)

    red, blue, green = FRAMINGS
    blue_exterior = exterior.copy()
    blue_exterior.dehn_fill([red, (0, 0), green])
    blue_exterior = blue_exterior.filled_triangulation()
    green_exterior = exterior.copy()
    green_exterior.dehn_fill([red, blue, (0, 0)])
    green_exterior = green_exterior.filled_triangulation()
    assert isometric(blue_exterior, snappy.Manifold(blue_name))
    assert isometric(
        green_exterior,
        snappy.Manifold(str(HERE / "friends" / f"{green_file}.tri")),
    )

    filled = exterior.copy()
    filled.dehn_fill(FRAMINGS)
    assert filled.homology().elementary_divisors() == []
    assert all(is_unlink(link.sublink([i])) for i in range(3))
    assert is_hopf_link(link.sublink([0, 1]))
    assert is_hopf_link(link.sublink([0, 2]))
    linking = matrix(link.linking_matrix())
    linking[0, 0] = -1
    assert int(linking.det()) == -1
    print(f"PASS: {name} is a 1-special RBG link for the stated pair")


def main():
    for name, (base_name, n) in FRIENDS.items():
        verify_friend(name, base_name, n)
    for name, (blue_name, green_file) in RBG_LINKS.items():
        verify_rbg(name, blue_name, green_file)


if __name__ == "__main__":
    main()
