"""Regenerate the SnapPy files for the preliminary n-friend examples."""

from pathlib import Path

import snappy
from plink import LinkManager


HERE = Path(__file__).resolve().parent

# These isosigs are the exteriors returned by find_common_n_surgery_via_words
# with max_len=3.  The recovered diagram is mirrored when needed so that its
# positive n-surgery, rather than its negative n-surgery, matches the input.
EXAMPLES = {
    "6_2_1_friend_K13n3596": {
        "isosig": "iLvLQQcbegfhghghhutuignno_BaaB",
        "base": "6_2",
        "n": 1,
    },
    "6_2_3_friend_K14n10164_mirror": {
        "isosig": "kvLALAQkdcfihggjijjssapaubtstt_baab",
        "base": "6_2",
        "n": 3,
    },
    "conway_1_friend_A": {
        "isosig": "pLvvLAPLQQQcghjjlnomknmolonpegbbuhqlrepdewr_BbaB",
        "base": "K11n34",
        "n": 1,
    },
    "conway_1_friend_B": {
        "isosig": (
            "rvLvAvMLMQQQccifkkkoplnmqonpqpqvpxovodiflsuqaeubk_abBa"
        ),
        "base": "K11n34",
        "n": 1,
    },
}


def save_link(link, path):
    """Save a spherogram link as a SnapPy Link Editor projection."""
    manager = LinkManager()
    link.view(manager)
    path.write_text(manager.SnapPea_projection_file())


def positive_surgery_matches(link, base_name, n):
    target = snappy.ManifoldHP(base_name)
    target.dehn_fill((n, 1))
    filled = snappy.ManifoldHP(link.exterior())
    filled.dehn_fill((n, 1))
    if abs(target.volume() - filled.volume()) > 1e-20:
        return False
    for _ in range(10):
        try:
            return bool(target.is_isometric_to(filled))
        except RuntimeError:
            target.randomize()
            filled.randomize()
    return False


def recover_short_diagram(data, tries=8):
    best = None
    for _ in range(tries):
        exterior = snappy.Manifold(data["isosig"])
        link = exterior.exterior_to_link()
        link.simplify("global")
        if not positive_surgery_matches(link, data["base"], data["n"]):
            link = link.mirror()
        assert positive_surgery_matches(link, data["base"], data["n"])
        if best is None or len(link.crossings) < len(best.crossings):
            best = link
    return best


def main():
    output = HERE / "friends"
    output.mkdir(exist_ok=True)
    for name, data in EXAMPLES.items():
        link = recover_short_diagram(data)
        save_link(link, output / f"{name}.lnk")
        link.exterior().save(str(output / f"{name}.tri"))
        print(name, len(link.crossings), link.exterior().volume())


if __name__ == "__main__":
    main()
