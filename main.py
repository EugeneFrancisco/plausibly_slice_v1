import snappy
from plink import LinkManager

from src.knot import Knot
from src.n_rbg import find_n_special_rbg_link
from tqdm import tqdm

FULL_SEARCH = False

if FULL_SEARCH:

    knot = Knot("6_2")
    friends = knot.find_n_friends(1)

    count = 0
    for friend in tqdm(friends):
        identifier = friend[3]
        knot_friend = Knot(identifier)
        link = Knot.search_for_n_rbg_link(knot, knot_friend, 1)
        if link is not None:
            print("yay")
else:
    manager = LinkManager()
    with open("6_2_and_K13n3596_example/diagram.lnk") as stream:
        manager._from_string(stream.read())
    expected = snappy.Link(manager.PD_code()).exterior()

    link = find_n_special_rbg_link(
        snappy.Manifold("6_2"), snappy.Manifold("K13n3596"), 1,
        target_exterior=expected,
    )
    if link is not None and link.exterior.is_isometric_to(expected):
        print("Recovered the saved 6_2 / K13n3596 RBG diagram.")
    else:
        print("Wasn't able to.")