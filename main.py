from src.knot import Knot
from tqdm import tqdm

knot = Knot("K11n34")
friends = knot.find_n_friends(1)

count = 0
for friend in tqdm(friends):
    identifier = friend[3]
    knot_friend = Knot(identifier)
    link = Knot.search_for_n_rbg_link(knot, knot_friend, 1)
    if link is not None:
        print("yay")
