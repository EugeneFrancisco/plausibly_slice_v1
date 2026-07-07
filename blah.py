from knot import Knot
import snappy as sn

six_two = sn.Manifold("6_2")

knot_6_2 = Knot(six_two)

friends = knot_6_2.find_n_friends(n=3, max_len=6)
print(friends)
