"""

Here, we deal by with a few special alternating knots, either by
showing they are 2-bridge in a range dealt with by Casson-Goron in one
of their original papers, or by showing they are concordant to -2*7a1

"""

import snappy
import plausible_knots
from spherogram.links.bands import normalize_crossing_labels
from sage.all import is_square
import casson_gordon


alt_18s = {'18ah_3327857':'0d1403_0_0',
           '18ah_4025786':'000d0a_1_0',
           '18ah_4099296':'2c1512_1_0',
           '18ah_4099297':'391801_1_0',

}

two_bridge = {'K14a12741':(25, 2),
              '16a350194':(49, 4),
              '17ah_0000055':(121, 48),
              '18ah_0000122':(169, 48),
              '19ah_00000457':(289, 190)}


targets = ['K14a12741', '16a350194', '17ah_0000055', '18ah_0000122',
           '18ah_3327857', '18ah_4025786', '18ah_4099296', '18ah_4099297',
           '19ah_00000457']

assert set(targets) == set(alt_18s.keys()) | set(two_bridge.keys())


def orient_pres_isometric(A, B):
    assert A.num_cusps() == 1
    for iso in A.is_isometric_to(B, True):
        if iso.cusp_maps()[0].determinant() == 1:
            return True
    return False


def test_two_bridge():
    for name, desc in two_bridge.items():
        p, q = desc
        M = snappy.PlausibleKnots[name]
        K = snappy.RationalTangle(p, q).numerator_closure()
        match = M.is_isometric_to(K.exterior())
        success, m = is_square(p, True)
        obs = casson_gordon.could_be_ribbon(m, q)
        print(name, match, obs)
         

def concordant_to_minus_2_7a1(knot, band):
    """
    The knot K7a1 = 7_7 has infinite order in the topological
    concordance group by [Livingston-Naik, JDG 1999].  We check that
    band gives a ribbon concordance to -2 * K7a1.
    """
    T = snappy.Link('K7a1').mirror()
    E = T.exterior()

    K = knot.copy()
    normalize_crossing_labels(K)

    L = K.add_band(band)
    L.simplify('global')
    summands = L.deconnect_sum()
    if len(L.link_components) == 1 and len(summands) == 2:
        exteriors = [S.exterior() for S in summands]
        if all(orient_pres_isometric(S, E) for S in exteriors):
           return True

    return False


def test_eighteens():
    for name, band in alt_18s.items():
        K = snappy.PlausibleKnots[name].link()
        print(name, concordant_to_minus_2_7a1(K, band))


if __name__ == '__main__':
    test_two_bridge()
    test_eighteens()
