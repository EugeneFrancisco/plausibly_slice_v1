"""
We implement the Casson-Gordon ribbon obstruction from their 1986
paper "Cobordism of classical knots".
"""

from sage.all import (QQ, gcd, divisors, OrderedPartitions,
                      continued_fraction, is_square)
import snappy
import plausible_knots

def count_points(a, b):
    """
    Returns four times the number of lattice points in the triangle
    with vertices (0, 0), (0, a), (a, b/a), counted as in [CG].
    """
    assert b % a != 0
    count = 2*a - 1  # entire contribution from the bottom side
    for x in range(0, a + 1):
        y = 1 
        while a*a*y <= b*x:
            if a*a*y == b*x or x == a:
                count += 2
            else:
                count += 4
            y += 1
    return count

def could_be_ribbon(m, q):
    """
    Computes the [CG] ribbon obstruction for the 2-bridge knot q/m^2
    """
    assert m > 1 and m % 2 == 1 and q % 2 == 0 and 1 < q < m*m
    return all(2*q*r*r - count_points(m*r, q*r*r) in {-1, 1}
               for r in range(1, m))

def could_be_ribbon_full(p, q=None):
    if q is None:
        p, q = p
    sq, m = is_square(p, True)
    return False if not sq else could_be_ribbon(m, q)

def likely_ribbon_params(m):
    ans = set()
    for k in range(1, m):
        if gcd(k, m) == 1:
            ans.update({k*m + 1, k*m - 1})
    for e in {-1, 1}:
        for d in divisors(2*m + e):
            ans.add((m - e)*d)
        for d in divisors(m + e):
            if d % 2 == 1:
                r = (m + e)//d
                ans.update({(m + e)*d, (2*m - e)*r})
    return sorted(q for q in ans if q % 2 == 0 and q < m*m)


def test_conjecture(m):
    """
    This has been done for all m < 5000 by [Eisermann-Lamm 2009], but
    its a good way to check that all our formula are correct.
    """
    qs = [q for q in range(2, m*m, 2)
          if gcd(q, m) == 1 and could_be_ribbon(m, q)]
    return qs == likely_ribbon_params(m)


def plausibly_slice_two_bridge_knots(num_cross):
    for x in OrderedPartitions(num_cross):
        if x[-1] != 1:
            r = QQ(continued_fraction(x))
            p, q = r.numerator(), r.denominator()
            assert p > q
            if p % 2 == 1:
                K = snappy.RationalTangle(p, q).numerator_closure()
                M = K.exterior()
                if K.signature() == 0 and M.fox_milnor_test():
                    K.two_bridge = (p, q)
                    K.name = f'K({p}, {q})'
                    yield K


def get_plausible_knot_names(num_cross):
    seen = dict()
    for K in plausibly_slice_two_bridge_knots(num_cross):
        E = K.exterior()
        M = snappy.PlausibleKnots.identify(E)
        name = M.name()
        p, q = K.two_bridge
        if name in seen and seen[name][1] % 2 == 1 and q % 2 == 0:
            seen[name] = (p, q)
        if name not in seen:
            seen[name] = (p, q)

    if num_cross >= 17:
        names = sorted(seen.keys())
    else:
        names = sorted(seen.keys(), key=lambda x:int(x.split('a')[-1]))
    return {name:seen[name] for name in names}

        
    
