"""Check whether the 1-surgery of 6_2 equals the 1-surgery of K13n3596."""

import os
import sys

import snappy

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'code'))
from find_0_friends import safe_is_isometric_to
from find_n_friends import n_surgery_slope


A = snappy.Manifold("6_2")
A.dehn_fill((1, 1))


B = snappy.Manifold('K13n3596')
B.dehn_fill((1, 1))

print(A.volume())
print(B.volume())
print(A.is_isometric_to(B))
