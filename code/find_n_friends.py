"""
Finding n-friends of a hyperbolic knot: knots whose n-surgeries agree.

This generalizes ``find_0_friends.py`` from the 0-surgery to the integer
n-surgery for any n >= 0.  For a knot K, let Z_K^{(n)} be the n-surgery.
Two knots K, K' are n-friends when Z_K^{(n)} = Z_{K'}^{(n)}.

The search mirrors the 0-surgery case exactly, with two changes:

    * We fill the n-surgery slope n*meridian + longitude instead of the
      longitude (0-surgery); see ``n_surgery_slope``.
    * H_1(Z_K^{(n)}) = Z/n is finite for n > 0 (it is Z for n = 0), so the
      geodesic we drill must be a *generator* of that homology; we test
      this with the full abelianization instead of the free one.  The dual
      of an n-friend is always such a generator.

Because a homology generator can still fill back along a rational or
opposite-sign slope, we additionally verify that the integer n-surgery of
the recovered knot really recovers Z_K^{(n)} before recording it (skipped
when n = 0, where the longitude is the only candidate slope).

Example usage::

  sage: run find_n_friends.py
  sage: ans = find_common_n_surgery_via_words('K11n34', 0, 3)   # 0-friends
  sage: ans = find_common_n_surgery_via_words('K11n34', 2, 3)   # 2-friends
"""

import sys
from math import gcd

import snappy
from snappy.snap import nsagetools
from tqdm import tqdm

from find_0_friends import (
    is_knot_exterior,
    safe_is_isometric_to,
    safe_length_spectrum,
)


def _word_is_homology_generator(phi, word) -> bool:
    """
    True if ``word`` generates H_1 of the closed n-surgery.

    ``phi`` is a MapToAbelianization of pi_1(Z_K^{(n)}), whose H_1 is Z/n
    (cyclic; Z when n = 0, trivial when n = 1).  A class generates a cyclic
    group Z/d iff its exponent is coprime to d, and gcd(e, 0) = |e| recovers
    the "generator of Z" test (|e| == 1) used for 0-friends.
    """
    divisors = phi.elementary_divisors
    if not divisors:            # H_1 trivial (n = 1): every class generates.
        return True
    assert len(divisors) == 1, f"expected cyclic H_1, got {divisors}"
    exponent = int(phi(word).exponents()[0])
    return gcd(exponent, int(divisors[0])) == 1


def n_surgery_slope(exterior, n: int) -> tuple[int, int]:
    """
    The integer n-surgery slope of a knot exterior: n*meridian + longitude.

    Uses the homological (Seifert) longitude, so the slope is correct
    regardless of how the exterior's peripheral basis is oriented (the
    meridian is (1, 0)).  n = 0 gives the longitude, i.e. the 0-surgery.
    """
    a, b = exterior.homological_longitude()
    return (n + a, b)


def _n_surgery_recovers(exterior, n: int, target) -> bool:
    """
    True if the integer n-surgery on ``exterior`` is isometric to ``target``.

    Guards against a homology generator whose dual slope is a rational or
    (-n)-surgery rather than the integer n-surgery we are after.
    """
    F = snappy.ManifoldHP(exterior)           # match target's precision
    F.dehn_fill(n_surgery_slope(exterior, n))
    return bool(safe_is_isometric_to(target, F))


def find_common_n_surgery_via_words(
        manifold: str | snappy.Manifold,
        n: int,
        max_len: float = 3.0,
        min_len: float = 0.0,
    ) -> list[tuple[str, complex, float, str]] | None:
    """
    Find n-friends of a knot whose n-surgery is hyperbolic.

    Generalizes ``find_common_zero_surgery_via_words`` to any n >= 0; the
    two agree at n = 0.  For a knot K, let Z_K^{(n)} be the n-surgery.
    Drills the short geodesics of Z_K^{(n)} and keeps those whose exterior
    is again a knot exterior K' with n-surgery equal to Z_K^{(n)}.

    Returns a list of 4-tuples encoding each n-friend K':

    1. Word in pi_1(Z_K^{(n)}) for the geodesic that was drilled.
    2. The complex length of said geodesic.
    3. The volume of the exterior of K'.
    4. A triangulation isosig of the exterior of K'.
    """
    if sys.getrecursionlimit() < 50000:
        sys.setrecursionlimit(50000)
    if isinstance(manifold, str):
        manifold = snappy.Manifold(manifold)

    E = snappy.ManifoldHP(manifold)
    M = E.copy()
    M.dehn_fill(n_surgery_slope(E, n))        # the n-surgery Z_K^{(n)}
    if M.solution_type(enum=True) != 1:       # need a hyperbolic n-surgery
        return

    G = M.fundamental_group(False, False, False)
    phi = nsagetools.MapToAbelianization(G)
    geodesics = safe_length_spectrum(M, max_len)
    if geodesics is None:
        return

    other_knots = []
    print(f"Testing {len(geodesics)} geodesics for {n} friends.")
    for g in tqdm(geodesics, desc=f"Testing geodesics for {n} friends"):
        if g.length.real() < min_len:
            continue
        if not _word_is_homology_generator(phi, g.word):
            continue
        F = M.drill_word(g.word).filled_triangulation()
        if F.solution_type(enum=True) not in [1, 2]:
            continue
        if safe_is_isometric_to(E, F):        # the geodesic recovers K itself
            continue
        slope = is_knot_exterior(F)
        if not slope:
            continue
        # Reframe F as the exterior of K' (meridian = its S^3 slope).
        F.dehn_fill(slope)
        F.set_peripheral_curves('fillings')
        F.dehn_fill((0, 0))
        if n != 0 and not _n_surgery_recovers(F, n, M):
            continue                          # rational-surgery coincidence
        F.randomize()
        other_knots.append((g.word,
                            complex(g.length),
                            F.volume(),
                            F.triangulation_isosig()))
    return other_knots
