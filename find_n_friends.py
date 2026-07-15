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

def is_knot_exterior(manifold):
    """
    Here, we assume any S^3 slope has length < 4 rather than the
    proven bound of 6.  (The longest known such slope has length 3.3)
    """
    assert manifold.num_cusps() == 1
    if manifold.homology().elementary_divisors() != [0]:
        return False

    try:
        slopes = manifold.short_slopes(4.0)[0]
    except RuntimeError:
        slopes = []

    (a, b) = manifold.homological_longitude()
    slopes = [s for s in slopes if abs(a*s[1] - b*s[0]) == 1]
    for s in slopes:
        M = snappy.Triangulation(manifold)
        M.dehn_fill(s)
        if is_three_sphere(M):
            return s
    return False


def safe_is_isometric_to(A, B):
    """
    Work around when SnapPy is indecisive about when two manifolds are
    the same.
    """
    A, B = A.copy(), B.copy()
    if A.solution_type(enum=True) == 1 and B.solution_type(enum=True) == 1:
        vol_A, vol_B = A.volume(), B.volume()
        RR = vol_A.parent()
        assert vol_A.prec() == vol_B.prec()
        if abs(A.volume() - B.volume()) > RR(2)**(-RR(0.75)*RR.prec()):
            return False
        for i in range(10):
            try:
                return A.is_isometric_to(B)
            except RuntimeError:
                A.randomize(), B.randomize()
        return False


def safe_length_spectrum(manifold, max_len):
    """
    Do a little randomization for when SnapPy struggles to find the
    Dirichlet domain.
    """
    M = manifold.copy()
    ans = None
    for i in range(5):
        try:
            ans = M.length_spectrum(max_len, grouped=False, include_words=True)
            break
        except RuntimeError:
            M.randomize()
    return ans

Slope = tuple[int, int]


def _word_is_homology_generator(
        phi: nsagetools.MapToAbelianization, word: str) -> bool:
    """
    True if ``word`` generates H_1 of the closed n-surgery.

    ``phi`` is a MapToAbelianization of pi_1(Z_K^{(n)}), whose H_1 is Z/n
    (cyclic; Z when n = 0, trivial when n = 1).  A class generates a cyclic
    group Z/d iff its exponent is coprime to d, and gcd(e, 0) = |e| recovers
    the "generator of Z" test (|e| == 1) used for 0-friends.

    Args:
        phi: Abelianization map of pi_1(Z_K^{(n)}).
        word: Group element to test, as a word in the generators.

    Returns:
        True iff ``word``'s class generates H_1.
    """
    divisors = phi.elementary_divisors
    if not divisors:            # H_1 trivial (n = 1): every class generates.
        return True
    assert len(divisors) == 1, f"expected cyclic H_1, got {divisors}"
    exponent = int(phi(word).exponents()[0])
    return gcd(exponent, int(divisors[0])) == 1


def n_surgery_slope(exterior: snappy.Manifold, n: int) -> Slope:
    """
    The integer n-surgery slope of a knot exterior: n*meridian + longitude.

    Uses the homological (Seifert) longitude, so the slope is correct
    regardless of how the exterior's peripheral basis is oriented (the
    meridian is (1, 0)).  n = 0 gives the longitude, i.e. the 0-surgery.

    Args:
        exterior: A one-cusped knot exterior.
        n: The surgery coefficient.

    Returns:
        The filling slope (n + a, b) for homological longitude (a, b).
    """
    a, b = exterior.homological_longitude()
    return (n + a, b) # (0, 1) -> (n + a, b)


def _closed_isometric(A: snappy.ManifoldHP, B: snappy.ManifoldHP) -> bool:
    """
    True if the two closed hyperbolic manifolds are isometric.

    Unlike ``safe_is_isometric_to`` (which only compares cusped exteriors
    that solve to solution_type 1), a Dehn filling routinely lands at
    solution_type 2, and occasionally at a degenerate one where
    ``is_isometric_to`` throws.  is_isometric_to handles type 2 fine, so we
    reject only on a genuine volume mismatch and randomize past throws.

    Args:
        A, B: The closed manifolds to compare (left unmodified).

    Returns:
        True iff ``A`` and ``B`` are isometric.
    """
    A, B = A.copy(), B.copy()
    vol_A, vol_B = A.volume(), B.volume()
    RR = vol_A.parent()
    if abs(vol_A - vol_B) > RR(2)**(-RR(0.75)*RR.prec()):
        return False
    for _ in range(10):
        try:
            return bool(A.is_isometric_to(B))
        except RuntimeError:
            A.randomize(), B.randomize()
    return False


def _geometric_triangulation(
        M: snappy.ManifoldHP, tries: int = 30) -> snappy.ManifoldHP | None:
    """
    Randomize the closed manifold until it has a geometric (solution_type 1)
    triangulation, which ``drill_word`` requires; return None if unreachable.

    A Dehn filling often lands at solution_type 2 (negatively oriented
    tetrahedra).  Many such surgeries randomize to a type 1 triangulation of
    the same manifold; some (e.g. 10_125's 3-surgery) never do, so we give up
    after ``tries`` and the caller skips that surgery.

    Args:
        M: The closed manifold to retriangulate (modified in place).
        tries: Maximum number of randomization attempts.

    Returns:
        ``M`` once it is solution_type 1, or None if unreachable.
    """
    for _ in range(tries):
        if M.solution_type(enum=True) == 1:
            return M
        M.randomize()
    return None


def _n_surgery_recovers(
        exterior: snappy.Manifold, n: int, target: snappy.ManifoldHP) -> bool:
    """
    True if the integer (+/-)n-surgery on ``exterior`` is isometric to ``target``.

    Rejects a homology generator whose dual slope is a *rational* rather than
    an integer surgery.  Both integer signs are accepted: SnapPy stores one
    fixed chirality per knot, so a genuine n-friend may recover the target
    only through its -n surgery (its mirror's +n surgery), as happens for the
    10_132/10_125 and 6_2/K14n10164 examples.

    Args:
        exterior: The candidate friend's knot exterior.
        n: The surgery coefficient.
        target: The closed n-surgery Z_K^{(n)} to match.

    Returns:
        True iff the +n or -n surgery on ``exterior`` recovers ``target``.
    """
    for slope in (n_surgery_slope(exterior, n), n_surgery_slope(exterior, -n)):
        F = snappy.ManifoldHP(exterior)       # match target's precision
        F.dehn_fill(slope)
        if _closed_isometric(target, F):
            return True
    return False


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

    Args:
        manifold: The knot K, as a name or SnapPy manifold.
        n: The surgery coefficient.
        max_len: Only drill geodesics up to this complex length.
        min_len: Skip geodesics shorter than this length.

    Returns:
        One 4-tuple (drilled word, complex length, exterior volume, exterior
        isosig) per n-friend K' found, or None if Z_K^{(n)} has no drillable
        (geometric) triangulation.
    """
    if sys.getrecursionlimit() < 50000:
        sys.setrecursionlimit(50000)
    if isinstance(manifold, str):
        manifold = snappy.Manifold(manifold)

    E = snappy.ManifoldHP(manifold)
    M = E.copy()
    M.dehn_fill(n_surgery_slope(E, n))        # the n-surgery Z_K^{(n)}
    M = _geometric_triangulation(M)           # drill_word needs solution_type 1
    if M is None:                             # non-hyperbolic or no geometric
        return                                # triangulation reachable

    G = M.fundamental_group(False, False, False)
    phi = nsagetools.MapToAbelianization(G)
    geodesics = safe_length_spectrum(M, max_len)
    if geodesics is None:
        return

    other_knots = []
    print(f"Testing {len(geodesics)} geodesics for {n} friends.")
    count = 0
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
        count += 1
        if n != 0 and not _n_surgery_recovers(F, n, M):
            continue                          # rational-surgery coincidence
        F.randomize()
        other_knots.append((g.word,
                            complex(g.length),
                            F.volume(),
                            F.triangulation_isosig()))
    return other_knots
