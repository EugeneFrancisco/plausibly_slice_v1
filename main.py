import snappy as sn

# One of the 101 pairs of 0-friends from known_0_friends.py.  Their
# common 0-surgery is hyperbolic, so we can certify the match with
# SnapPy's isometry check.
KNOT_A = 'K12n309'
KNOT_B = 'K14n14254'

def zero_surgery(knot: sn.Manifold):
    """
    Return the 0-surgery of the named knot as a closed 3-manifold.
    """
    knot.dehn_fill((0, 1))
    return knot


def are_same_manifold(A, B):
    """
    True if A and B are (verified) isometric, hence homeomorphic by
    Mostow rigidity.  Retries with randomization since SnapPy's isometry
    check is occasionally indecisive.
    """
    for _ in range(5):
        try:
            return A.is_isometric_to(B)
        except RuntimeError:
            A.randomize()
            B.randomize()
    return False

def verify_0_friends(knot_a: str, knot_b: str) -> bool:
    """
    If returns true, then knot_a and knot_b are definitely friends.
    Otherwise, they might not be.

    ``knot_a`` and ``knot_b`` are knot names given as strings that
    SnapPy can look up via ``sn.Manifold(name)`` -- e.g. Rolfsen names
    like ``'K12n309'`` or census names like ``'K14n14254'``.
    """
    knot_a_manifold = sn.Manifold(knot_a)
    knot_b_manifold = sn.Manifold(knot_b)
    zero_a = zero_surgery(knot_a_manifold)
    zero_b = zero_surgery(knot_b_manifold)

    print(f'0-surgery of {knot_a}: {zero_a}')
    print(f'  volume = {zero_a.volume()},  homology = {zero_a.homology()}')
    print(f'0-surgery of {knot_b}: {zero_b}')
    print(f'  volume = {zero_b.volume()},  homology = {zero_b.homology()}')

    if are_same_manifold(zero_a, zero_b):
        return True
    else:
        return False

if __name__ == "__main__":
    main()
