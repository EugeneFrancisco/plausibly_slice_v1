"""
Worked examples for the n-friend search (``Knot.find_n_friends``).

Two knots are n-friends when their integer n-surgeries are the same closed
manifold -- the surgery generalization of the 0-friends the paper uses.
See ``code/find_n_friends.py`` and ``Knot.find_n_friends`` /
``Knot.is_n_friend`` / ``Knot.n_surgery``.

Run in the ``sage`` conda env:

    conda activate sage
    python test.py
"""

from __future__ import annotations

from knot import Knot


def example_zero_friend():
    """
    n = 0 is exactly the paper's 0-friend search.

    find_n_friends(0) runs the same code path as the original 0-surgery
    search, recovering the known 0-friend of K11n34 and certifying it.
    """
    print("=" * 68)
    print("Example 1: n = 0 recovers the paper's 0-friend of K11n34")
    print("=" * 68)
    K = Knot("K11n34")
    K.find_n_friends(0, 3)
    friends = K.n_friends.get(0, [])
    print(f"    0-friends found: {len(friends)}")
    for f in friends:
        print(f"      certified 0-friend? {K.is_n_friend(f, 0)}   "
              f"vol(exterior) = {f.exterior().volume():.5f}")
    return K


def example_is_n_friend(K):
    """Certify the discovered pair at several n directly with is_n_friend."""
    print("\n" + "=" * 68)
    print("Example 2: is_n_friend distinguishes a 0-friend from an n-friend")
    print("=" * 68)
    if not K.n_friends.get(0):
        print("    (no 0-friend was found this run; skipping)")
        return
    f = K.n_friends[0][0]
    for n in (0, 1, 2, 3):
        # The pair shares the 0-surgery, but generally not the n-surgery.
        print(f"    K11n34 and its 0-friend share their n={n}-surgery? "
              f"{K.is_n_friend(f, n)}")


def example_positive_n_search(names=("K11n34", "K11n42"), ns=(2, 3)):
    """
    Run the generalized search at n > 0.

    Integer n-surgeries are far more rigid than the 0-surgery, so short
    n-friends are rare; the search reports whatever it certifies.  Every
    reported friend passes the homology-generator filter and the explicit
    "integer n-surgery recovers the original" check, so is_n_friend is True.
    """
    print("\n" + "=" * 68)
    print("Example 3: the generalized search at n > 0")
    print("=" * 68)
    for name in names:
        for n in ns:
            K = Knot(name)
            K.find_n_friends(n, 3)
            friends = K.n_friends.get(n, [])
            print(f"    {name}, n = {n}: {len(friends)} n-friend(s) found", end="")
            if not friends:
                print()
            for f in friends:
                print(f"  ->  certified? {K.is_n_friend(f, n)}   "
                      f"vol = {f.exterior().volume():.5f}")
    print("    (no short n-friends here: unlike the 0-surgery, integer")
    print("     n-surgeries are rigid, so genuine n-friends are rare.)")


if __name__ == "__main__":
    example_is_n_friend(example_zero_friend())
    example_positive_n_search()
