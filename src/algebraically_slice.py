"""
The question of whether a knot is algebraically slice is
algorithmically decidable.  See e.g.

[L] Chuck Livingston, "The algebraic concordance order of a knot"
    https://arxiv.org/abs/0806.3068v2

Most of this goes back to Levine's original work on algebraic
concordance, which reduces the question primarily to computing whether
certain quadratic forms over p-adic fields are 0 in the Witt group.
The whole thing is a bit involved, and no one seems to have completely
implemented it.

A probabilistic approach to finding certificates that a knot is
algebraically slice is described in Remark 16 of

https://arxiv.org/abs/1509.07634

and implemented by Lewark in

https://github.com/LLewark/galg-computations

However, we had no success applying Lewark's code to 19-crossing
plausibly slice knots.

The code here uses Lewark's idea as a starting point but goes in a
different direction.  We also implemented the first step of Levine's
approach, which is to replace a singular Seifert matrix by smaller
non-singular one that is Witt equivalent to it.
"""

import snappy, plausibly_slice
import random
from sage.quadratic_forms.qfsolve import qfsolve
from sage.all import (identity_matrix, ZZ, QQ, Integer, matrix, block_matrix,
                      block_diagonal_matrix, gcd, Matrix, vector, PariError,
                      extend_to_primitive, Permutation, elementary_matrix)


def test_idea(n):
    """
    Double-checking that the "block" version of Levine's lemma can
    be thought of a recursively applying the case where there is
    just one row of zeros.
    """
    V = matrix(ZZ, 3*n, 3*n)
    for i in range(n):
        V[i + n, i] = 1
    for i in range(n, 3*n):
        for j in range(n, 3*n):
            V[i, j] = 9

    sigma = Permutation(list(range(1, 2*n, 2)) + list(range(2, 2*n + 1, 2)) + list(range(2*n + 1, 3*n + 1)))
    P = sigma.to_matrix()
    W = P * V * P.transpose()

    for i in range(n):
        target_col = identity_matrix(W.nrows()).column(1)
        if not (W.row(0) == 0 and W.column(0) == target_col):
            return False
        W = W[2:, 2:]

    return True


def invertible_representative(seifert_matrix):
    """
    The first step in the approach in [L] is to replace the initial
    Seifert matrix with a smaller one that is invertable using the
    approach on page 3.

    To keep the sizes of the entries of the final matrix under
    control, we use LLL at various points.  We also use a "block
    matrix" form of Levine's lemma to reduce the amount of recursion.
    """
    V = seifert_matrix
    n = V.nrows()
    assert V.base_ring() == ZZ
    assert (V - V.transpose()).det().abs() == 1
    assert (V + V.transpose()).rank() == n

    # Start with the sublattice L that is left-orthogonal to
    # everything, i.e. the v with v*V = 0.

    L = V.left_kernel_matrix()
    k = L.nrows()
    if k == 0:
        return V
    L = L.LLL()  # The basis is the rows.

    # Now look at the left-orthogonal complement of L, that is,
    #
    #    U = {u with u * V * v = 0 for all v in L}

    U = matrix(V * L.transpose()).left_kernel_matrix()
    assert U * V * L.transpose() == 0

    # L is a submodule of U, and we now find a basis of U that extends
    # the one for L.

    L_in_U = U.solve_left(L).change_ring(ZZ)
    E = extend_to_primitive(L_in_U).change_ring(ZZ)
    assert E.det().abs() == 1

    # Another basis for U, this time with the basis for L for the
    # first vectors

    U = E * U
    assert U[:k, :] == L
    # cleanup
    U = L.stack(U[k: , :].LLL())

    # Extend U to a Z-basis of Z^n where the new vectors give the
    # standard dual basis of Hom(L, Z) via a -> (x -> a*V*x).

    R = extend_to_primitive(U).change_ring(ZZ)
    A = R[U.nrows(): , :].LLL()
    F = A * V * L.transpose()
    assert F.det().abs() == 1
    A = F.inverse().change_ring(ZZ) * A

    C = L.stack(A).stack(U[k: , :])
    assert C.det().abs() == 1

    # We apply the block matrix generalization of the lemma in Levine's
    # Inventiones 1969 paper on page 101, to see that B is equivalent
    # in G^Q and hence G^Z. The block version follows from the
    # original, see the "test_idea" function for more.

    W = C * V * C.transpose()
    assert W[:k, :] == 0 and W[k:2*k, 0:k] == 1 and W[2*k: , :k] == 0
    B = W[2*k:, 2*k:]
    assert (B - B.transpose()).det().abs() == 1
    return invertible_representative(B)


def random_elementary_matrix(size):
    indices = list(range(size))
    i = random.choice(indices)
    indices.pop(i)
    j = random.choice(indices)
    s = random.choice([-1, 1])
    return elementary_matrix(ZZ, size, row1=i, row2=j, scale=s)


def random_walk_SLnZ(size, max_rand_steps):
    A = identity_matrix(size)
    L = random.randrange(max_rand_steps)
    for i in range(L):
        A = A * random_elementary_matrix(size)
    return A


def randomized_qfsolve(Q, rand_steps, max_failures=10):
    # Similar to https://github.com/LLewark/galg-computations, but we
    # pick the random C in a different way that makes the entries much
    # smaller.
    for i in range(max_failures):
        C = random_walk_SLnZ(Q.nrows(), rand_steps)
        W = C.transpose() * Q * C
        try:
            sol = W.__pari__().qfsolve()
        except PariError:
            continue
        if sol.type() == 't_INT':
            continue
        if sol.type() == 't_COL':
            v = vector(QQ, sol)
        elif sol.type() == 't_MAT':
            v = matrix(QQ, sol).column(0)
        return v / gcd(v)


def algebraically_slice(knot, tries=10000):
    """
    Quick and dirty test that can certify a knot as algebraically slice.
    When it returns True, the knot is algebraically slice but False is
    inconclusive.

    If Q = V + V^t and T = V^-1 * V^t, then being algebraically slice
    is equivalent to there being a metabolizer for Q which is
    T-invariant.

    Our strategy is to pick a random vector v with v^t * Q * v = 0
    and look at its T-orbit.
    """
    V = invertible_representative(knot.seifert_matrix())
    if V.nrows() == 0:
        return True
    Q = V + V.transpose()
    T = V.inverse() * V.transpose()
    k = Q.nrows()//2
    failures = 0
    for i in range(tries):
        v = randomized_qfsolve(Q, 15)
        if v is None:
            break
        W = matrix([T**n * v for n in range(k)])
        if W * Q * W.T == 0:
            w = T * W.row(-1)
            R = W.row_space()
            if w in R and R.dimension() == k:
                return True

    return False


def test_for_assertion_failures(n):
    for i in range(n):
        M = snappy.PlausibleKnots.random()
        K = M.link()
        V = K.seifert_matrix()
        W = invertible_representative(V)
        print(M.name(), V.nrows(), W.nrows(), max(abs(a) for a in W.list()))


def test_false_positives():
    for M in snappy.HTLinkExteriors(cusps=1):
        if M.solution_type().startswith('all'):
            if algebraically_slice(M.link()):
                ident = snappy.PlausibleKnots.identify(M)
                print(M, ident)
                assert ident


def test_for_new_failures(n):
    for i in range(n):
        M = snappy.PlausibleKnots.random()
        K = M.link()
        print(M.name(), algebraically_slice(K))
