"""
Searching for slice knots by adding T # -T, where T is trefoil, and
then two bands.  See Section 2.6 and Figure 8 for more.

The main function is ``try_paired_twists``.
"""

import spherogram
import snappy
import plausible_knots
from spherogram.links.bands.search import ribbon_concordant_links
from spherogram.links.bands.merge_links import link_isotopy_classes

def normalize_crossing_labels(link):
    for i, C in enumerate(link.crossings):
        C.label = i

def extend_chain_forward(C):
    if C.sign == 1:
        A, a = C.adjacent[1]
        B, b = C.adjacent[2]
        if A == B and A.sign == 1:
            return A
    if C.sign == -1:
        A, a = C.adjacent[2]
        B, b = C.adjacent[3]
        if A == B and A.sign == -1:
            return A
    return False

def can_extend_backwards(C):
    if C.sign == 1:
        A, a = C.adjacent[3]
        B, b = C.adjacent[0]
        if A == B and A.sign == 1:
            return True
    if C.sign == -1:
        A, a = C.adjacent[0]
        B, b = C.adjacent[1]
        if A == B and A.sign == -1:
            return True
    return False

def forward_chain(crossing):
    chain = [crossing]
    while A := extend_chain_forward(chain[-1]):
        chain.append(A)
    return chain

def oriented_bigon_chains(diagram):
    """
    returns a tuple. The first is a list of the positive chains, each
    of which is a list of crossings. The second is the list of negative
    chains, each of which is a list of crossings.
    >>> L = snappy.Link('K12a990')
    >>> oriented_bigon_chains(L)
    ([[10, 11]], [[2, 3]])
    """
    pos_bigons, neg_bigons = [], []
    for C in diagram.crossings:
        if not can_extend_backwards(C):
            if chain := forward_chain(C):
                if len(chain) > 1:
                    if C.sign == 1:
                        pos_bigons.append(chain)
                    else:
                        neg_bigons.append(chain)

    return pos_bigons, neg_bigons


def insert_pos_twist(link, crossing):
    L, C = link, crossing
    N = spherogram.Crossing()
    link.crossings.append(N)
    if C.sign == 1:
        A, a = C.adjacent[3]
        B, b = C.adjacent[0]
        assert A != C and B != C
        N[2], N[1] = C[3], C[0]
        N[3], N[0] = A[a], B[b]
    else:
        A, a = C.adjacent[0]
        B, b = C.adjacent[1]
        assert A != C and B != C
        N[2], N[1] = C[0], C[1]
        N[3], N[0] = A[a], B[b]


def insert_neg_twist(link, crossing):
    L, C = link, crossing
    N = spherogram.Crossing()
    link.crossings.append(N)
    if C.sign == 1:
        A, a = C.adjacent[3]
        B, b = C.adjacent[0]
        assert A != C and B != C
        N[3], N[2] = C[3], C[0]
        N[0], N[1] = A[a], B[b]
    else:
        A, a = C.adjacent[0]
        B, b = C.adjacent[1]
        assert A != C and B != C
        N[3], N[2] = C[0], C[1]
        N[0], N[1] = A[a], B[b]


def insert_paired_twists(link, pos_chain, neg_chain, num_twists):
    """
    The example from the HKL paper, demonstrating that K12a990 is
    slice without immediately producing a ribbon disk.

    >>> K = snappy.Link('K12a990')
    >>> pos_chains, neg_chains = oriented_bigon_chains(K)
    >>> L = insert_paired_twists(K, pos_chains[0], neg_chains[0], 3)
    >>> ans = L.simplify('global'); L
    <Link K12a990: 0 comp; 0 cross>
    >>> L.unlinked_unknot_components
    3

    WARNING: We do not bother to ensure that the resulting components
    are all oriented compatibly with the original link.

    """
    L = link.copy()
    p = link.crossings.index(pos_chain[0])
    n = link.crossings.index(neg_chain[0])
    P, N = L.crossings[p], L.crossings[n]
    assert P.sign == 1 and N.sign == -1

    for i in range(num_twists):
        insert_neg_twist(L, P)
        insert_pos_twist(L, N)

    L._rebuild(same_components_and_orientations=False)
    return L


def try_paired_twists(knot):
    """
    >>> K = snappy.Link('K12a990')
    >>> 'unknot' in try_paired_twists(K)
    True
    """
    assert len(knot.link_components) == 1 and knot.unlinked_unknot_components == 0
    knot = knot.copy()
    normalize_crossing_labels(knot)
    pos_chains, neg_chains = oriented_bigon_chains(knot)
    ans = dict()
    new_links = []
    for pos_chain in pos_chains:
        for neg_chain in neg_chains:
            max_twists = max(len(pos_chains), len(neg_chains)) + 3
            for num_twists in range(3, max_twists + 1, 2):
                link_with_twists = insert_paired_twists(knot, pos_chain, neg_chain, num_twists)
                if len(link_with_twists.link_components) == 3:
                    link_with_twists.simplify('global')
                    # cap off any trivial components
                    link_with_twists.unlinked_unknot_components = 0
                    if (len(link_with_twists.link_components) == 0 or
                        bands.is_unlink_exterior(link_with_twists.exterior())):
                        ans['unknot'] = [knot.PD_code(),
                                          num_twists,
                                          [C.label for C in pos_chain],
                                          [C.label for C in neg_chain]]
                        return ans
                    elif len(link_with_twists.link_components) == 1:
                        id1 = snappy.PlausibleKnots.identify(link_with_twists.exterior())
                        id2 = link_with_twists.exterior().identify()
                        if id2 and not id1:
                            print('MOVED OUT OF SAMPLE?!??!')
                        ans[id1.name()] = [knot.PD_code(),
                                              num_twists,
                                              [C.label for C in pos_chain],
                                              [C.label for C in neg_chain]]
                    else:
                        link_with_twists.twist_spec = [knot.PD_code(),
                                                       num_twists,
                                                       [C.label for C in pos_chain],
                                                       [C.label for C in neg_chain]]
                        new_links.append(link_with_twists)

    if 'unknot' in ans:
        return ans

    new_links = [L for L in new_links if len(L.crossings) < len(knot.crossings)]
    links_w_mflds = [(L, L.exterior()) for L in new_links if bands.could_be_strongly_slice(L)]
    new_links_w_mflds = link_isotopy_classes(links_w_mflds)
    new_links = sorted([L for L, _ in new_links_w_mflds], key=lambda x:len(x.crossings))
    for D in new_links:
        ans = ribbon_concordant_links(D,
                                      max_bands=2,
                                      max_twists=2,
                                      max_band_len=7,
                                      paths='shortest',
                                      stop_at_unlink=True,
                                      certify=True,
                                      print_progress=False)
        if 'unknot' in ans:
            ans = {'unknot':D.twist_spec + ans['unknot']}
    return ans


if __name__ == '__main__':
    import doctest
    print(doctest.testmod())
