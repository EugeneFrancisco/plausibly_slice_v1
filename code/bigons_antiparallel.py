"""
Searching for slice knots by adding F # F, where F = K4a1 = 4_1,
and then two bands.  See Section 2.6 and Figure 9 for more.

The main function is ``try_4_1s``.  See the very end of the file for
usage examples.
"""

import spherogram
import snappy
from spherogram.links.bands.search import ribbon_concordant_links, could_be_strongly_slice


def normalize_crossing_labels(link):
    for i, C in enumerate(link.crossings):
        C.label = i


def pos_01_bigon(C):
    """
    if the 0 and 1 ends of a positive crossing C are the 1 and 0 ends of another positive crossing
    this returns that crossing. Otherwise, it returns None
    """
    if C.sign == 1:
        A, a = C.adjacent[1]
        B, b = C.adjacent[0]
        if A == B and A.sign == 1 and a == 0:
            return A
    return None


def pos_23_bigon(C):
    """
    if the 2 and 3 ends of a positive crossing C are the 3 and 2 ends of another positive crossing
    this returns that crossing. Otherwise, it returns None
    """
    if C.sign == 1:
        A, a = C.adjacent[2]
        B, b = C.adjacent[3]
        if A == B and A.sign == 1 and a == 3:
            return A
    return None


def neg_12_bigon(C):
    """
    if the 1 and 2 ends of a negative crossing C are the 2 and 1 ends of another negative crossing
    this returns that crossing. Otherwise, it returns None
    """
    if C.sign == -1:
        A, a = C.adjacent[1]
        B, b = C.adjacent[2]
        if A == B and A.sign == -1 and a == 2:
            return A
    return None


def neg_03_bigon(C):
    """
    if the 0 and 3 ends of a negative crossing C are the 3 and 0 ends of another negative crossing
    this returns that crossing. Otherwise, it returns None
    """
    if C.sign == -1:
        A, a = C.adjacent[0]
        B, b = C.adjacent[3]
        if A == B and A.sign == -1 and a == 3:
            return A
    return None
    

def antiparallel_bigon_ends(diagram):
    """
    returns a tuple. The first is a list of the positive chains, each
    of which is a list of crossings. The second is the list of negative
    chains, each of which is a list of crossings.
    >>> L = snappy.Link('K12a990')
    >>> antiparallel_bigon_ends(L)
    ([[10, 11]], [[2, 3]])
    """
    pos_bigon_ends_01, pos_bigon_ends_23, neg_bigon_ends_12, neg_bigon_ends_03= [], [], [], []
    for C in diagram.crossings:
        if C.sign == 1:
            if pos_01_bigon(C) is not None and pos_23_bigon(C) is None:
                pos_bigon_ends_01.append(C)
            if pos_23_bigon(C) is not None and pos_01_bigon(C) is None:
                pos_bigon_ends_23.append(C)

        if C.sign == -1:
            if neg_03_bigon(C) is not None and neg_12_bigon(C) is None:
                neg_bigon_ends_03.append(C)
            if neg_12_bigon(C) is not None and neg_03_bigon(C) is None:
                neg_bigon_ends_12.append(C)
    
    return pos_bigon_ends_01 + pos_bigon_ends_23, neg_bigon_ends_12 + neg_bigon_ends_03


def insert_4_1_to_neg(link, crossing, s_over, s_under):
    """
    insert a 4_1 to a crossing that is part of a negative antiparallel bigon
    where the insertion is at strands s_over and s_under (as in the diagram)
    the insertion is so that s_over and s_under end up connected to a positive crossing
    in the added 4_1

    s_over and s_under should be one incoming and one outgoing strand
    """
    L = link.copy()
    crossing_index = link.crossings.index(crossing)
    C = L.crossings[crossing_index]

    orig_num_crossings = len(link.crossings)
    i = orig_num_crossings
    A, a = C.adjacent[s_over]
    B, b = C.adjacent[s_under]

    assert A != C and B != C


    c1 = spherogram.Crossing()
    c1.label = i
    c2 = spherogram.Crossing()
    c2.label =i+1
    c3 = spherogram.Crossing()
    c3.label = i+2
    c4 = spherogram.Crossing()
    c4.label = i+3
    
    for c in [c1, c2, c3, c4]:
        L.crossings.append(c)

    c1[1], c1[0] = C[s_over], C[s_under]
    c1[2], c1[3] = c2[3], c2[2]
    c2[0], c2[1] = c4[3], c3[0]
    c3[1], c3[2] = c4[2], c4[1]
    c4[0] = A[a]
    c3[3] = B[b]

    L._rebuild()
    return(L)


def insert_4_1_to_neg_both(link, crossing):
    return (insert_4_1_to_neg(link, crossing, 1, 2), insert_4_1_to_neg(link, crossing, 3, 0))


def insert_4_1_to_pos_both(link, crossing):
    """
    returns two links each of which consist of the original link with 4_1
    inserted near near crossing (which should be a positive crossing that is
    an end of an anti-parallel bigon.
    the insertion is at a negative antiparallel bigon of 4_1
    and is at two strands of crossing one of which is incoming and one outgoing
    """
    
    L = link.copy()
    L = L.mirror()
    crossing_index = link.crossings.index(crossing)
    C = L.crossings[crossing_index]
    mirror_inserted1 = insert_4_1_to_neg(L, C, 1, 2)
    mirror_inserted2 = insert_4_1_to_neg(L, C, 3, 0)

    return (mirror_inserted1.mirror(), mirror_inserted2.mirror())



def insert_6_1_to_neg1(link, crossing, s_over, s_under):
    """
    Inserts the 6_1 along the two half-twist part. Inserts one where the
    two half-twists have positive sign.
    """
    
    L = link.copy()
    crossing_index = link.crossings.index(crossing)
    C = L.crossings[crossing_index]


    orig_num_crossings = len(link.crossings)
    i = orig_num_crossings
    A, a = C.adjacent[s_over]
    B, b = C.adjacent[s_under]

    assert A != C and B != C


    c1 = spherogram.Crossing()
    c1.label = i
    c2 = spherogram.Crossing()
    c2.label =i+1
    c3 = spherogram.Crossing()
    c3.label = i+2
    c4 = spherogram.Crossing()
    c4.label = i+3
    c5 = spherogram.Crossing()
    c5.label = i+4
    c6 = spherogram.Crossing()
    c6.label = i+5

    
    for c in [c1, c2, c3, c4, c5, c6]:
        L.crossings.append(c)

    c1[1], c1[0] = C[s_over], C[s_under]
    c1[2], c1[3] = c2[3], c2[2]
    c2[0], c2[1] = c3[3], c6[0]
    
    c3[1], c3[2] = c4[2], c4[1]
    c4[0], c4[3] = c5[3], c5[0]
    c5[1], c5[2] = c6[2], c6[1]
    
    c3[0] = A[a]
    c6[3] = B[b]


    L._rebuild()
    return(L)


def insert_6_1_to_pos1(link, crossing):
    L = link.copy()
    L = L.mirror()
    crossing_index = link.crossings.index(crossing)
    C = L.crossings[crossing_index]

    mirror_inserted1 = insert_6_1_to_neg1(L, C, 1, 2)
    mirror_inserted2 = insert_6_1_to_neg1(L, C, 3, 0)

    return (mirror_inserted1.mirror(), mirror_inserted2.mirror())


def insert_6_1_to_neg1_both(link, crossing):
    return (insert_6_1_to_neg1(link, crossing, 1, 2), insert_6_1_to_neg1(link, crossing, 3, 0))


def try_4_1s(knot):
    assert len(knot.link_components) == 1 and knot.unlinked_unknot_components == 0
    knot = knot.copy()
    normalize_crossing_labels(knot)
    pos_bigon_ends0, neg_bigon_ends0 = antiparallel_bigon_ends(knot)

    pos_bigon_ends = [(end, 1) for end in pos_bigon_ends0]
    neg_bigon_ends = [(end, -1) for end in neg_bigon_ends0]
    bigon_ends = pos_bigon_ends + neg_bigon_ends

    
    certificates = {knot:[knot.PD_code()]}

    
    for (C1, sign1) in bigon_ends:
        for (C2, sign2) in bigon_ends:
            C2_index = knot.crossings.index(C2)

            if sign1 == -1:
                L1, L2 = insert_4_1_to_neg_both(knot, C1)
            else:
                L1, L2 = insert_4_1_to_pos_both(knot, C1)
            assert L1.crossings[C2_index].label == C2_index
            assert L2.crossings[C2_index].label == C2_index


            if sign2 == 1:
                L3, L4 = insert_4_1_to_pos_both(L1, L1.crossings[C2_index])
                L5, L6 = insert_4_1_to_pos_both(L2, L2.crossings[C2_index])

            else:
                L3, L4 = insert_4_1_to_neg_both(L1, L1.crossings[C2_index])
                L5, L6 = insert_4_1_to_neg_both(L2, L2.crossings[C2_index])


            for L in [L3, L4, L5, L6]:
                if len(L.link_components) != 3:
                    continue

                L_PD_code = L.PD_code()

                #link_with_4_1s.view(show_crossing_labels=True)
                L.simplify('global')
                if len(L.crossings) > len(knot.crossings):
                    continue

                L.unlinked_unknot_components = 0
                if len(L.link_components) == 0:
                    certificates['unknot'] = [knot.PD_code(), 'inserted 4_1s at '+str(C1) + ' and ' + str(C2), L_PD_code, 'unknot']
                    return certificates

                if not could_be_strongly_slice(L):
                    continue
                print("Adding bands to " + str(L))
                ans = ribbon_concordant_links(L,
                                              max_bands=1,
                                              max_twists=1,
                                              max_band_len=5,
                                              paths='shortest',
                                              stop_at_unlink=True,
                                              R1_only=False,
                                              certify=True,
                                              print_progress=False)
                if 'unknot' in ans:
                    certificates['unknot'] = [knot.PD_code(), 'inserted 4_1s at '+str(C1) + ' and ' + str(C2), L_PD_code, 'unknot']
                    return certificates

    return certificates


if __name__ == '__main__':
    # This is 19nh_205383898
    K = snappy.Link('DT: sasHJILMRDAGBpCFSqnkoE.1000101011110101011')
    print(antiparallel_bigon_ends(K))
    
    K1 = snappy.Link("6_1")
    K1 = snappy.Link("4_1").connected_sum(snappy.Link("4_1"))
    print(try_4_1s(K1))

    for knot_name in ["K6a3", "K8a4", "K8a16", "K8n1"]:
        print("Trying: " + knot_name)
        K1 = snappy.Link(knot_name)
        print(try_4_1s(K1))
