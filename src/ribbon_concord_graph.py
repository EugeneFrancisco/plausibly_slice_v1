"""

This file uses the ribbon concordance graph Gamma_{PS_19} of Section 7
to propogate information and complete the proof of Theorem 1.9.

The function main at the bottom runs all tests and is automatically
invoked if you run this Python file as script::

  python ribbon_concord_graph.py

"""

import pandas as pd
import networkx as nx
import snappy
from spherogram.links.bands.core import uncompress_band_spec


def examine_ribbon_cert(cert):
    if isinstance(cert, str):
        cert = eval(cert)
    assert len(cert) % 2 == 1

    if cert[-1] != 'unknot':
        R = snappy.RibbonLinks[cert[-1]]
        cert = cert[:-1] + R.ribbon_cert

    num_bands = len(cert) // 2
    crossings = [len(cert[i]) for i in range(0, len(cert), 2)]
    its_bands = [uncompress_band_spec(cert[i])
                 for i in range(1, len(cert), 2)]

    return (num_bands,
            crossings,
            max([len(b[0]) for b in its_bands]),
            max([abs(b[-1]) for b in its_bands]))


def examine_ribbon(df):
    """
    Builds the ribbon concordance directed graph Gamma_{PS_19} of
    Section 7 of the paper from the data in the ``concord_to`` and
    ``ribbon_cert`` columns.

    Returns Gamma_{PS_19} and its connected components.
    """
    df.ribbon_to = df.ribbon_to.fillna('[]')
    knot_sample = set(df.name)
    knot_sample.add('unknot')

    concordances = nx.DiGraph()

    # Add the edges with the larger number of bands first so that each
    # edge is labled by the smallest number of bands we've seen.


    for i, row in df[df.ribbon_cert.notnull()].iterrows():
        knot = row['name']
        num_bands = examine_ribbon_cert(row['ribbon_cert'])[0]
        concordances.add_edge(knot, 'unknot', bands=num_bands)

    def add_edges(row):
        for im in eval(row['ribbon_to']):
            concordances.add_edge(row['name'], im, bands=1)

    df.apply(add_edges, axis=1)
    assert set(concordances.nodes()).issubset(knot_sample)

    print(f'Nodes: {concordances.number_of_nodes()} '
          f'Edges: {concordances.number_of_edges()}')
    # One of Gordon's conjectures, since proved by Ian Agol.
    assert nx.is_directed_acyclic_graph(concordances)

    # Check that connected components have a single sink as a directed
    # graph.

    components = list(nx.connected_components(concordances.to_undirected()))
    for comp in components:
        C = concordances.subgraph(comp)
        u = list(nx.topological_sort(C))[-1]
        ancestors = nx.ancestors(C, u)
        if len(ancestors) + 1 != C.number_of_nodes():
            print(f'Component without unique sink of size {C.number_of_nodes()}')
            print(f'Expected minimal element was {u}, but these are not above it:')
            print(set(C.nodes()) - set(ancestors) - {u})

    # Update the dataframe
    ribbon = nx.ancestors(concordances, 'unknot')
    def num_bands(K):
        path = nx.shortest_path(concordances, K, 'unknot',
                                weight='bands', method='bellman-ford')
        num_bands = 0
        for i in range(len(path) - 1):
            num_bands += concordances[path[i]][path[i+1]]['bands']
        return num_bands

    ribbon_bounds = {K:num_bands(K) for K in ribbon}
    to_update = df.name.isin(ribbon_bounds)
    df.loc[to_update, 'ribbon_bound'] = df.name[to_update].map(ribbon_bounds)
    return concordances, components


def check_concord_decreases(dataframe, concordances, numerical_column):
    volumes = dict(zip(dataframe.name, dataframe[numerical_column]))
    diffs, ratios = [], []
    min_diff_pair, max_ratio_pair = [(100000, None, None), (0, None, None)]
    for K1, K2 in concordances.edges:
        if K2 != 'unknot':
            V1, V2 = volumes[K1], volumes[K2]
            if V2 != -1:
                assert V1 > V2
                diff, ratio = V1 - V2, V2/V1
                diffs.append(diff)
                ratios.append(ratio)
                if min_diff_pair[0] > diff:
                    min_diff_pair = (diff, K1, K2, V1, V2)
                if ratio > max_ratio_pair[0]:
                    max_ratio_pair = (ratio, K1, K2, V1, V2)

    return min_diff_pair, max_ratio_pair


def check_no_concord_alt_to_non(dataframe, concordances):
    alt = dict(zip(dataframe.name, dataframe.alt))
    for K1, K2 in concordances.edges:
        if K2 != 'unknot':
            if alt[K1] == 1 and alt[K2] == 0:
                raise ValueError('Look for flying pigs')


def check_no_fibered_to_nonfibered(dataframe, concordances):
    fibered = dict(zip(dataframe.name, dataframe.fibered))
    for K1, K2 in concordances.edges:
        if K2 != 'unknot':
            if fibered[K1] == 1 and fibered[K2] == 0:
                raise ValueError(f'Edge {K1} --> {K2} is interesting')


def count_concord_from_alt(dataframe, concordances):
    count = 0
    alt = dict(zip(dataframe.name, dataframe.alt))
    for K1, K2 in concordances.edges:
        if alt[K1] == 1:
            count += 1
    return count


def save_computational_columns(dataframe):
    """
    Save the slice, ribbon, and top_slice columns for later
    comparison.
    """
    for col in ['slice', 'ribbon', 'top_slice']:
        dataframe[col + '_old'] = dataframe[col]


def update_slice_basic(dataframe):
    """
    The slice, ribbon, and top_slice columns are "computational", to
    be filled in using data from the other columns.  In particular,
    this function wipes their contents.
    """
    df = dataframe

    # (column containing invariant, value obstrucing slicing)
    slice_obs = [('s_0', 0), ('s_2', 0), ('s_3', 0), ('s_Z', '((0,),(0,))'),
                 ('Sq1', '(0,0,0,0)'),  ('Sq1_odd', '(0,0,0,0)'), ('CLS_red', '(0,0)'),
                 ('tau', 0), ('nu', 0), ('epsilon', 0),
                 ('sl_3_0', 0), ('sl_3_2', 0), ('sl_3_3', 0)]

    # Remove bifact obs specific to alt knots as is no longer needed.
    not_slice = ((df.HKL_basic.notnull())|
                 (df.HKL_fancy.notnull())|
                 (df.HKL_direct.notnull()))
    for col, val in slice_obs:
        not_slice = not_slice | (df[col].notnull()&(df[col]!=val))

    # Reset cols
    for col in ['slice', 'ribbon']:
        df[col] = 0
        df.loc[not_slice, col] = -1

    known_ribbon = df.ribbon_bound > 0
    for col in  ['slice', 'ribbon']:
         assert all(df[known_ribbon][col] == 0)
         df.loc[known_ribbon, col] = 1

    ### Now the knots that we just know are slice, but not nec. ribbon ###

    concord_slice = df.concord_to.apply(lambda x:isinstance(x, str) and
                                     (x.find('unknot') > -1 or x.find('K6a3') > -1))

    # Results from adding 4_1's
    # raw_output = open('whistler/adding_4_1s_output/adding_4_1_results_summary.txt')
    # more_slice = df.name.isin({line.strip() for line in raw_output.readlines()[::2]})
    just_slice = concord_slice
    # assert all(df[just_slice].slice != -1) and all(df[just_slice].ribbon != -1)

    progress = sum(just_slice&(df.slice!=1))
    assert progress == 0
    print('Originally we only knew slice for', sum(just_slice),
          'knots but all are now known to be ribbon')

    df.loc[just_slice, 'slice'] = 1

    ### Now the knots that we just know *not ribbon* ###

    non_ribbon = df.HKL_ribbon_obs.notnull()
    assert all(df[non_ribbon].ribbon != 1)
    update = non_ribbon & (df.ribbon == 0)
    print('Just nonribbon:', sum(update))
    df.loc[update, 'ribbon'] = -1

    ### Finally, some special cases ###
    #
    #  * K11n34: Lisa's Annals paper.
    #
    #  * 17ns_29: described as Fig8[1] by Burton; should be the (2, 1)
    #    cable on the fig 8 https://arxiv.org/abs/2207.14187
    #
    #  * The knot 19nh_000143796 is handled in maggie/maggie_obs.py
    #
    #  * The other 24 knots are our work, using the Narakmura-Ren
    #     enhanced version of the Piccirillo principle.
    #
    # Of these, only K11n34 and K13n866 are parts of nontrivial
    # components of our ribbon concordance graph, of sizes 1673 and
    # 142, respectively.

    not_smooth_slice  = {'K11n34', '17ns_29',
                         '19nh_000143796',
                         'K13n866', 'K15n25044', '16n74539',
                         '16n180537','17nh_0001844', '17nh_0002715',
                         '17nh_0212094', '17nh_0212095', '18nh_00010270',
                         '18nh_00166702', '18nh_00610378', '18nh_00610381',
                         '19nh_000002588', '19nh_000003154', '19nh_000003570',
                         '19nh_000032808', '19nh_000076489', '19nh_000018991',
                         '19nh_000066839', '19nh_000066841', '19nh_000177115',
                         '19nh_000177116', '19nh_001336127', '19nh_002457201'}

    update = df.name.isin(not_smooth_slice)
    assert all(df.loc[update, 'slice'] == 0)
    assert all(df.loc[update, 'ribbon'] == 0)
    df.loc[update, 'slice'] = -1
    df.loc[update, 'ribbon'] = -1

    # Now we turn to the topological setting

    not_slice = (df.HKL_basic.notnull()|
                 df.HKL_fancy.notnull()|
                 df.HKL_direct.notnull())
    df['top_slice'] = 0
    df.loc[not_slice, 'top_slice'] = -1

    # The last one is from Theorem 5.14.

    known_top_slice = ((df.slice==1)|(df.ribbon==1)|
                       (df.alex=='1')|(df.name=='18nh_00000601'))
    assert all(df[known_top_slice]['top_slice'] == 0)
    df.loc[known_top_slice, 'top_slice'] = 1

    # These are were shown to be nonslice in "special_alt_knots.py"

    not_top_slice = {'K14a12741', '16a350194', '17ah_0000055',
                     '18ah_0000122', '18ah_3327857', '18ah_4025786',
                     '18ah_4099296', '18ah_4099297', '19ah_00000457'}

    update = df.name.isin(not_top_slice)
    for col in ['top_slice', 'slice', 'ribbon']:
        assert all(df.loc[update, col] == 0)
        df.loc[update, col] = -1


def find_extra_top_slice(df, ribbon_digraph, connected_components):
    # Used in user function update_slice_from_concordance.
    HKL_cols = ['HKL_basic', 'HKL_fancy', 'HKL_direct']
    extra_slice = set()
    extra_nonslice = set()
    for C in connected_components:
        dc = df.loc[list(C - {'unknot'})]
        top_slice = set(dc.top_slice)
        if top_slice == {1} or top_slice == {-1}:
            continue

        # Find the root knot
        D = ribbon_digraph.subgraph(C)
        u = list(nx.topological_sort(D))[-1]
        ancestors = nx.ancestors(D, u)
        assert len(ancestors) + 1 == D.number_of_nodes()

        if top_slice == {0, 1}:
            new_slice = set(dc.index[dc.top_slice==0])
            extra_slice.update(new_slice)
            print(f'Slice comp of size {len(C)} with sink {u}')
            if new_slice:
                print(f'   Gained {len(new_slice)} topologically slice knots')
                cause = dc[dc.top_slice==1].iloc[0]
                print(f"   Due to {cause.name} where Alex is {cause.alex}")
        elif top_slice == {0, -1}:
            new_nonslice = set(dc.index[dc.top_slice==0])
            extra_nonslice.update(new_nonslice)
            print(f'Nonslice comp of size {len(C)} with sink {u}')
            if new_nonslice:
                print(f'   Gained {len(new_nonslice)} topologically nonslice knots')
                cause = dc[dc.top_slice==-1].iloc[0]
                HKL = [cause[col] for col in HKL_cols]
                print(f"   Due to {cause.name} where HKL is {HKL}")
        elif top_slice == {0}:
            print('Unknown comp of size', len(C), 'with sink', u)

    return extra_slice, extra_nonslice

    # Put new info into main dataframe
    extra_nonslice.update(known_extra_not_slice)
    to_update = dataframe.name.isin(extra_nonslice)
    dataframe.loc[to_update, 'top_slice'] = -1
    to_update = dataframe.name.isin(extra_slice)
    dataframe.loc[to_update, 'top_slice'] = 1
    return extra_slice, extra_nonslice


def find_extra_smooth_slice(df, ribbon_digraph, connected_components):
    # Used in user function update_slice_from_concordance.
    extra_nonslice = set()
    for C in connected_components:
        dc = df.loc[list(C - {'unknot'})]
        slice = set(dc.slice)
        ribbon = set(dc.ribbon)
        if slice == ribbon == {1} or slice == ribbon == {-1}:
            continue

        # We have a few non-ribbon knots of unknown slicing, but the
        # only components where they appear are of this form:
        if slice != ribbon:
            assert slice == {0, -1} and ribbon == {-1}

        # Find the root knot
        D = ribbon_digraph.subgraph(C)
        u = list(nx.topological_sort(D))[-1]
        ancestors = nx.ancestors(D, u)
        assert len(ancestors) + 1 == D.number_of_nodes()

        # Do deductions
        if slice == {0, 1}:
            raise ValueError('Problem with slice component')
        elif slice == {0, -1}:
            new_nonslice = set(dc.index[dc.slice==0])
            extra_nonslice.update(new_nonslice)
            print(f'Nonslice comp of size {len(C)} with sink {u}')
            if new_nonslice:
                print(f'   Gained {len(new_nonslice)} smoothly nonslice knots')
        elif slice == {0}:
            print('Unknown comp of size', len(C), 'with sink', u)
        else:
            raise ValueError

    return extra_nonslice


def update_slice_from_concordances(dataframe, ribbon_digraph, connected_components=None):
    HKL_cols = ['HKL_basic', 'HKL_fancy', 'HKL_direct']
    cols = ['name', 'slice', 'top_slice', 'ribbon', 'alex'] + HKL_cols
    df = dataframe[cols].set_index('name', drop=True)
    if connected_components is None:
        G = nx.Graph(ribbon_digraph)
        connected_components = nx.connected_components(G)

    extra_top_slice, extra_top_nonslice = find_extra_top_slice(df, ribbon_digraph, connected_components)

    # Put new info into main dataframe
    to_update = dataframe.name.isin(extra_top_nonslice)
    for col in ['top_slice', 'slice', 'ribbon']:
        dataframe.loc[to_update, col] = -1

    to_update = dataframe.name.isin(extra_top_slice)
    dataframe.loc[to_update, 'top_slice'] = 1
    df = dataframe[cols].set_index('name', drop=True)

    extra_nonslice = find_extra_smooth_slice(df, ribbon_digraph, connected_components)
    to_update = dataframe.name.isin(extra_nonslice)
    dataframe.loc[to_update, 'slice'] = -1
    dataframe.loc[to_update, 'ribbon'] = -1
    return extra_top_slice, extra_top_nonslice, extra_nonslice


def main():
    """

    Load the database, build the ribbon concordance graph, and check
    everything matches.

    """
    print('Loading data from file...')
    df = pd.read_csv('../data/plausibly_slice.csv.bz2')

    print('Saved slice, ribbon, and top_slice columns for comparison at the end.')
    save_computational_columns(df)

    print('\nComputing ribbon concordance graph and its components...')
    rib_concord_graph, components = examine_ribbon(df)

    print('\nStart with the direct information...')
    update_slice_basic(df)

    print('\n Propogating information...')
    updated = update_slice_from_concordances(df, rib_concord_graph, components)

    print('\n Checking results against the data in the file...')
    for col in ['slice', 'ribbon', 'top_slice']:
        matches = all(df[col] == df[col + '_old'])
        print(f'    {col} matches: {matches}')
        
            

          

    
    

if __name__ == '__main__':
    main()
