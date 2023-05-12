import os

os.environ["QT_QPA_PLATFORM"] = "offscreen"  # because we might not have a X server

import ete3
import pandas as pd
import warnings


def load_tree(netwik_file):
    return ete3.Tree(netwik_file, quoted_node_names=True, format=1)


def root_tree_by_phyla(T, phyla):
    """Root the tree next to the phylum that is as far apart as possible from the other phyla"""
    phylum_LCA = {}

    for p in phyla.unique():
        phylum_LCA[p] = T.get_common_ancestor(*tuple(phyla.index[phyla == p].values))

    Dist = pd.DataFrame()
    for p1, lca1 in phylum_LCA.items():
        for p2, lca2 in phylum_LCA.items():
            Dist.loc[p1, p2] = T.get_distance(lca1, lca2)

    furthest_phylum = Dist.mean().idxmax()
    outgroup = phylum_LCA[furthest_phylum]

    if not outgroup == T:
        T.set_outgroup(outgroup)


def layout_black_circles(node):
    # If node is a leaf
    if node.is_leaf():
        node.img_style["fgcolor"] = "k"
    else:
        node.img_style["size"] = 0


def render_tree(T, out):
    from ete3 import TreeStyle

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.mode = "c"
    ts.scale = 200
    ts.show_scale = False

    T.render(out, tree_style=ts, layout=layout_black_circles)
