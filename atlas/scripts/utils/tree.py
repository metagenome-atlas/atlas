import os
os.environ['QT_QPA_PLATFORM']='offscreen' # because we might not have a X server

import ete3
from ete3 import TreeStyle

import parsers_checkm


def load_tree(netwik_file):
    return ete3.Tree(netwik_file,quoted_node_names=True,format=1)

def root_tree_by_phyla(T,phyla):
    """ Root the tree next to the least frequent phyla if possible

    """


    Freq_pyla= phyla.value_counts()

    for p in reversed(Freq_pyla.index):
        LCA = T.get_common_ancestor(*tuple(phyla.index[phyla==p].values))

        if not T== LCA:
            T.set_outgroup(LCA)
            print(f"set {p} as outgroup for Tree rooting")
            break


    T.unroot()

def layout_black_circles(node):
    # If node is a leaf
    if node.is_leaf():
        node.img_style["fgcolor"]='k'
    else:
        node.img_style["size"]=0

def render_tree(T,out):

    ts = TreeStyle()
    ts.show_leaf_name= False
    ts.mode = "c"
    ts.scale=200
    ts.show_scale=False

    T.render(out,tree_style=ts,layout=layout_black_circles)


if __name__ == "__main__":


    T= load_tree(snakemake.input.tree)
    phyla= parsers_checkm.load_checkm_tax(snakemake.input.taxonomy).phylum

    root_tree_by_phyla(T,phyla)

    T.write(outfile=snakemake.output.tree)

    render_tree(T,snakemake.output.svg)
    render_tree(T,snakemake.output.pdf)
