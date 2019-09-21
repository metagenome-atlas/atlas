import os
os.environ['QT_QPA_PLATFORM']='offscreen' # because we might not have a X server

import ete3
import pandas as pd
import warnings

#from . import parsers_checkm
def load_checkm_tax(checkm_taxonomy_file):

    checkmTax= pd.read_table(checkm_taxonomy_file,index_col=0)

    checkmTax = checkmTax['Taxonomy (contained)']

    if checkmTax.isnull().any():
        warnings.warn("Some samples have no taxonomy asigned based on checkm. Samples:\n"+ \
                    ', '.join(checkmTax.index[checkmTax.isnull()])
                    )
        checkmTax= checkmTax.dropna().astype(str)

    checkmTax= pd.DataFrame(list(  checkmTax.apply(lambda s: s.split(';'))),
                       index=checkmTax.index)

    checkmTax.columns=['kindom','phylum','class','order','family','genus','species']
    return checkmTax


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
    from ete3 import TreeStyle

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
