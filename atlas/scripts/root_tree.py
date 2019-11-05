import sys
import ete3


sys.stdout= open(snakemake.log[0],"w")

T= ete3.Tree(snakemake.input.tree,quoted_node_names=True,format=1)
T.unroot()
T.set_outgroup(T.get_midpoint_outgroup())
T.write(outfile=snakemake.output.tree)
