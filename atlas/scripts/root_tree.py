import sys
import ete3

with open(snakemake.log[0], "w") as log:
    sys.stderr = sys.stdout = log

    T = ete3.Tree(snakemake.input.tree, quoted_node_names=True, format=1)
    T.unroot()
    if len(T) > 2:
        T.set_outgroup(T.get_midpoint_outgroup())
    T.write(outfile=snakemake.output.tree)
