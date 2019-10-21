from utils import tree,taxonomy
import sys



sys.stdout= open(snakemake.log[0],"w")

T= tree.load_tree(snakemake.input.tree)


if snakemake.params.taxonomy=='checkm':
    phyla= taxonomy.load_checkm_tax(snakemake.input.taxonomy).phylum
elif snakemake.params.taxonomy=='gtdb':
    phyla= taxonomy.load_gtdb_tax(snakemake.input.taxonomy).phylum

tree.root_tree_by_phyla(T,phyla)

T.write(outfile=snakemake.output.tree)
