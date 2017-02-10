import click
import logging
import multiprocessing
import os
from atlas import __version__
from atlas.conf import make_config
from atlas.parsers import cazy_parser, cog_parser, eggnog_parser, enzyme_parser, refseq_parser
from atlas.tables import merge_tables, count_tables
from atlas.workflows import assemble, download


logging.basicConfig(level=logging.INFO, datefmt="%Y-%m-%d %H:%M", format="[%(asctime)s] %(message)s")


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    """ATLAS"""


@cli.command("cazy", short_help="process blast hits for cazy (dbcan) reference")
@click.argument("tsv", type=click.Path(exists=True))
@click.argument("namemap", type=click.Path(exists=True))
@click.argument("output", type=click.File("w", atomic=True))
@click.option("-s", "--summary-method", type=click.Choice(["majority", "best"]), default="best", show_default=True, help="summary method for annotating ORFs; when majority and there is no majority, best is used")
@click.option("--min-identity", type=int, default=60, show_default=True, help="minimum allowable percent ID of BLAST hit")
@click.option("--min-bitscore", type=int, default=0, show_default=True, help="minimum allowable bitscore of BLAST hit; 0 disables")
@click.option("--min-length", type=int, default=60, show_default=True, help="minimum allowable BLAST alignment length")
@click.option("--max-evalue", type=float, default=0.000001, show_default=True, help="maximum allowable e-value of BLAST hit")
@click.option("--top-fraction", type=float, default=1, show_default=True, help="filters ORF BLAST hits before finding majority by only keep hits within this fraction, e.g. 0.98, of the highest bitscore; this is recommended over --max-hits")
@click.option("--max-hits", type=int, default=10, show_default=True, help="maximum number of BLAST hits to consider when summarizing ORFs as a majority")
@click.option("--table-name", default="dbcan", help="table name within namemap database; expected columns are listed above")
def run_cazy_parser(tsv, namemap, output, summary_method, min_identity, min_bitscore, min_length, max_evalue, top_fraction, max_hits, table_name):
    """Parse BLAST hits from CAZy reference database.

    The BLAST hits are assumed to be sorted by query with decreasing bitscores (best alignment first):

        \b
        sort -k1,1 -k12,12rn tsv > sorted_tsv

    Expected columns in the CAZy database:

        \b
        cazy_gene
        cazy_family
        cazy_class
        cazy_ec

    For a given gene match, all possible ECs are returned in a single line separated by '|'.

    """
    cazy_parser(tsv, namemap, output, summary_method, min_identity, min_bitscore, min_length, max_evalue, top_fraction, max_hits, table_name)


@cli.command("cog", short_help="process blast hits for COG reference")
@click.argument("tsv", type=click.Path(exists=True))
@click.argument("namemap", type=click.Path(exists=True))
@click.argument("output", type=click.File("w", atomic=True))
@click.option("-s", "--summary-method", type=click.Choice(["majority", "best"]), default="best", show_default=True, help="summary method for annotating ORFs; when majority and there is no majority, best is used")
@click.option("--min-identity", type=int, default=60, show_default=True, help="minimum allowable percent ID of BLAST hit")
@click.option("--min-bitscore", type=int, default=0, show_default=True, help="minimum allowable bitscore of BLAST hit; 0 disables")
@click.option("--min-length", type=int, default=60, show_default=True, help="minimum allowable BLAST alignment length")
@click.option("--max-evalue", type=float, default=0.000001, show_default=True, help="maximum allowable e-value of BLAST hit")
@click.option("--top-fraction", type=float, default=1, show_default=True, help="filters ORF BLAST hits before finding majority by only keep hits within this fraction, e.g. 0.98, of the highest bitscore; this is recommended over --max-hits")
@click.option("--max-hits", type=int, default=10, show_default=True, help="maximum number of BLAST hits to consider when summarizing ORFs as a majority")
@click.option("--table-name", default="cog", help="table name within namemap database; expected columns are listed above")
def run_cog_parser(tsv, namemap, output, summary_method='best', min_identity=60, min_bitscore=0, min_length=60, max_evalue=0.000001, top_fraction=1, max_hits=10, table_name="cog"):
    """Parse BLAST hits from COG reference database.

    The BLAST hits are assumed to be sorted by query with decreasing bitscores (best alignment first):

        \b
        sort -k1,1 -k12,12rn tsv > sorted_tsv

    Expected columns in the COG database:

        \b
        cog_protein_key (unique and matches sequence name in reference fasta)
        cog_protein_id
        cog_id
        cog_functional_class
        cog_annotation
        cog_functional_class_description
    """
    cog_parser(tsv, namemap, output, summary_method, min_identity, min_bitscore, min_length, max_evalue, top_fraction, max_hits, table_name)


@cli.command("eggnog", short_help="process blast hits for eggnog reference")
@click.argument("tsv", type=click.Path(exists=True))
@click.argument("namemap", type=click.Path(exists=True))
@click.argument("output", type=click.File("w", atomic=True))
@click.option("-s", "--summary-method", type=click.Choice(["majority", "best"]), default="best", show_default=True, help="summary method for annotating ORFs; when majority and there is no majority, best is used")
@click.option("--min-identity", type=int, default=60, show_default=True, help="minimum allowable percent ID of BLAST hit")
@click.option("--min-bitscore", type=int, default=0, show_default=True, help="minimum allowable bitscore of BLAST hit; 0 disables")
@click.option("--min-length", type=int, default=60, show_default=True, help="minimum allowable BLAST alignment length")
@click.option("--max-evalue", type=float, default=0.000001, show_default=True, help="maximum allowable e-value of BLAST hit")
@click.option("--top-fraction", type=float, default=1, show_default=True, help="filters ORF BLAST hits before finding majority by only keep hits within this fraction, e.g. 0.98, of the highest bitscore; this is recommended over --max-hits")
@click.option("--max-hits", type=int, default=10, show_default=True, help="maximum number of BLAST hits to consider when summarizing ORFs as a majority")
@click.option("--table-name", default="eggnog", help="table name within namemap database; expected columns are listed above")
def run_eggnog_parser(tsv, namemap, output, summary_method, min_identity, min_bitscore, min_length, max_evalue, top_fraction, max_hits, table_name):
    """Parse BLAST hits from EGGNOG.

    The BLAST hits are assumed to be sorted by query with decreasing bitscores (best alignment first):

        \b
        sort -k1,1 -k12,12rn tsv > sorted_tsv

    Expected columns in the EggNOG database:

        \b
        uniprot_ac
        eggnog_ssid_b
        eggnog_species_id
        uniprot_id
        cog_func_id
        cog_id
        cog_product
        cog_level1_code
        cog_level1_name
        cog_level2_name
        ko_id
        ko_level1_name
        ko_level2_name
        ko_level3_id
        ko_level3_name
        ko_gene_symbol
        ko_product
        ko_ec

    """
    eggnog_parser(tsv, namemap, output, summary_method, min_identity, min_bitscore, min_length, max_evalue, top_fraction, max_hits, table_name)


@cli.command("enzyme", short_help="process blast hits for ENZYME reference")
@click.argument("tsv", type=click.Path(exists=True))
@click.argument("namemap", type=click.Path(exists=True))
@click.argument("output", type=click.File("w", atomic=True))
@click.option("-s", "--summary-method", type=click.Choice(["majority", "best"]), default="best", show_default=True, help="summary method for annotating ORFs; when majority and there is no majority, best is used")
@click.option("--min-identity", type=int, default=60, show_default=True, help="minimum allowable percent ID of BLAST hit")
@click.option("--min-bitscore", type=int, default=0, show_default=True, help="minimum allowable bitscore of BLAST hit; 0 disables")
@click.option("--min-length", type=int, default=60, show_default=True, help="minimum allowable BLAST alignment length")
@click.option("--max-evalue", type=float, default=0.000001, show_default=True, help="maximum allowable e-value of BLAST hit")
@click.option("--top-fraction", type=float, default=1, show_default=True, help="filters ORF BLAST hits before finding majority by only keep hits within this fraction, e.g. 0.98, of the highest bitscore; this is recommended over --max-hits")
@click.option("--max-hits", type=int, default=10, show_default=True, help="maximum number of BLAST hits to consider when summarizing ORFs as a majority")
@click.option("--table-name", default="enzyme", help="table name within namemap database; expected columns are listed above")
def run_enzyme_parser(tsv, namemap, output, summary_method, min_identity, min_bitscore, min_length, max_evalue, top_fraction, max_hits, table_name):
    """Parse BLAST hits from ENZYME reference database.

    The BLAST hits are assumed to be sorted by query with decreasing bitscores (best alignment first):

        \b
        sort -k1,1 -k12,12rn tsv > sorted_tsv

    Expected columns in the ENZYME database:

        \b
        uniprot_entry
        uniparc_entry
        enzyme_ec
        enzyme_name

    For a given UniParc match, all possible ECs and ENZYME recommended names are returned in a
    single line separated by '|'.

    """
    enzyme_parser(tsv, namemap, output, summary_method, min_identity, min_bitscore, min_length, max_evalue, top_fraction, max_hits, table_name)


@cli.command("refseq", short_help="enables tree based LCA and LCA star methods")
@click.argument("tsv", type=click.Path(exists=True))
@click.argument("namemap", type=click.Path(exists=True))
@click.argument("treefile", type=click.Path(exists=True))
@click.argument("output", type=click.File("w", atomic=True))
@click.option("-s", "--summary-method", type=click.Choice(["lca", "majority", "best"]), default="lca", show_default=True, help="summary method for annotating ORFs; when using LCA, it's recommended that one limits the number of hits using --top-fraction though function will be assigned per the best hit; 'best' is fastest")
@click.option("-a", "--aggregation-method", type=click.Choice(["lca", "lca-majority", "majority"]), default="lca-majority", show_default=True, help="summary method for aggregating ORF taxonomic assignments to contig level assignment; 'lca' will result in most stringent, least specific assignments")
@click.option("--majority-threshold", type=float, default=0.51, show_default=True, help="constitutes a majority fraction at tree node for 'lca-majority' ORF aggregation method")
@click.option("--min-identity", type=int, default=60, show_default=True, help="minimum allowable percent ID of BLAST hit")
@click.option("--min-bitscore", type=int, default=0, show_default=True, help="minimum allowable bitscore of BLAST hit; 0 disables")
@click.option("--min-length", type=int, default=60, show_default=True, help="minimum allowable BLAST alignment length")
@click.option("--max-evalue", type=float, default=0.000001, show_default=True, help="maximum allowable e-value of BLAST hit")
@click.option("--max-hits", type=int, default=10, show_default=True, help="maximum number of BLAST hits to consider when summarizing ORFs; can drastically alter ORF LCA assignments if too high without further limits")
@click.option("--table-name", default="refseq", help="table name within namemap database; expected columns are 'name', 'function', and 'taxonomy'")
@click.option("--top-fraction", type=float, default=1, show_default=True, help="filters ORF BLAST hits by only keep hits within this fraction of the highest bitscore; this is recommended over --max-hits")
def run_refseq_parser(tsv, namemap, treefile, output, summary_method, aggregation_method, majority_threshold, min_identity, min_bitscore, min_length, max_evalue, max_hits, table_name, top_fraction):
    """Parse TSV (tabular BLAST output [-outfmt 6]), grabbing taxonomy metadata from ANNOTATION to
    compute LCAs.

    The BLAST hits are assumed to be sorted by query with decreasing bitscores (best alignment first):

        \b
        sort -k1,1 -k12,12rn tsv > sorted_tsv

    Annotation file should include your BLAST subject sequence ID, a function, a taxonomy name,
    the taxonomy ID, and the parent taxonomy ID. This file is generated from `prepare-refseq`:

        \b
        gi|507051347|ref|WP_016122340.1| two-component sensor histidine kinase Bacillus cereus 1396 86661
        gi|507052147|ref|WP_016123132.1| two-component sensor histidine kinase Bacillus cereus 1396 86661
        gi|507053266|ref|WP_016124222.1| two-component sensor histidine kinase Bacillus cereus 1396 86661

    The RefSeq function is always derived from the best BLAST hit.

    The output will give contig, ORF ID, the lineage assigned to the contig based on
    --aggregation-method, the probability of error (erfc), taxonomy assigned to the ORF, the
    best hit's product, the best hit's evalue, and the best hit's bitscore:

        \b
        contig     orf          taxonomy erfc orf_taxonomy refseq               refseq_evalue refseq_bitscore
        k121_52126 k121_52126_1 root     1.0  root         hypothetical protein 1.0e-41       176.8

    """
    refseq_parser(tsv, namemap, treefile, output, summary_method, aggregation_method, majority_threshold, min_identity, min_bitscore, min_length, max_evalue, max_hits, table_name, top_fraction)


@cli.command("merge-tables", short_help="merge tables on 'contig' and 'orf' keys")
@click.argument("tables", type=click.File("r"), nargs=-1)
@click.argument("output", type=click.File("w", atomic=True))
def run_merge_tables(tables, output):
    """Takes the output from parsers and combines them into a single TSV table.

    Headers are required and should contain 'contig' and 'orf' column labels.
    """
    merge_tables(tables, output)


@cli.command("counts", short_help="aggregate read counts across annotation combinations")
@click.argument("prefix")
@click.argument("merged", click.Path(exists=True, dir_okay=False, resolve_path=True))
@click.argument("counts", click.Path(exists=True, dir_okay=False, resolve_path=True))
@click.argument("combinations")
@click.option("--suffix", default=".tsv", show_default=True, help="output file suffix")
def run_counts(prefix, merged, counts, combinations, suffix=".tsv"):
    """Aggregate and integrate count data from `counts` with annotation data in `merged`. The
    merged data is the result of `merge-tables`. Count data is a TSV formatted with a header:

        \b
        gene  count
        orf1  10
        orf2  7
        orf3  9

    `combinations` are specified as a JSON string with key to values pairs, e.g.:

        \b
        '{"KO":["ko_id", "ko_gene_symbol", "ko_product", "ko_ec"], "KO_Product":["ko_product"]}'

    Counts get aggregated (summed) across all values, such that the above example gives two files:

        \b
        <prefix>_KO.tsv

            \b
            ko_id   ko_gene_symbol  ko_product                         ko_ec      count
            K00784  rnz             ribonuclease Z                     3.1.26.11  72
            K01006  ppdK            pyruvate, orthophosphate dikinase  2.7.9.1    177
            K01187  malZ            alpha-glucosidase                  3.2.1.20   91

        \b
        <prefix>_KO_product.tsv

            \b
            ko_product                  count
            alpha-glucosidase           91
            beta-galactosidase          267
            cell division protein FtsQ  8
    """
    count_tables(prefix, merged, counts, combinations, suffix)


@cli.command("make-config", short_help="prepopulate a configuration file with samples and defaults")
@click.argument("config")
@click.argument("path")
@click.option("--data-type", default="metagenome", type=click.Choice(["metagenome", "metatranscriptome"]),
              show_default=True, help="sample data type")
@click.option("--database-dir", default="databases", show_default=True,
              help="location of formatted databases (from `atlas download`)")
@click.option("--threads", default=None, type=int,
              help="number of threads to use per multi-threaded job")
@click.option("--assembler", default="megahit", type=click.Choice(["megahit", "spades"]),
              show_default=True, help="contig assembler")
def run_make_config(config, path, data_type, database_dir, threads, assembler):
    """Write the file `config` and complete the sample names and paths for all FASTQ files in
    `path`.

    `path` is traversed recursively and adds any file with '.fastq' or '.fq' extension with the
    file name as the sample ID. Any single-end (non-interleaved) FASTQs under `path` will cause
    errors if left in the configuration file.
    """
    make_config(config, path, data_type, database_dir, threads, assembler)


@cli.command("assemble", context_settings=dict(ignore_unknown_options=True), short_help="assembly workflow")
@click.argument("config")
@click.option("-j", "--jobs", default=multiprocessing.cpu_count(), type=int, show_default=True, help="use at most this many cores in parallel; total running tasks at any given time will be jobs/threads")
@click.option("-o", "--out-dir", default=os.path.realpath("."), show_default=True, help="results output directory")
@click.option("--dryrun", is_flag=True, default=False, show_default=True, help="do not execute anything")
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def run_assemble(config, jobs, out_dir, dryrun, snakemake_args):
    assemble(os.path.realpath(config), jobs, out_dir, dryrun, snakemake_args)


@cli.command("download", context_settings=dict(ignore_unknown_options=True), short_help="download reference files")
@click.option("-j", "--jobs", default=multiprocessing.cpu_count(), type=int, show_default=True, help="use at most this many cores in parallel; total running tasks at any given time will be jobs/threads")
@click.option("-o", "--out-dir", default=os.path.join(os.path.realpath("."), "databases"), show_default=True, help="database download directory")
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def run_download(jobs, out_dir, snakemake_args):
    download(jobs, out_dir, snakemake_args)


if __name__ == "__main__":
    cli()
