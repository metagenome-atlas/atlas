import os, sys
import logging, traceback

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )


# Install exception handler
sys.excepthook = handle_exception


import pandas as pd
from utils.parsers_bbmap import parse_bbmap_log_file


def parse_map_stats(sample_data, out_tsv):
    stats_df = pd.DataFrame()
    for sample in sample_data.keys():
        df = pd.read_csv(sample_data[sample]["contig_stats"], sep="\t")
        assert df.shape[0] == 1, "Assumed only one row in file {}; found {}".format(
            sample_data[sample]["contig_stats"], df.iloc[0]
        )
        df = df.iloc[0]
        df.name = sample
        genes_df = pd.read_csv(sample_data[sample]["gene_table"], index_col=0, sep="\t")
        df["N_Predicted_Genes"] = genes_df.shape[0]
        used_reads, mapped_reads = parse_bbmap_log_file(
            sample_data[sample]["mapping_log"]
        )
        df["Assembled_Reads"] = mapped_reads
        df["Percent_Assembled_Reads"] = mapped_reads / used_reads * 100

        stats_df = stats_df.append(df)
    stats_df = stats_df.loc[:, ~stats_df.columns.str.startswith("scaf_")]
    stats_df.columns = stats_df.columns.str.replace("ctg_", "")
    stats_df.to_csv(out_tsv, sep="\t")
    return stats_df


def main(samples, contig_stats, gene_tables, mapping_logs, combined_stats):
    sample_data = {}
    for sample in samples:
        sample_data[sample] = {}
        for c_stat in contig_stats:
            # underscore version was for simplified local testing
            # if "%s_" % sample in c_stat:
            if "%s/" % sample in c_stat:
                sample_data[sample]["contig_stats"] = c_stat
        for g_table in gene_tables:
            # if "%s_" % sample in g_table:
            if "%s/" % sample in g_table:
                sample_data[sample]["gene_table"] = g_table
        for mapping_log in mapping_logs:
            # if "%s_" % sample in mapping_log:
            if "%s/" % sample in mapping_log:
                sample_data[sample]["mapping_log"] = mapping_log

    parse_map_stats(sample_data, combined_stats)


if __name__ == "__main__":

    main(
        samples=snakemake.params.samples,
        contig_stats=snakemake.input.contig_stats,
        gene_tables=snakemake.input.gene_tables,
        mapping_logs=snakemake.input.mapping_logs,
        combined_stats=snakemake.output.combined_contig_stats,
    )
