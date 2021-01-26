import datetime
import shutil
import sys

with open(snakemake.log[0],'w') as log:

    sys.stderr = sys.stdout = log
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d-%X')

    def get_read_stats(fraction, params_in):

        subfolder = os.path.join(snakemake.params.folder, fraction)
        tmp_file=os.path.join(subfolder,"read_stats.tmp")
        shell("""
                mkdir -p {subfolder} 2> {log}

                reformat.sh {params_in} \
                bhist={subfolder}/base_hist.txt \
                qhist={subfolder}/quality_by_pos.txt \
                lhist={subfolder}/readlength.txt \
                gchist={subfolder}/gc_hist.txt \
                gcbins=auto \
                bqhist={subfolder}/boxplot_quality.txt \
                threads={threads} \
                overwrite=true \
                unpigz=t β
                -Xmx{mem}G \
                2> >(tee -a {log} {tmp_file} )
             """.format(subfolder=subfolder, params_in=params_in, log=log,
                        threads=threads, mem=snakemake.resources.java_mem,tmp_file=tmp_file))

        content = open(tmp_file).read()
        pos = content.find('Input:')
        if pos == -1:
            raise Exception("Didn't find read number in file:\n\n" + content)
        else:

            content[pos:].split()[1:4]
                    # Input:    123 reads   1234 bases
            n_reads, _, n_bases = content[pos:].split()[1:4]

            os.remove(tmp_file)
        return int(n_reads), int(n_bases)


    if PAIRED_END:
        n_reads_pe, n_bases_pe = get_read_stats('pe', "in1={0} in2={1}".format(*input))
        n_reads_pe = n_reads_pe / 2
        headers = ['Sample', 'Step', 'Total_Reads', 'Total_Bases',
                   'Reads_pe', 'Bases_pe', 'Reads_se', 'Bases_se',
                   'Timestamp']

        if os.path.exists(snakemake.params.single_end_file):
            n_reads_se, n_bases_se = get_read_stats('se', "in=" + snakemake.params.single_end_file)
        else:
            n_reads_se, n_bases_se = 0, 0

        values=[n_reads_pe + n_reads_se, n_bases_pe + n_bases_se,
                n_reads_pe, n_bases_pe,
                n_reads_se, n_bases_se]
    else:
        headers= ['Sample', 'Step', 'Total_Reads', 'Total_Bases', 'Reads',
                  'Bases', 'Timestamp']
        values = 2 * get_read_stats('', "in=" + input[0])

    with open(snakemake.output.read_counts, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        f.write('\t'.join([snakemake.wildcards.sample, snakemake.wildcards.step] + [str(v) for v in values] + [timestamp]) + '\n')

    shutil.make_archive(snakemake.params.folder, 'zip', snakemake.params.folder)
    shutil.rmtree(snakemake.params.folder)
