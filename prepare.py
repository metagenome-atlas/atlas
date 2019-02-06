import os
import logging
import pandas as pd
from collections import defaultdict


def get_sample_files(path,outfile='samples.tsv'):
    samples = defaultdict(dict)
    seen = set()
    for dir_name, sub_dirs, files in os.walk(os.path.abspath(path)):
        for fname in files:

            if ".fastq" in fname or ".fq" in fname:

                sample_id = fname.split(".fastq")[0].split(".fq")[0]

                sample_id = sample_id.replace("_R1", "").replace("_r1", "").replace("_R2", "").replace("_r2", "")
                sample_id = sample_id.replace("_", "-").replace(" ", "-")

                fq_path = os.path.join(dir_name, fname)

                if fq_path in seen: continue

                if "_R2" in fname or "_r2" in fname:

                    if 'R2' in samples[sample_id]:
                        logging.error(f"Duplicate sample {sample_id} was found after renaming; skipping... \n Samples: \n{samples}")

                    samples[sample_id]['R2'] = fq_path
                else:
                    if 'R1' in samples[sample_id]:
                        logging.error(f"Duplicate sample {sample_id} was found after renaming; skipping... \n Samples: \n{samples}")

                    samples[sample_id]['R1'] = fq_path


    samples= pd.DataFrame(samples).T

    if samples.isna().any().any():
        logging.error(f"Missing files:\n {samples}")

    if os.path.exists(outfile):
        logging.error(f"Output file {outfile} already exists I don't date to overwrite it.")
    else:
        samples.to_csv(outfile,sep='\t')


    return samples


if __name__ == '__main__':
    import sys
    get_sample_files(sys.argv[1])
