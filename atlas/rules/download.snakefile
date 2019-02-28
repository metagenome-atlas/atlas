import hashlib
import os


ZENODO_ARCHIVE = "1134890"
EGGNOG_VERSION = "4.5.1"


def md5(fname):
    # https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
    hash_md5 = hashlib.md5()
    if not os.path.exists(fname):
        return None
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()




# this values are incuded in the snakefile
DBDIR = os.path.realpath(config["database_dir"])
CHECKMDIR = os.path.join(DBDIR, "checkm")
CHECKM_ARCHIVE = "checkm_data_v1.0.9.tar.gz"
CAT_DIR= os.path.join(DBDIR,'CAT')
CAT_flag_downloaded = os.path.join(CAT_DIR,'downloaded')
EGGNOG_DIR = DBDIR

CONDAENV = "../envs"

# note: saving OG_fasta.tar.gz in order to not create secondary "success" file
FILES = {"adapters.fa": "ae839dc79cfb855a1b750a0d593fe01e",
         "phiX174_virus.fa": "82516880142e8c89b466bc6118696c47",
         "refseq.db": "42b8976656f2cfd661b8a299d6e24c19",
         #"refseq.dmnd": "c01facc7e397270ccb796ea799a09108",
         #"refseq.tree": "469fcbeb15dd0d4bf8f1677682bde157",
         "silva_rfam_all_rRNAs.fa": "f102e35d9f48eabeb0efe9058559bc66",
         "OG_fasta": "8fc6ce2e055d1735dec654af98a641a4",
         "eggnog.db": "e743ba1dbc3ddc238fdcc8028968aacb",
         "eggnog_proteins.dmnd": "5efb0eb18ed4575a20d25773092b83b9",
         "og2level.tsv": "d35ffcc533c6e12be5ee8e5fd7503b84",
         CHECKM_ARCHIVE: "631012fa598c43fdeb88c619ad282c4d"}

localrules: download, transfer_files

rule download:
    input:
        expand("{dir}/{filename}", dir=DBDIR, filename=list(FILES.keys())),
        "%s/taxon_marker_sets.tsv" % CHECKMDIR





rule transfer_files:
    output:
        "%s/{filename}" % DBDIR
    run:
        eggnog_files = {"OG_fasta": "OG_fasta.tar.gz",
            "eggnog.db": "eggnog.db.gz",
            "og2level.tsv": "og2level.tsv.gz",
            "eggnog_proteins.dmnd": "eggnog_proteins.dmnd.gz",
        }

# OG_fasta.tar.gz will be unzipped to a directory !!

        if wildcards.filename in eggnog_files:
            dl_filename = eggnog_files[wildcards.filename]
            dl_output = os.path.join(os.path.dirname(output[0]), dl_filename)
            shell("curl 'http://eggnogdb.embl.de/download/emapperdb-%s/%s' -s > %s" % (EGGNOG_VERSION, dl_filename, dl_output))
            # validate the download
            if not FILES[wildcards.filename] == md5(dl_output):
                raise OSError(2, "Invalid checksum", dl_output)
            # handle extraction/decompression
            if dl_output.endswith(".tar.gz"):
                shell("tar -zxf %s --directory {DBDIR}" % dl_output)
            else:
                shell("gunzip %s" % dl_output)
        else:
            shell("curl 'https://zenodo.org/record/%s/files/{wildcards.filename}' -s > {output}" % ZENODO_ARCHIVE)
            if not FILES[wildcards.filename] == md5(output[0]):
                raise OSError(2, "Invalid checksum", output[0])


rule download_checkm_data:
    input:
        "%s/%s" % (DBDIR, CHECKM_ARCHIVE)
    output:
        "%s/taxon_marker_sets.tsv" % CHECKMDIR,
        "%s/selected_marker_sets.tsv" % CHECKMDIR,
        "%s/pfam/tigrfam2pfam.tsv" % CHECKMDIR,
        "%s/pfam/Pfam-A.hmm.dat" % CHECKMDIR,
        "%s/img/img_metadata.tsv" % CHECKMDIR,
        "%s/hmms_ssu/SSU_euk.hmm" % CHECKMDIR,
        "%s/hmms_ssu/SSU_bacteria.hmm" % CHECKMDIR,
        "%s/hmms_ssu/SSU_archaea.hmm" % CHECKMDIR,
        "%s/hmms_ssu/createHMMs.py" % CHECKMDIR,
        "%s/hmms/phylo.hmm.ssi" % CHECKMDIR,
        "%s/hmms/phylo.hmm" % CHECKMDIR,
        "%s/hmms/checkm.hmm.ssi" % CHECKMDIR,
        "%s/hmms/checkm.hmm" % CHECKMDIR,
        "%s/genome_tree/missing_duplicate_genes_97.tsv" % CHECKMDIR,
        "%s/genome_tree/missing_duplicate_genes_50.tsv" % CHECKMDIR,
        "%s/genome_tree/genome_tree.taxonomy.tsv" % CHECKMDIR,
        "%s/genome_tree/genome_tree_reduced.refpkg/phylo_modelJqWx6_.json" % CHECKMDIR,
        "%s/genome_tree/genome_tree_reduced.refpkg/genome_tree.tre" % CHECKMDIR,
        "%s/genome_tree/genome_tree_reduced.refpkg/genome_tree.log" % CHECKMDIR,
        "%s/genome_tree/genome_tree_reduced.refpkg/genome_tree.fasta" % CHECKMDIR,
        "%s/genome_tree/genome_tree_reduced.refpkg/CONTENTS.json" % CHECKMDIR,
        "%s/genome_tree/genome_tree.metadata.tsv" % CHECKMDIR,
        "%s/genome_tree/genome_tree_full.refpkg/phylo_modelEcOyPk.json" % CHECKMDIR,
        "%s/genome_tree/genome_tree_full.refpkg/genome_tree.tre" % CHECKMDIR,
        "%s/genome_tree/genome_tree_full.refpkg/genome_tree.log" % CHECKMDIR,
        "%s/genome_tree/genome_tree_full.refpkg/genome_tree.fasta" % CHECKMDIR,
        "%s/genome_tree/genome_tree_full.refpkg/CONTENTS.json" % CHECKMDIR,
        "%s/genome_tree/genome_tree.derep.txt" % CHECKMDIR,
        "%s/.dmanifest" % CHECKMDIR,
        "%s/distributions/td_dist.txt" % CHECKMDIR,
        "%s/distributions/gc_dist.txt" % CHECKMDIR,
        "%s/distributions/cd_dist.txt" % CHECKMDIR
    params:
        path = CHECKMDIR
    shell:
        "tar -zxf {input} --directory {params.path}"


rule download_cat_db:
    output:
        touch(CAT_flag_downloaded)
    params:
        db_folder=CAT_DIR
    resources:
        mem= config.get('diamond_mem',10)
    threads:
        config.get('diamond_threads',10)
    conda:
        "%s/cat.yaml" % CONDAENV
    shell:
        " CAT prepare -d {params.db_folder} -t {params.db_folder} --existing --nproc {threads}"


onsuccess:
    print("All databases have downloaded and validated successfully")


onerror:
    print("An error occurred while downloading reference databases.")
    print("ATLAS databases can be manually downloaded from: https://zenodo.org/record/%s" % ZENODO_ARCHIVE)
    print("eggNOG databases can be manually downloaded from: http://eggnogdb.embl.de/download/emapperdb-%s" % EGGNOG_VERSION)
