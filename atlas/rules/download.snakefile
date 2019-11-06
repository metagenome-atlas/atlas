import hashlib
import os


ZENODO_ARCHIVE = "1134890"
EGGNOG_VERSION = "5"


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
EGGNOG_DIR = os.path.join(DBDIR,'EggNOGV2')
GTDBTK_DATA_PATH=os.path.join(DBDIR,"GTDB-TK")
CONDAENV = "../envs"

# note: saving OG_fasta.tar.gz in order to not create secondary "success" file
FILES = {"adapters.fa": "ae839dc79cfb855a1b750a0d593fe01e",
         "phiX174_virus.fa": "82516880142e8c89b466bc6118696c47",
         "refseq.db": "42b8976656f2cfd661b8a299d6e24c19",
         "refseq.dmnd": "c01facc7e397270ccb796ea799a09108",
         "refseq.tree": "469fcbeb15dd0d4bf8f1677682bde157",
         "silva_rfam_all_rRNAs.fa": "f102e35d9f48eabeb0efe9058559bc66",
         "OG_fasta": "8fc6ce2e055d1735dec654af98a641a4",
         "eggnog.db": "e743ba1dbc3ddc238fdcc8028968aacb",
         "eggnog_proteins.dmnd": "5efb0eb18ed4575a20d25773092b83b9",
         "og2level.tsv": "d35ffcc533c6e12be5ee8e5fd7503b84",
         CHECKM_ARCHIVE: "631012fa598c43fdeb88c619ad282c4d"}


CHECKMFILES=[   "%s/taxon_marker_sets.tsv" % CHECKMDIR,
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
        ]



localrules: download,download_eggNOG_files,download_atlas_files,unpack_checkm_data
ruleorder: download_eggNOG_fastas > download_eggNOG_files > download_atlas_files

rule download:
    input:
        expand("{dir}/{filename}", dir=DBDIR,
               filename=["adapters.fa","phiX174_virus.fa"]),
        expand("{dir}/{filename}", dir=EGGNOG_DIR,
               filename=["eggnog.db","eggnog_proteins.dmnd"]),
        CHECKMFILES,
        os.path.join(GTDBTK_DATA_PATH,'downloaded_success')



rule download_eggNOG_files:
    output:
        f"{EGGNOG_DIR}/eggnog.db",
        f"{EGGNOG_DIR}/eggnog_proteins.dmnd"
    threads:
        1
    run:
        shell(f"download_eggnog_data.py -yf --data_dir {EGGNOG_DIR} " )
        # validate the download
        for file in input:
            if not FILES[os.path.basename(file)] == md5(file):
                raise OSError(2, "Invalid checksum", file)
        # check if old eggNOG dir exists
        old_eggnogdir=EGGNOG_DIR.replace('V2','')
        if os.path.exists(old_eggnogdir):
            logger.info("The eggnog database form the olf v1 was found on your system."
                        f"You can savely remove this folder {old_eggnogdir}")




rule download_atlas_files:
    output:
        f"{DBDIR}/{{filename}}"
    threads:
        1
    run:
        shell("wget -O {output} 'https://zenodo.org/record/{ZENODO_ARCHIVE}/files/{wildcards.filename}' ")
        if not FILES[wildcards.filename] == md5(output[0]):
            raise OSError(2, "Invalid checksum", output[0])


rule unpack_checkm_data:
    input:
        os.path.join(DBDIR, CHECKM_ARCHIVE)
    output:
        CHECKMFILES
    params:
        path = CHECKMDIR
    shell:
        "tar -zxf {input} --directory {params.path}"

localrules: initialize_checkm
rule initialize_checkm:
    input:
        ancient(CHECKMFILES)
    output:
        touched_output = touch("logs/checkm_init.txt")
    params:
        database_dir = CHECKMDIR,
    conda:
        "%s/checkm.yaml" % CONDAENV
    log:
        "logs/initialize_checkm.log"
    shell:
        "checkm data setRoot {params.database_dir} &> {log} "

rule download_cat_db:
    output:
        touch(CAT_flag_downloaded)
    params:
        db_folder=CAT_DIR
    resources:
        mem= config.get('large_mem',250)
    threads:
        config.get('large_threads',16)
    conda:
        "%s/cat.yaml" % CONDAENV
    shell:
        " CAT prepare -d {params.db_folder} -t {params.db_folder} --existing --nproc {threads}"



rule download_gtdb:
    output:
        touch(os.path.join(GTDBTK_DATA_PATH,'downloaded_success'))
    conda:
        "../envs/gtdbtk.yaml"
    threads:
        1
    shell:
        "GTDBTK_DATA_PATH={GTDBTK_DATA_PATH} ;  "
        "download-db.sh ;"

onsuccess:
    print("All databases have downloaded and validated successfully")


onerror:
    print("An error occurred while downloading reference databases.")
    print("ATLAS databases can be manually downloaded from: https://zenodo.org/record/%s" % ZENODO_ARCHIVE)
    print("eggNOG databases can be manually downloaded from: http://eggnogdb.embl.de/download/emapperdb-%s" % EGGNOG_VERSION)
    print("CAT databases can be manually downloaded from: https://github.com/dutilh/CAT")
