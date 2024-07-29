import hashlib
import os
from pathlib import Path

# this values are included in the snakefile
DBDIR = Path(config["database_dir"]).resolve()

GUNCDIR = DBDIR / "gunc_database"
BUSCODIR = DBDIR / "busco_lineages"

ZENODO_ARCHIVE = "1134890"
EGGNOG_VERSION = "5"
EGGNOG_DIR = DBDIR / ("EggNOG_V" + EGGNOG_VERSION)

CONDAENV = "../envs"


## GTDBTk

GTDB_VERSION = "V09_R220"
GTDB_DATA_URL = "https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package"
GTDBTK_DATA_PATH = DBDIR / ("GTDB_" + GTDB_VERSION)


def all_partial_gtdb_tarbals(
    wildcards,
    GTDB_REFSEQ_VERSION=220,
    GTDB_PATIAL_SUFFIXES=["a" + i for i in "abcdefghijk"],
):
    return expand(
        GTDBTK_DATA_PATH / "gtdbtk_r{gtdb_refseq_version}_data.tar.gz.part_{suffix}",
        gtdb_refseq_version=GTDB_REFSEQ_VERSION,
        suffix=GTDB_PATIAL_SUFFIXES,
    )


localrules:
    download_partial_gtdb,
    extract_gtdb,


rule download_partial_gtdb:
    output:
        temp(
            GTDBTK_DATA_PATH
            / "gtdbtk_r{gtdb_refseq_version}_data.tar.gz.part_{suffix}"
        ),
    threads: 1
    params:
        url=lambda wc, output: f"{GTDB_DATA_URL}/split_package/{ Path(output[0]).name}",
    resources:
        time_min=60 * int(config.get("runtime", {"long": 10})["long"]),
    log:
        "logs/download/gtdbtk_r{gtdb_refseq_version}_part_{suffix}.log",
    shell:
        " wget --no-check-certificate {params.url} -O {output} &> {log} "


rule extract_gtdb:
    input:
        all_partial_gtdb_tarbals,
    output:
        touch(os.path.join(GTDBTK_DATA_PATH, "downloaded_success")),
    threads: 1
    resources:
        time_min=60 * int(config.get("runtime", {"long": 10})["long"]),
    log:
        stdout="logs/download/gtdbtk_untar.log",
        stderr="logs/download/gtdbtk_untar.err",
    shell:
        '( cat {input} | tar -xzvf - -C "{GTDBTK_DATA_PATH}" --strip 1 ) 2> {log.stderr} > {log.stdout} '


### end GTDBTk


def md5(fname):
    # https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
    hash_md5 = hashlib.md5()
    if not os.path.exists(fname):
        return None
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


# note: saving OG_fasta.tar.gz in order to not create secondary "success" file
FILES = {
    "adapters.fa": "ae839dc79cfb855a1b750a0d593fe01e",
    "phiX174_virus.fa": "82516880142e8c89b466bc6118696c47",
    "silva_rfam_all_rRNAs.fa": "f102e35d9f48eabeb0efe9058559bc66",
}


def get_eggnog_db_file():
    return ancient(
        expand(
            "{path}/{files}",
            path=EGGNOG_DIR,
            files=["eggnog.db", "eggnog_proteins.dmnd"],
        )
    )


localrules:
    download,
    download_eggNOG_files,
    download_atlas_files,
    download_gunc,


ruleorder: download_eggNOG_files > download_atlas_files


rule download:
    input:
        expand(
            "{dir}/{filename}", dir=DBDIR, filename=["adapters.fa", "phiX174_virus.fa"]
        ),
        get_eggnog_db_file(),
        f"{DBDIR}/CheckM2",
        GTDBTK_DATA_PATH / "downloaded_success",


rule download_eggNOG_files:
    output:
        f"{EGGNOG_DIR}/eggnog.db",
        f"{EGGNOG_DIR}/eggnog_proteins.dmnd",
    threads: 1
    conda:
        "../envs/eggNOG.yaml"
    shell:
        f"download_eggnog_data.py -yf --data_dir {EGGNOG_DIR} "


rule download_atlas_files:
    output:
        f"{DBDIR}/{{filename}}",
    threads: 1
    wildcard_constraints:
        filename="[A-Za-z0-9_.]+",
    run:
        shell(
            "wget -O {output} 'https://zenodo.org/record/{ZENODO_ARCHIVE}/files/{wildcards.filename}' "
        )
        if not FILES[wildcards.filename] == md5(output[0]):
            raise OSError(2, "Invalid checksum", output[0])


rule checkm2_download_db:
    output:
        directory(f"{DBDIR}/CheckM2"),
    conda:
        "../envs/checkm2.yaml"
    threads: 1
    log:
        "logs/download/checkm2.log",
    resources:
        time_min=60 * int(config.get("runtime", {"long": 10})["long"]),
    shell:
        " checkm2 database --download --path {output} "
        " &>> {log}"


rule download_gunc:
    output:
        os.path.join(GUNCDIR, "{gunc_database}"),
    conda:
        "../envs/gunc.yaml"
    threads: 1
    resources:
        time_min=60 * int(config.get("runtime", {"default": 5})["default"]),
        mem_mb=config.get("simplejob_mem", 1) * 1000,
        tmpdir=config.get("tmpdir", "."),  # you can store the file in the main working folder if you want
    log:
        "logs/downloads/gunc_download_{gunc_database}.log",
    shell:
        "gunc download_db {resources.tmpdir} -db {wildcards.gunc_database} &> {log} ;"
        "mv {resources.tmpdir}/gunc_db_{wildcards.gunc_database}*.dmnd {output} 2>> {log}"


rule download_busco:
    output:
        directory(BUSCODIR),
    conda:
        "../envs/busco.yaml"
    threads: 1
    resources:
        time_min=60 * int(config.get("runtime", {"default": 5})["default"]),
        mem_mb=config.get("simplejob_mem", 1) * 1000,
    log:
        "logs/busco_lineages.log",
    shell:
        "busco -q --download_path {output} --download prokaryota &> {log}"


onsuccess:
    print("All databases have downloaded and validated successfully")


onerror:
    print("An error occurred while downloading reference databases.")
    # print(
    #     "ATLAS databases can be manually downloaded from: https://zenodo.org/record/%s"
    #     % ZENODO_ARCHIVE
    # )
    # print(
    #     "eggNOG databases can be manually downloaded from: http://eggnogdb.embl.de/download/emapperdb-%s"
    #     % EGGNOG_VERSION
    # )

