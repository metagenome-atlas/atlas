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


rule transfer_files:
    output:
        "%s/{filename}" % DBDIR
    run:
        eggnog_files = {"OG_fasta.tar.gz": "OG_fasta.tar.gz",
            "eggnog.db": "eggnog.db.gz",
            "og2level.tsv": "og2level.tsv.gz",
            "eggnog_proteins.dmnd": "eggnog_proteins.dmnd.gz",
        }

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


rule extract_checkm_data:
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


onsuccess:
    print(("All databases have downloaded and validated successfully.\nWhen generating your "
           "configuration file, use '--database-dir %s'") % config["db_dir"])


onerror:
    print("An error occurred while downloading reference databases.")
    print("ATLAS databases can be manually downloaded from: https://zenodo.org/record/%s" % ZENODO_ARCHIVE)
    print("eggNOG databases can be manually downloaded from: http://eggnogdb.embl.de/download/emapperdb-%s" % EGGNOG_VERSION)
