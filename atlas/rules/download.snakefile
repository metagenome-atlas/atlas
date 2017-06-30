import hashlib
import os


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
        "%s/{filename}" % config["db_dir"]
    run:
        shell("curl 'https://zenodo.org/record/804435/files/{wildcards.filename}' -s > {output}")
        if not FILES[wildcards.filename] == md5(output[0]):
            raise OSError(2, "Invalid checksum", output[0])


onsuccess:
    print(("All databases have downloaded and validated successfully.\nWhen generating your "
           "configuration file, use '--database-dir %s'") % config["db_dir"])


onerror:
    print("An error occurred while downloading reference databases.")
    print("Databases can be manually downloaded from: https://zenodo.org/record/804435")
