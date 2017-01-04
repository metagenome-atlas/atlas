import hashlib
import os


USERNAME = "observers"
PASSWORD = "noVa9lib"


def md5(fname):
    # https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
    hash_md5 = hashlib.md5()
    if not os.path.exists(fname):
        return None
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def touch(fname, times=None):
    # https://stackoverflow.com/questions/1158076/implement-touch-using-python
    with open(fname, 'a'):
        os.utime(fname, times)


rule transfer_files:
    output:
        "%s/{filename}" % OUTDIR
    run:
        shell("curl 'ftp://{USERNAME}:{PASSWORD}@ftp.pnl.gov/outgoing/atlas/{wildcards.filename}' -s > {output}")
        shell("curl 'ftp://{USERNAME}:{PASSWORD}@ftp.pnl.gov/outgoing/atlas/{wildcards.filename}.md5' -s > {output}.md5")
        with open(output[0] + ".md5") as fh:
            for line in fh:
                known_md5 = line.strip().split("  ")[0]
                break
        if not known_md5 == md5(output[0]):
            raise OSError(2, "Invalid checksum", output[0])


onsuccess:
    print(("All databases have downloaded and validated successfully.\nWhen generating your "
           "configuration file, use '--database-dir %s'") % OUTDIR)
