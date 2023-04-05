import snakemake
from . import _version
import sys
from pathlib import Path



def import_utils():

    scripts_folder = Path(__file__).parent/ "workflow"/"scripts"

    sys.path.append(str(scripts_folder))
    # now we can import the utils
    import utils

import_utils()


TAX_LEVELS = ["superkingdom", "phylum", "class",
              "order", "family", "genus", "species"]
BLAST6 = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


__version__ = _version.get_versions()["version"]
