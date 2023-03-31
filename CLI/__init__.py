#import utils
import snakemake
from . import _version
import os
import sys

sys.path.append(os.path.join(os.path.dirname(
    __file__), "../workflow", "scripts"))


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
