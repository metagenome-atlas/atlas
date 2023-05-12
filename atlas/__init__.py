import snakemake
from . import _version
import os

from .workflow.scripts import utils


TAX_LEVELS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
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
