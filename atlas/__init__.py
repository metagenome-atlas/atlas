__version__ = "2.8.0"

import os,sys
sys.path.append(os.path.join(os.path.dirname(__file__),"workflow","scripts"))

import utils

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
