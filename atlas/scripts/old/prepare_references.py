import click
import csv
import gzip
import logging
import re
import sys
from collections import defaultdict
from itertools import groupby


logging.basicConfig(
    level=logging.INFO, datefmt="%Y-%m-%d %H:%M", format="[%(asctime)s] %(message)s"
)


def print_fasta_record(name, seq, out_handle=sys.stdout, wrap=100):
    """Print a fasta record accounting for line wraps in the sequence.

    Args:
        name (str): name or header for fasta entry
        seq (str): sequence
        out_handle (Optional): open file handle in which to write or stdout
        wrap (Optional[int]) : line width of fasta sequence; None is supported for no wrapping

    """
    print(">", name, sep="", file=out_handle)
    if wrap:
        for i in range(0, len(seq), wrap):
            print(seq[i : i + wrap], file=out_handle)
    else:
        print(seq, file=out_handle)


def format_fasta_record(name, seq, wrap=100):
    """Fasta __str__ method.

    Convert fasta name and sequence into wrapped fasta format.

    Args:
        name (str): name of the record
        seq (str): sequence of the record
        wrap (int): length of sequence per line

    Yields:
        tuple: name, sequence

    >>> format_fasta_record("seq1", "ACTG")
    ">seq1\nACTG"
    """
    record = ">" + name + "\n"
    if wrap:
        for i in range(0, len(seq), wrap):
            record += seq[i : i + wrap] + "\n"
    else:
        record += seq + "\n"
    return record.strip()


def read_fasta(fh):
    """Fasta iterator.

    Accepts file handle of .fasta and yields name and sequence.

    Args:
        fh (file): Open file handle of .fasta file

    Yields:
        tuple: name, sequence

    Example:
        >>> import os
        >>> from itertools import groupby
        >>> f = open("test.fasta", 'w')
        >>> f.write("@seq1\nACTG")
        >>> f.close()
        >>> f = open("test.fastq")
        >>> for name, seq in read_fastq(f):
                assert name == "seq1"
                assert seq == "ACTG"
        >>> f.close()
        >>> os.remove("test.fasta")

    """
    for header, group in groupby(fh, lambda line: line[0] == ">"):
        if header:
            line = next(group)
            name = line[1:].strip()
        else:
            seq = "".join(line.strip() for line in group)
            yield name, seq


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option("0.0.0")
@click.pass_context
def cli(obj):
    "Methods to prepare databases for parsing BLAST tabular output."


@cli.command("prepare-metacyc", short_help="prepares metacyc mapping files")
@click.argument("fasta", type=click.Path(exists=True))
@click.argument("seqids", type=click.Path(exists=True))
@click.argument("reactions", type=click.Path(exists=True))
@click.argument("pathwaylinks", type=click.Path(exists=True))
@click.argument("unirefclusters", type=click.Path(exists=True))
@click.argument("outmap", type=click.Path())
@click.argument("outfasta", type=click.Path())
def prepare_metacyc_reference(
    fasta, seqids, reactions, pathwaylinks, unirefclusters, outmap, outfasta
):
    """
    # via uniprot
    fasta uniref100.fasta.gz
    # via metacyc
    seqids 20.5/data/uniprot-seq-ids.dat
    reactions 20.5/data/reactions.data
    pathwaylinks 20.5/data/pathway-links.dat

    uniref mappings where cluster n>1 (uniref search): count:[2 TO *] AND identity:1.0

    needs cluster id and cluster members
    """

    logging.info("Parsing sequence IDs from %s" % seqids)
    seqids_obj = ""
    with open(seqids) as fh:
        for line in fh:
            if line.startswith(";;"):
                continue
            if not line.strip():
                continue
            seqids_obj += (
                " " + line.strip() if line.strip().startswith('"') else line.strip()
            )

    # strip the first and last parentheses
    seqids_str = seqids_obj[1:-1]

    # fix the parentheses within entries
    for reg in re.findall(r"\|(.*?)\|", seqids_obj[1:-1]):
        # remove parenthesis between pipes
        seqids_str = seqids_str.replace(reg, reg.replace("(", "").replace(")", ""))
    seqids_str = seqids_str.replace("|", "")

    logging.info("Parsing pathway links from %s" % pathwaylinks)
    # read in pathway links
    pathway_links = {}
    with open(pathwaylinks) as fh:
        for line in fh:
            # deal with header lines
            if line.startswith("#"):
                continue
            toks = line.strip().split("\t")
            pathway_str = toks[1]
            # formaldehyde oxidation V (H<sub>4</sub>MPT pathway)
            # superpathway of anthocyanin biosynthesis (from delphinidin 3-<i>O</i>-glucoside)
            for http_format in re.findall(r"(\<.*?\>)", toks[1]):
                pathway_str = pathway_str.replace(http_format, "")
            # &beta; --> beta
            pathway_str = pathway_str.replace("&beta;", "beta")
            # ignores synonyms in toks[2:]
            pathway_links[toks[0]] = pathway_str

    logging.info("Parsing reactions of %s" % reactions)
    # read in reactions for multiple ECs and pathways
    reaction_links = {}
    with open(reactions, "r", encoding="ISO-8859-1") as fh:
        pathways = []
        ec_numbers = []
        for line in fh:
            # EC-NUMBER - |EC-1.14.19.ar| --> EC-NUMBER - EC-1.14.19.ar
            line = line.replace("|", "")
            if "UNIQUE-ID" in line:
                name = line.strip().partition(" - ")[-1]
            if "EC-NUMBER" in line:
                ec_numbers.append(line.strip().partition(" - ")[-1].replace("EC-", ""))
            if "IN-PATHWAY" in line:
                pathways.append(line.strip().partition(" - ")[-1])
            if line.startswith("//"):
                # store current record
                if not name:
                    sys.exit("{ec_numbers}; {pathways}".format(**locals()))
                pathway_names = []
                for pathway in pathways:
                    try:
                        pathway_names.append(pathway_links[pathway])
                    except KeyError:
                        # some pathways do not have a description
                        pass
                reaction_links[name] = {
                    "ec": ec_numbers,
                    "pathways": pathways,
                    "pathway_names": pathway_names,
                }
                name = ""
                ec_numbers = []
                pathways = []

    logging.info("Annotating UniProt IDs across maps")
    uniprot_to_reactions = defaultdict(set)
    uniprot_to_ecs = defaultdict(set)
    uniprot_to_pathways = defaultdict(set)

    for metacyc_entry in re.findall(r"\((.*?)\)", seqids_str):
        # space delimited string
        toks = [t.replace('"', "") for t in metacyc_entry.split()]

        for uid in toks[2:]:
            uniprot_to_reactions[uid].add(toks[0])
            try:
                for ec_id in reaction_links[toks[0]]["ec"]:
                    uniprot_to_ecs[uid].add(ec_id)
                for pathway_id in reaction_links[toks[0]]["pathways"]:
                    uniprot_to_pathways[uid].add(pathway_id)
            except KeyError:
                # reaction that has uniprot seqs, but lacks further documentation in reactions.dat
                # using inline EC number
                uniprot_to_ecs[uid].add(toks[1])
                uniprot_to_pathways[uid].add("")

    # despite cluster translation we still expect a few to be missing
    clusters = {}
    with gzip.open(unirefclusters, "rt") as fh:
        for line in fh:
            found = False
            toks = line.strip().split("\t")
            uniref_id = toks[0].partition("_")[2]
            uids = toks[2].split("; ")
            for uid in uids:
                if uid in uniprot_to_reactions:
                    found = True
            if found:
                clusters[uniref_id] = uids

    for uniref_id, uids in clusters.items():
        combined_reactions = set()
        combined_ecs = set()
        combined_pathways = set()
        for uid in uids:
            combined_reactions.update(uniprot_to_reactions[uid])
            combined_ecs.update(uniprot_to_ecs[uid])
            combined_pathways.update(uniprot_to_pathways[uid])
            uniprot_to_reactions.pop(uid)
            uniprot_to_ecs.pop(uid)
            uniprot_to_pathways.pop(uid)

        uniprot_to_reactions[uniref_id] = combined_reactions
        uniprot_to_ecs[uniref_id] = combined_ecs
        uniprot_to_pathways[uniref_id] = combined_pathways

    # print metacyc map file
    with open(outmap, "w") as fo:
        for uid in uniprot_to_reactions.keys():
            pathway_descriptions = set()
            for pathway_id in uniprot_to_pathways[uid]:
                try:
                    pathway_descriptions.add(pathway_links[pathway_id])
                except KeyError:
                    # pathway is defined in reactions.dat, but not in present in pathway-links.dat
                    pass

            print(
                uid,
                "|".join(uniprot_to_reactions[uid]),
                "|".join(uniprot_to_ecs[uid]),
                "|".join(uniprot_to_pathways[uid]),
                "|".join(pathway_descriptions),
                sep="\t",
                file=fo,
            )

    logging.info("%d unique UniProt entries in the map" % len(uniprot_to_reactions))

    logging.info("Grabbing reference sequences from %s" % fasta)
    # >UniRef100_Q6GZX4 Putative transcription factor 001R n=2 Tax=Frog virus 3 TaxID=10493 RepID=001R_FRG3G
    counter = 0
    with gzip.open(fasta, "rt") as fh, open(outfasta, "w") as fo:
        for i, (name, seq) in enumerate(read_fasta(fh)):
            uniprot_id = name.partition(" ")[0].partition("_")[-1]
            if uniprot_id in uniprot_to_reactions:
                counter += 1
                print_fasta_record(uniprot_id, seq, fo)
        logging.info(
            "%d unique sequences in the reference fasta out of %d input" % (counter, i)
        )


@cli.command("prepare-refseq", short_help="prepares refseq mapping files")
@click.argument("fasta", type=click.Path(exists=True))
@click.argument("namesdmp", type=click.Path(exists=True))
@click.argument("nodesdmp", type=click.Path(exists=True))
@click.argument("namemap", type=click.File("w"))
@click.argument("tree", type=click.File("w"))
def prepare_refseq_reference(fasta, namesdmp, nodesdmp, namemap, tree):
    """Takes the RefSeq FASTA which will look like:

        \b
        >WP_016124266.1 two-component sensor histidine kinase [Bacillus cereus]
        MNFNKRLIIQFIMQHVFVLVTLLIAVVAAFTYLIFLLTSTLYEPNIPDSDSFTISRYISSEDGHISLQSEVQDLIKEKND

    And from NCBI's taxonomy dump, names.dmp, which starts like:

        \b
        1 | all      |                         | synonym         |
        1 | root     |                         | scientific name |
        2 | Bacteria | Bacteria <prokaryotes>  | scientific name |

    And nodes.dmp:

        \b
        1 | 1      | no rank      |...
        2 | 131567 | superkingdom |...
        6 | 335928 | genus        |...

    To create a TSV map (NAMEMAP) of FASTA entry ID to function, taxonomy name, taxonomy ID, and parent
    taxonomy ID:

        \b
        WP_016124266.1<tab>two-component sensor histidine kinase<tab>Bacillus cereus

    As well as TREE:

        \b
        Bacillus cereus<tab>1396<tab>86661<tab>species

    """

    # most names, for translating reference fasta from tax name to tax id
    name_to_tax = {}
    tax_to_name = {}
    # for best names, for translating nodes.dmp tax id to tax name
    tax_to_scientific_name = {}

    logging.info("Reading in %s" % namesdmp)
    with open(namesdmp) as dmp:
        for tax_id, group in groupby(
            dmp, key=lambda x: [i.strip() for i in x.strip().split("|")][0]
        ):
            scientific_name = ""
            synonym = ""
            for line in group:
                toks = [
                    x.strip().replace("'", "").replace('"', "")
                    for x in line.strip().split("|")
                ]
                if toks[-2] == "scientific name":
                    tax_to_scientific_name[toks[0]] = toks[1]
                elif toks[-2] == "misspelling":
                    name_to_tax[toks[1]] = toks[0]
                    continue
                elif toks[-2] == "synonym":
                    synonym = toks[1]
                else:
                    if synonym:
                        tax_to_name[toks[0]] = synonym
                    else:
                        tax_to_name[toks[0]] = toks[1]
                name_to_tax[toks[1]] = toks[0]

    logging.info("Iterating over %s" % fasta)
    with gzip.open(fasta, "rt") as fa:
        for name, seq in read_fasta(fa):
            name_parts = name.partition(" ")
            # biotin--[acetyl-CoA-carboxylase] synthetase [Aeromonas allosaccharophila]
            function = name_parts[2].rpartition("[")[0].strip()
            taxonomy_name = (
                name_parts[2]
                .rpartition("[")[2]
                .strip("]")
                .replace("'", "")
                .replace('"', "")
            )
            # phosphoribosyltransferase [[Haemophilus] parasuis]
            if "[[" in name_parts[2]:
                taxonomy_name = "[%s" % taxonomy_name
            # https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?lvl=0&id=1737424
            if taxonomy_name == "Blautia sp. GD8":
                taxonomy_name = "Blautia sp. GD9"
            if taxonomy_name == "unclassified Actinobacteria (miscellaneous)":
                taxonomy_name = "unclassified Actinobacteria (class) (miscellaneous)"
            try:
                taxonomy_id = name_to_tax[taxonomy_name]
            except KeyError:
                taxonomy_name = taxonomy_name.replace("]", "").replace("[", "")
                taxonomy_id = name_to_tax[taxonomy_name]

            print(name_parts[0], function, taxonomy_id, sep="\t", file=namemap)

    logging.info("Reading in %s" % nodesdmp)
    with open(nodesdmp) as dmp:
        for line in dmp:
            toks = [x.strip() for x in line.strip().split("|")]
            try:
                taxonomy_name = tax_to_scientific_name[toks[0]]
            except KeyError:
                taxonomy_name = tax_to_name[toks[0]]
            print(toks[0], taxonomy_name, toks[1], toks[2], sep="\t", file=tree)

    logging.info("Process complete")


@cli.command("prepare-eggnog", short_help="prepares eggnog mapping and reference files")
@click.argument("fasta", type=click.File("r"))
@click.argument("namemap", type=click.Path(exists=True))
@click.argument("outputfasta", type=click.File("w"))
@click.argument("outputmap", type=click.Path())
def prepare_eggnog_reference(fasta, namemap, outputfasta, outputmap):
    names = set()
    # IDs for which we have uniprot AC
    with open(namemap, "r", encoding="ISO-8859-1") as nm:
        reader = csv.reader(nm, delimiter="\t")
        # skip the header
        next(reader)
        for toks in reader:
            names.add(toks[1])
    logging.info("Total names in %s: %d" % (namemap, len(names)))
    # IDs for which we have sequences
    fasta_names = set()
    fasta_counter = 0
    for name, seq in read_fasta(fasta):
        fasta_counter += 1
        name = name.partition(".")[-1]
        # we were unable to translate this to meaningful metadata
        if not name in names:
            continue
        # save the entry
        fasta_names.add(name)
        # print this nonredundant fasta
        print_fasta_record(name, seq, outputfasta)
    logging.info(
        "Total unique matching fasta entries: %d (from %d entries)"
        % (len(fasta_names), fasta_counter)
    )
    # the union of IDs
    with open(namemap, "r", encoding="ISO-8859-1") as nm, open(outputmap, "w") as ofh:
        reader = csv.reader(nm, delimiter="\t")
        logging.info("Finding union of name map and fasta")
        for toks in reader:
            if toks[1] in fasta_names:
                print(*toks, sep="\t", file=ofh)


@cli.command("prepare-cazy", short_help="prepares cazy (dbcan) reference")
@click.argument("faminfo", click.Path(exists=True))
@click.argument("cazydb_fasta", click.Path(exists=True))
@click.argument("out_map")
@click.argument("out_fasta")
def prepare_cazy(faminfo, cazydb_fasta, out_map, out_fasta):
    """

    faminfo:

        \b
        http://csbl.bmb.uga.edu/dbCAN/download/FamInfo.txt

    cazydb_fasta is:

        \b
        http://csbl.bmb.uga.edu/dbCAN/download/CAZyDB.07152016.fa

    """
    click.echo("parsing %s" % faminfo)
    fam_info = {}
    with open(faminfo) as fh:
        next(fh)
        for line in fh:
            toks = line.strip().split("\t")
            if "Unclassified" in toks[0]:
                toks[0] = toks[0].replace("Unclassified-", "")
            if toks[0] in fam_info:
                print(line)
                sys.exit("Non-unique family present in Family Info")
            # unsure how useful cazy_activities is here as it includes everything across all families
            fam_info[toks[0]] = toks[2]

    click.echo("parsing %s" % cazydb_fasta)
    cazy_fasta_map = {}
    with open(cazydb_fasta) as fh:
        for name, seq in read_fasta(fh):
            name_parts = name.split("|", 2)
            cazy_fasta_map[name_parts[0]] = {
                "seq": seq,
                "ecs": "" if len(name_parts) == 2 else name_parts[2],
                "cazy_family": name_parts[1],
            }

    with open(out_map, "w") as omap, open(out_fasta, "w") as ofa:
        for name, meta in cazy_fasta_map.items():
            # removes masked sequences
            if meta["seq"].islower():
                continue
            print(format_fasta_record(name, meta["seq"]), file=ofa)
            # cazy_gene, cazy_family, cazy_class, cazy_ecs

            if not meta["cazy_family"] in fam_info:
                cazy_class = re.findall(r"([A-Z]+)", meta["cazy_family"])[0]
                if not cazy_class:
                    cazy_class = "NA"
            else:
                cazy_class = fam_info[meta["cazy_family"]]

            print(
                name, meta["cazy_family"], cazy_class, meta["ecs"], sep="\t", file=omap
            )


@cli.command("prepare-enzyme", short_help="prepares EC reference from ENZYME")
@click.argument("enzyme_dat", click.Path(exists=True))
@click.argument("uniparc_map", click.Path(exists=True))
@click.argument("uniparc_fasta", click.Path(exists=True))
@click.argument("out_map")
@click.argument("out_fasta")
def prepare_enzyme(enzyme_dat, uniparc_map, uniparc_fasta, out_map, out_fasta):
    """

    enzyme_dat:

        \b
        ftp://ftp.expasy.org/databases/enzyme/enzyme.dat

    uniparc_map is the compressed version of was 'all' from:

        \b
        http://www.uniprot.org/uniparc/

    uniparc_fasta is:

        \b
        ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniparc/uniparc_active.fasta.gz

    """
    enzyme = set()
    uniprot_to_uniparc = {}

    click.echo("parsing %s" % uniparc_map)
    with gzip.open(uniparc_map, "rt") as fh:
        # Entry    Organisms    UniProtKB    First seen    Last seen    Length
        next(fh)
        for line in fh:
            toks = line.strip().split("\t")
            for uniprot_entry in toks[2].split(";"):
                uniprot_entry = uniprot_entry.strip()

                if not uniprot_entry or "obsolete" in uniprot_entry:
                    continue
                if uniprot_entry in uniprot_to_uniparc:
                    print(
                        "Uniprot Entry", uniprot_entry, "is not unique", "\nLINE:", line
                    )
                    sys.exit(1)

                uniprot_to_uniparc[uniprot_entry] = toks[0]

    click.echo("parsing %s" % enzyme_dat)
    uniparc_mappings = defaultdict(list)
    with open(enzyme_dat) as fh:
        for key, group in groupby(fh, key=lambda x: x.startswith("//")):
            if not key:
                uniprot_entries = []
                recommended_name = "NA"
                ec_id = "NA"

                for line in group:
                    toks = line.strip().partition("   ")

                    if line.startswith("ID"):
                        ec_id = toks[2]

                    if line.startswith("DE"):
                        if recommended_name == "NA":
                            recommended_name = toks[2].strip(".")
                        # multi-line name
                        else:
                            if recommended_name.endswith("-"):
                                recommended_name += toks[2].strip(".")
                            else:
                                recommended_name += " " + toks[2].strip(".")

                    if line.startswith("DR"):
                        for entry_and_name in toks[2].split(";"):
                            uniprot_entry = entry_and_name.split(",")[0].strip()

                            if uniprot_entry:
                                uniprot_entries.append(uniprot_entry)

                if ec_id and uniprot_entries:

                    for uniprot_entry in uniprot_entries:
                        uniparc_id = uniprot_to_uniparc[uniprot_entry]

                        uniparc_mappings[uniparc_id].append(
                            [uniprot_entry, ec_id, recommended_name]
                        )
                        enzyme.add(uniparc_id)

    with open(out_map, "w") as ofh:
        for uniparc_id, meta in uniparc_mappings.items():
            uniprots = []
            ecs = []
            names = []
            for name_list in meta:
                if name_list[0] not in uniprots:
                    uniprots.append(name_list[0])

                if name_list[1] not in ecs:
                    ecs.append(name_list[1])
                    # ECs are coupled with a name
                    names.append(name_list[2])

            print(
                uniparc_id,
                "|".join(uniprots),
                "|".join(ecs),
                "|".join(names),
                sep="\t",
                file=ofh,
            )

    click.echo("parsing %s" % uniparc_fasta)
    with gzip.open(uniparc_fasta, "rt") as fh, open(out_fasta, "w") as ofh:
        for name, seq in read_fasta(fh):
            name = name.partition(" ")[0]
            if name in enzyme:
                print(format_fasta_record(name, seq), file=ofh)


@cli.command("prepare-cog", short_help="prepares COG mapping files")
@click.argument("fasta", type=click.Path(exists=True))
@click.argument("namemap", type=click.Path(exists=True))
@click.argument("funcdef", type=click.Path(exists=True))
@click.argument("namedef", type=click.Path(exists=True))
@click.argument("out_fasta", type=click.File("w"))
@click.argument("out_map", type=click.File("w"))
def prepare_refseq_reference(fasta, namemap, funcdef, namedef, out_fasta, out_map):
    """
    COG data downloaded from: ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data

    fasta: prot2003-2014.fa.gz

    namemap: cog2003-2014.csv

        \b
        158333741,Acaryochloris_marina_MBIC11017_uid58167,158333741,432,1,432,COG0001,0,
        158339504,Acaryochloris_marina_MBIC11017_uid58167,158339504,491,1,491,COG0001,0,
        379012832,Acetobacterium_woodii_DSM_1030_uid88073,379012832,430,1,430,COG0001,0,

    funcdef: fun2003-2014.tab

        \b
        # Code    Name
        J         Translation, ribosomal structure and biogenesis
        A         RNA processing and modification

    namedef: cognames2003-2014.tab

        \b
        # COG    func    name
        COG0001  H       Glutamate-1-semialdehyde aminotransferase
        COG0002  E       N-acetyl-gamma-glutamylphosphate reductase

    """

    class_to_description = {}
    with open(funcdef) as fh:
        next(fh)
        for line in fh:
            toks = line.strip().split("\t")
            assert len(toks) == 2
            assert toks[0] not in class_to_description
            class_to_description[toks[0]] = toks[1]

    cog_to_name = {}
    with open(namedef, encoding="ISO-8859-1") as fh:
        next(fh)
        for line in fh:
            toks = line.strip().split("\t")
            assert len(toks) == 3
            assert toks[0] not in cog_to_name
            cog_to_name[toks[0]] = {"function": toks[1], "name": toks[2]}

    cog_map = defaultdict(dict)
    with open(namemap) as fh:
        for line in fh:
            # 158333741,Acaryochloris_marina_MBIC11017_uid58167,158333741,432,1,432,COG0001,0,
            toks = line.strip().split(",")

            protein_id = toks[0]
            domain_start = int(toks[4])
            domain_stop = int(toks[5])
            cog_id = toks[6]
            try:
                cog_functional_class = cog_to_name[cog_id]["function"]
                cog_annotation = cog_to_name[cog_id]["name"]
            except KeyError:
                # we're only keeping current definitions
                # ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/static/lists/listCOGs.html
                continue

            cog_functional_class_description = "; ".join(
                [class_to_description[i] for i in cog_functional_class]
            )

            # update existing COG start and stop
            if protein_id in cog_map and cog_id in cog_map[protein_id]:

                if domain_start < cog_map[protein_id][cog_id]["domain_start"]:
                    cog_map[protein_id][cog_id]["domain_start"] = domain_start
                if domain_stop > cog_map[protein_id][cog_id]["domain_stop"]:
                    cog_map[protein_id][cog_id]["domain_stop"] = domain_stop
            else:
                cog_map[protein_id][cog_id] = {
                    "cog_functional_class": cog_functional_class,
                    "cog_annotation": cog_annotation,
                    "domain_start": domain_start,
                    "domain_stop": domain_stop,
                    "cog_functional_class_description": cog_functional_class_description,
                }

    with gzip.open(fasta, mode="rt") as fh:
        for name, seq in read_fasta(fh):
            protein_id = name.split("|")[1]
            if protein_id in cog_map:
                for i, (cog_id, cog_data) in enumerate(cog_map[protein_id].items()):
                    cog_seq = seq[
                        cog_data["domain_start"] - 1 : cog_data["domain_stop"]
                    ]
                    # print the fasta entry for a cog sequence
                    print_fasta_record("%s_%d" % (protein_id, i), cog_seq, out_fasta)
                    # print the metadata for a cog; unique primary key is first column
                    print(
                        "%s_%d" % (protein_id, i),
                        protein_id,
                        cog_id,
                        cog_data["cog_functional_class"],
                        cog_data["cog_annotation"],
                        cog_data["cog_functional_class_description"],
                        sep="\t",
                        file=out_map,
                    )


if __name__ == "__main__":
    cli()
