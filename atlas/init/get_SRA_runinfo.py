#!/usr/bin/env python
import os
import requests
import csv, io
import traceback
import xmltodict

try:
    from tqdm import tqdm
except ImportError:
    tqdm = list


class SRAUtils:
    """
    Prefixes can be found here http://www.ddbj.nig.ac.jp/sub/prefix.html for DDBJ, ENA/EBI (ERA), and NCBI(SRA)

    Some interesting information here http://trace.ddbj.nig.ac.jp/dra/submission_e.html#Organization_of_metadata_objects
    """

    @staticmethod
    def getInfoTableFromSearchTerm(search):
        # return types and databases listed here
        # https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly
        #
        payload = {"save": "efetch", "db": "sra", "rettype": "runinfo", "term": search}
        r = requests.get(
            "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi", params=payload
        )
        if 200 == r.status_code:
            if r.text.isspace():
                raise Exception("Got blank string from " + str(r.url))
            else:
                reader_list = csv.DictReader(io.StringIO(r.text))
                infoRows = []
                for row in reader_list:
                    infoRows.append(row)
                if 0 == len(infoRows):
                    raise Exception(
                        'Found %d entries in SRA for "%s" when expecting at least 1'
                        % (len(infoRows), search)
                    )
                else:
                    return infoRows
                return infoRows
        else:
            raise Exception(
                "Error in downloading from "
                + str(r.url)
                + " got response code "
                + str(r.status_code)
            )

    @staticmethod
    def getDictTableFromSearchTerm(search):
        # return types and databases listed here
        # https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly
        # @todo need to add way to select only what is needed
        payload = {"save": "efetch", "db": "sra", "rettype": "full", "term": search}
        r = requests.get(
            "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi", params=payload
        )
        if 200 == r.status_code:
            if r.text.isspace():
                raise Exception("Got blank string from " + str(r.url))
            else:
                return xmltodict.parse(r.text)
                """
                reader_list = csv.DictReader(io.StringIO())
                infoRows = []
                for row in reader_list:
                    infoRows.append(row)
                if 0 == len(infoRows):
                    raise Exception('Found %d entries in SRA for "%s" when expecting at least 1' % (len(infoRows), search))
                else:
                    return infoRows
                return infoRows
                """
        else:
            raise Exception(
                "Error in downloading from "
                + str(r.url)
                + " got response code "
                + str(r.status_code)
            )

    @staticmethod
    def getRunAccsFromInfoTable(infoTab):
        runInfo = []
        for row in infoTab:
            runInfo.append(str(row.get("Run")))
        return runInfo

    @staticmethod
    def getSraUrlFromRunAccession(accesion):
        if not type(accesion) is str:
            raise Exception(
                "Error in getSraUrlFromRunAccession: accesion should be str, not "
                + str(type(accesion))
            )
        if len(accesion) < 7:
            raise Exception(
                "Error in getSraUrlFromRunAccession: accession should be least 7 character long, not "
                + str(len(accesion))
                + ", for "
                + str(accesion)
            )

        if not accesion.startswith("ERR") and not accesion.startswith("SRR"):
            raise Exception(
                "Error in getSraUrlFromRunAccession: accession should start with either ERR or SRR, not "
                + str(accesion[0:3])
            )

        template = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/{PREFIX}/{PREFIX}{PREFIXNUMS}/{ACCESION}/{ACCESION}.sra"
        return template.format(
            PREFIX=accesion[0:3], PREFIXNUMS=accesion[3:6], ACCESION=accesion
        )

    @staticmethod
    def getInfoFromRunAcc(run):
        if not run.startswith("ERR") and not run.startswith("SRR"):
            raise Exception("run should start with ERR or SRR, not: " + run)
        infoTab = SRAUtils.getInfoTableFromSearchTerm(run)
        runInfo = []
        for row in infoTab:
            if run == str(row.get("Run")):
                runInfo.append(row)
        if 0 == len(runInfo):
            raise Exception(
                'Found %d entries in SRA for "%s" when expecting at least 1'
                % (len(runInfo), run)
            )
        else:
            return runInfo

    @staticmethod
    def getInfoFromSubmissionAcc(submission):
        if not submission.startswith("ERA") and not submission.startswith("SRA"):
            raise Exception(
                "submission should start with ERR or SRR, not: " + submission
            )
        infoTab = SRAUtils.getInfoTableFromSearchTerm(submission)
        runInfo = []
        for row in infoTab:
            if submission == str(row.get("Submission")):
                runInfo.append(row)
        if 0 == len(runInfo):
            raise Exception(
                'Found %d entries in SRA for "%s" when expecting at least 1'
                % (len(runInfo), submission)
            )
        else:
            return runInfo

    @staticmethod
    def getInfoFromSampleAcc(sample):
        if not sample.startswith("ERS") and not sample.startswith("SRS"):
            raise Exception("sample should start with ERS or SRS, not: " + sample)
        infoTab = SRAUtils.getInfoTableFromSearchTerm(sample)
        runInfo = []
        for row in infoTab:
            if sample == str(row.get("Sample")):
                runInfo.append(row)
        if 0 == len(runInfo):
            raise Exception(
                'Found %d entries in SRA for "%s" when expecting at least 1'
                % (len(runInfo), sample)
            )
        else:
            return runInfo

    @staticmethod
    def getInfoFromProjectAcc(project):
        if not project.startswith("ERP") and not project.startswith("SRP"):
            raise Exception("project should start with ERP or SRP, not: " + project)
        infoTab = SRAUtils.getInfoTableFromSearchTerm(project)
        runInfo = []
        for row in infoTab:
            if project == str(row.get("SRAStudy")):
                runInfo.append(row)
        if 0 == len(runInfo):
            raise Exception(
                'Found %d entries in SRA for "%s" when expecting at least 1'
                % (len(runInfo), project)
            )
        else:
            return runInfo

    @staticmethod
    def getInfoFromExperimentAcc(experiment):
        if not experiment.startswith("ERX") and not experiment.startswith("SRX"):
            raise Exception(
                "experiment should start with ERX or SRX, not: " + experiment
            )
        infoTab = SRAUtils.getInfoTableFromSearchTerm(experiment)
        runInfo = []
        for row in infoTab:
            if experiment == str(row.get("Experiment")):
                runInfo.append(row)
        if 0 == len(runInfo):
            raise Exception(
                'Found %d entries in SRA for "%s" when expecting at least 1'
                % (len(runInfo), experiment)
            )
        else:
            return runInfo

    @staticmethod
    def getInfoFromBioProjectAcc(bioProject):
        # or
        if (
            not bioProject.startswith("PRJNA")
            and not bioProject.startswith("PRJEA")
            and not bioProject.startswith("PRJEB")
            and not bioProject.startswith("PRJDA")
        ):
            raise Exception(
                "bioProject should start with PRJNA, PRJEA, PRJEB, and PRJDA, not: "
                + bioProject
            )
        infoTab = SRAUtils.getInfoTableFromSearchTerm(bioProject)
        runInfo = []
        for row in infoTab:
            if bioProject == str(row.get("BioProject")):
                runInfo.append(row)
        if 0 == len(runInfo):
            raise Exception(
                'Found %d entries in SRA for "%s" when expecting at least 1'
                % (len(runInfo), bioProject)
            )
        else:
            return runInfo

    @staticmethod
    def getInfoFromBioSampleAcc(bioSample):
        """
        SAMD     DDBJ
        SAME     ENA/EBI
        SAMN     NCBI
        """
        if (
            not bioSample.startswith("SAME")
            and not bioSample.startswith("SAMD")
            and not bioSample.startswith("SAMN")
        ):
            raise Exception(
                "bioProject should start with SAME, SAMD, and SAMN, not: " + bioSample
            )
        infoTab = SRAUtils.getInfoTableFromSearchTerm(bioSample)
        runInfo = []
        for row in infoTab:
            if bioSample == str(row.get("BioSample")):
                runInfo.append(row)
        if 0 == len(runInfo):
            raise Exception(
                'Found %d entries in SRA for "%s" when expecting at least 1'
                % (len(runInfo), bioSample)
            )
        else:
            return runInfo

    @staticmethod
    def getInfoFromSRAIdentifier(identifier):
        if identifier[1:3] == "RX":
            return SRAUtils.getInfoFromExperimentAcc(identifier)
        elif identifier[1:3] == "RP":
            return SRAUtils.getInfoFromProjectAcc(identifier)
        elif identifier[1:3] == "RS":
            return SRAUtils.getInfoFromSampleAcc(identifier)
        elif identifier[1:3] == "RR":
            return SRAUtils.getInfoFromRunAcc(identifier)
        elif identifier[1:3] == "RA":
            return SRAUtils.getInfoFromSubmissionAcc(identifier)
        elif identifier.startswith("SAM"):
            # SAME, SAMD, and SAMN
            return SRAUtils.getInfoFromBioSampleAcc(identifier)
        elif identifier.startswith("PRJ"):
            # DDBJ archvie bioproject prefix PRJNA SAMEA2796165
            return SRAUtils.getInfoFromBioProjectAcc(identifier)
        else:
            raise Exception(
                "Error, unrecognized prefix for sra Identifier " + str(identifier)
            )


# main function


def get_runtable_from_ids(identifiers, output_file="SRA_runtable.tsv", overwrite=False):

    if type(identifiers) == str:
        identifiers = [identifiers]

    with open(output_file, "w") as outInfoFile:
        identifierCount = 0

        # don't show progress bar if only one elelment
        if len(identifiers) > 1:
            identifier_with_progressbar = tqdm(identifiers)
        else:
            identifier_with_progressbar = identifiers

        for identifier in identifier_with_progressbar:
            try:
                tab = SRAUtils.getInfoFromSRAIdentifier(identifier)
                # write header
                if 0 == identifierCount:
                    outInfoFile.write("\t".join(tab[0].keys()) + "\n")
                for row in tab:
                    outInfoFile.write("\t".join(row.values()) + "\n")
                identifierCount = identifierCount + 1

            except Exception as err:
                raise Exception(
                    "Failed  to get info for " + str(identifier) + ", mess: " + str(err)
                ) from err


def parse_arguments_from_terminal():

    ## Comand line interface
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "identifiers",
        type=str,
        help="A list of space separated SRA identifiers e.g. SRP046206 SRX188939 SRS807544 SRR1759594 PRJNA63661",
        nargs="+",
    )
    parser.add_argument(
        "-o",
        "--outprefix",
        type=str,
        help="Output prefix for the sra table",
        default="SRA",
    )
    parser.add_argument(
        "-f",
        "--overwrite",
        action="store_true",
        help="Overwrite table if it already exists",
    )

    # parse arguments from terminal
    args = parser.parse_args()

    # if output prefix is only a folder add 'SRA
    outprefix = args.outprefix
    if outprefix.endswith(os.path.sep):
        outprefix += "SRA"

    # Define output_file
    output_file = outprefix + "_runtable.tsv"

    if os.path.exists(output_file) and not args.overwrite:
        raise Exception(
            "File " + output_file + " already exists, use --overwrite to over write it"
        )

    return args.identifiers, output_file


if __name__ == "__main__":

    get_runtable_from_ids(*parse_arguments_from_terminal())
