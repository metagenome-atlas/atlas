import os
import time
import requests
import xml.etree.ElementTree as ET
import logging
import collections
import json
import itertools


try:
    from tqdm import tqdm
except ImportError:
    tqdm = list

import pandas as pd


# Define these constants so that they can be referred to in other classes
# without index errors.
STUDY_ACCESSION_KEY = 'study_accession'
RUN_ACCESSION_KEY = 'run'
BASES_KEY = 'bases'
SAMPLE_NAME_KEY = 'sample_name'
NCBI_API_KEY_ENV = 'NCBI_API_KEY'
BIOPROJECT_ACCESSION_KEY = 'bioproject'


#from bird_tool_utils import iterable_chunks
def iterable_chunks(iterable, n):
    '''Given an iterable, return it in chunks of size n. In the last chunk, the
    remaining space is replaced by None entries.
    '''
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=None)


class SraMetadata:
    def add_api_key(self, other_params):
        if NCBI_API_KEY_ENV in os.environ:
            other_params['api_key'] = os.environ[NCBI_API_KEY_ENV]
        return other_params

    def _retry_request(self, description, func):
        '''Retry a reqests.post or requests.get 3 times, returning the request
        when OK, otherwise raising an Exception'''

        num_retries = 3
        sleep_time = 60
        def retrying(i, num_retries=3):
            if i < num_retries-1:
                logging.warning("Retrying request (retry {} of {})".format(i+1, num_retries-1))
        
        for i in range(num_retries):
            try:
                this_res = func()
                if not this_res.ok:
                    logging.warning("Request not OK when {}: {}: {}".format(description, this_res, this_res.text))
                    logging.warning("Sleeping for {} seconds before retrying".format(sleep_time))
                    time.sleep(60)
                    retrying(i)
                else:
                    return this_res
            except Exception as e:
                logging.warning("Exception raised when {}: {}".format(description, e))
                logging.warning("Sleeping for {} seconds before retrying".format(sleep_time))
                time.sleep(60)
                retrying(i)
        raise Exception("Failed to {} after {} attempts".format(description, num_retries))


    def fetch_runs_from_bioprojects(self, bioproject_accessions):
        retmax = 10000
        query_string = " OR ".join(["{}[BioProject]".format(bioproject_accession) for bioproject_accession in bioproject_accessions])
        logging.debug("Querying with string: {}".format(query_string))
        res = requests.get(
            url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            params=self.add_api_key({
                "db": "sra",
                "term": query_string,
                "tool": "kingfisher",
                "email": "kingfisher@github.com",
                "retmax": retmax,
                "usehistory": "y",
                }),
            )
        if not res.ok:
            raise Exception("HTTP Failure when requesting search from bioproject: {}: {}".format(res, res.text))
        root = ET.fromstring(res.text)
        sra_ids = list([c.text for c in root.find('IdList')])
        if len(sra_ids) == retmax:
            logging.warning("Unexpectedly found the maximum number of results for this query, possibly some results will be missing")
        webenv = root.find('WebEnv').text

        # Now convert the IDs into runs
        metadata = self.efetch_metadata_from_ids(webenv, None, len(sra_ids))
        return metadata[RUN_ACCESSION_KEY].to_list()

    def efetch_metadata_from_ids(self, webenv, accessions, num_ids):
        data_frames = []

        retmax = num_ids+10
        logging.debug("Running efetch ..")
        res = self._retry_request(
            'efetch_from_ids',
            lambda: requests.get(
                url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                params=self.add_api_key({
                    "db": "sra",
                    "tool": "kingfisher",
                    "email": "kingfisher@github.com",
                    "webenv": webenv,
                    "query_key": 1
                    }),
                ))
        if not res.ok:
            raise Exception("HTTP Failure when requesting efetch from IDs: {}: {}".format(res, res.text))

        root = ET.fromstring(res.text)

        def try_get(func):
            try:
                return func()
            except AttributeError:
                return ''
            except KeyError:
                return None

        if root.find("ERROR") is not None:
            logging.error("Error when fetching metadata: {}".format(root.find("ERROR").text))

        # Some samples such as SAMN13241871 are linked to multiple runs e.g. SRR10489833
        accessions_set = None if accessions is None else set(accessions)

        for pkg in root.findall('EXPERIMENT_PACKAGE'):
            d = collections.OrderedDict()
            d['experiment_accession'] = try_get(lambda: pkg.find('./EXPERIMENT').attrib['accession'])
            d['experiment_title'] = try_get(lambda: pkg.find('./EXPERIMENT/TITLE').text)
            l = pkg.find('./EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR')
            d['library_name'] = try_get(lambda: l.find('LIBRARY_NAME').text)
            d['library_strategy'] = try_get(lambda: l.find('LIBRARY_STRATEGY').text)
            d['library_source'] = try_get(lambda: l.find('LIBRARY_SOURCE').text)
            d['library_selection'] = try_get(lambda: l.find('LIBRARY_SELECTION').text)
            d['library_layout'] = try_get(lambda: l.find('LIBRARY_LAYOUT')[0].tag)
            d['platform'] = try_get(lambda: pkg.find('./EXPERIMENT/PLATFORM')[0].tag)
            d['model'] = try_get(lambda: pkg.find('./EXPERIMENT/PLATFORM/')[0].text)
            d['submitter'] = ''
            for k, v in pkg.find('./SUBMISSION').attrib.items():
                if k not in ('accession','alias'):
                    if d['submitter'] == '':
                        d['submitter'] = v
                    else:
                        d['submitter'] = "{}, {}".format(d['submitter'], v)
            d[STUDY_ACCESSION_KEY] = try_get(lambda: pkg.find('./STUDY').attrib['accession'])
            d[BIOPROJECT_ACCESSION_KEY] = try_get(lambda: pkg.find('./STUDY/IDENTIFIERS/EXTERNAL_ID[@namespace="BioProject"]').text)
            d['study_alias'] = try_get(lambda: pkg.find('./STUDY').attrib['alias'])
            d['study_centre_project_name'] = try_get(lambda: pkg.find('./STUDY/DESCRIPTOR/CENTER_PROJECT_NAME').text)
            d['organisation'] = try_get(lambda: pkg.find('./Organization/Name').text)
            d['organisation_department'] = try_get(lambda: pkg.find('./Organization/Address/Department').text)
            d['organisation_institution'] = try_get(lambda: pkg.find('./Organization/Address/Institution').text)
            d['organisation_street'] = try_get(lambda: pkg.find('./Organization/Address/Street').text)
            d['organisation_city'] = try_get(lambda: pkg.find('./Organization/Address/City').text)
            d['organisation_country'] = try_get(lambda: pkg.find('./Organization/Address/Country').text)
            first_name = try_get(lambda: pkg.find('./Organization/Contact/Name/First').text)
            last_name = try_get(lambda: pkg.find('./Organization/Contact/Name/Last').text)
            d['organisation_contact_name'] = "{} {}".format(first_name, last_name)
            d['organisation_contact_email'] = try_get(lambda: pkg.find('./Organization/Contact').attrib['email'])
            # Remove carriage return from sample description e.g. for SRR1263047
            d['sample_description'] = try_get(lambda: pkg.find('./SAMPLE/DESCRIPTION').text.replace('\r',''))
            d['sample_alias'] = try_get(lambda: pkg.find('./SAMPLE').attrib['alias'])
            d['sample_accession'] = try_get(lambda: pkg.find('./SAMPLE').attrib['accession'])
            d['biosample'] = try_get(lambda: pkg.find('./SAMPLE/IDENTIFIERS/EXTERNAL_ID[@namespace="biosample"]').text)
            d['sample_title'] = try_get(lambda: pkg.find('./SAMPLE/TITLE').text)
            d['taxon_name'] = try_get(lambda: pkg.find('./SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME').text)
            sample_sample_name = None
            sample_title = None
            if pkg.find('./SAMPLE/SAMPLE_ATTRIBUTES'):
                for attr in pkg.find('./SAMPLE/SAMPLE_ATTRIBUTES'):
                    tag_el = attr.find('TAG')
                    value_el = attr.find('VALUE')

                    # Some samples like
                    # https://www.ncbi.nlm.nih.gov/sra?term=ERS1240061&report=FullXml
                    # have entries with a tag, but no value. Ignore these.
                    if tag_el is not None and value_el is not None:
                        tag = tag_el.text
                        value = value_el.text
                        if tag == 'Title':
                            sample_title = value
                        elif tag == 'sample name':
                            sample_sample_name = value
                        d[tag] = value
            if sample_sample_name is not None:
                d[SAMPLE_NAME_KEY] = sample_sample_name
            elif sample_title is not None:
                d[SAMPLE_NAME_KEY] = sample_title
            else:
                d[SAMPLE_NAME_KEY] = d['library_name'] #default, maybe there's always a title though?
            d['study_title'] = try_get(lambda: pkg.find('./STUDY/DESCRIPTOR/STUDY_TITLE').text)
            d['design_description'] = try_get(lambda: pkg.find('./EXPERIMENT/DESIGN/DESIGN_DESCRIPTION').text)
            d['study_abstract'] = try_get(lambda: pkg.find('./STUDY/DESCRIPTOR/STUDY_ABSTRACT').text)
            study_links_xrefs = try_get(lambda: pkg.find('./STUDY/STUDY_LINKS'))
            if study_links_xrefs is not None:
                # Convert db to lower case because otherwise have PUBMED and pubmed e.g. ERR1914274 and SRR9113719
                study_links = list([
                    {'db': x.find('DB').text.lower(), 'id': x.find('ID').text}
                    for x in study_links_xrefs.findall('STUDY_LINK/XREF_LINK')])
                # Record URL links like SRR7051324
                for x in study_links_xrefs.findall('STUDY_LINK/URL_LINK'):
                    study_links.append({
                        'label': x.find('LABEL').text,
                        'url': x.find('URL').text
                    })
                d['study_links'] = json.dumps(study_links)
            else:
                d['study_links'] = json.dumps([])
            
            # Account for the fact that multiple runs may be associated with
            # this sample
            d['number_of_runs_for_sample'] = len(pkg.findall('./RUN_SET/RUN'))

            for run in pkg.findall('./RUN_SET/RUN'):
                accession_here = run.attrib['accession']
                if accessions_set is None or accession_here in accessions_set:
                    d2 = d.copy()
                    d2['spots'] = try_get(lambda: int(run.attrib['total_spots']))
                    d2[BASES_KEY] = try_get(lambda: int(run.attrib['total_bases']))
                    d2['run_size'] = try_get(lambda: int(run.attrib['size']))
                    d2[RUN_ACCESSION_KEY] = try_get(lambda: run.attrib['accession'])
                    d2['published'] = try_get(lambda: run.attrib['published'])
                    stats = run.find('Statistics')
                    if stats is not None:
                        for (i, r) in enumerate(stats):
                            d2['read{}_length_average'.format(i+1)] = r.attrib['average']
                            d2['read{}_length_stdev'.format(i+1)] = r.attrib['stdev']
                    data_frames.append(d2)

        return pd.DataFrame(data_frames)

    def _print_xml(self, element, prefix):
        if prefix is None or prefix == '':
            p2 = ""
        else:
            p2 = "{}_".format(prefix)
        for k, v in element.attrib.items():
            print("\t".join(['{}{}'.format(p2,k),v]))
        if element.text:
            print("\t".join(['{}{}'.format(p2, element.tag), element.text]))
        for e in element:
            self.print_xml(e, '{}{}'.format(p2, e.tag))

    def efetch_sra_from_accessions(self, accessions):
        all_accessions = list(set(accessions))
        if len(all_accessions) == 0:
            return []

        metadata_chunks = []
        chunk_size = 500

        # Complicated calculation here.
        num_chunks = len(all_accessions)/chunk_size
        if num_chunks > int(num_chunks):
            num_chunks = int(num_chunks) + 1
        num_chunks = int(num_chunks)

        if len(all_accessions) > chunk_size:
            logging.info("Querying for {} accessions in {} chunks".format(len(all_accessions), num_chunks))
            accession_iter = tqdm(iterable_chunks(all_accessions, chunk_size), total=num_chunks)
        else:
            accession_iter = [all_accessions]

        for accessions_plus in accession_iter:
            accessions = [a for a in accessions_plus if a is not None]

            if num_chunks == 1:
                logging.info("Querying NCBI esearch for {} distinct accessions e.g. {}".format(
                    len(accessions), accessions[0]))
            sra_ids = []

            webenv = None
            request_term = ' OR '.join(["{}[accn]".format(acc) for acc in accessions])

            retmax = len(accessions)+10
            params=self.add_api_key({
                "db": "sra",
                "term": request_term,
                "tool": "kingfisher",
                "email": "kingfisher@github.com",
                "retmax": retmax,
                "usehistory": "y",
                })
            if webenv is None:
                params['WebEnv'] = webenv

            res = self._retry_request(
                "esearch from accessions", 
                lambda: requests.post(
                    url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                    data=params))

            root = ET.fromstring(res.text)
            if webenv is None:
                webenv = root.find('WebEnv').text
            id_list_node = root.find('IdList')
            sra_ids = list(set([c.text for c in id_list_node]))

            if len(sra_ids) == 0:
                logging.warning("Unable to find any accessions, from the list: {}".format(accessions))
                return None

            if num_chunks == 1:
                logging.info("Querying NCBI efetch for {} distinct IDs e.g. {}".format(
                    len(sra_ids), sra_ids[0]))
            metadata = self.efetch_metadata_from_ids(webenv, accessions, len(sra_ids))

            # Ensure all hits are found, and trim results to just those that are real hits
            if RUN_ACCESSION_KEY not in metadata.columns:
                raise Exception("No metadata could be retrieved")

            if len(metadata) != len(accessions):
                found_runs = set(metadata[RUN_ACCESSION_KEY].to_list())
                not_found = list([a for a in accessions if a not in found_runs])
                logging.warning("Unable to find all accessions. The {} missing ones were: {}".format(
                    len(not_found), not_found
                ))
            metadata_chunks.append(metadata)

        metadata = pd.concat(metadata_chunks)
        metadata.sort_values([STUDY_ACCESSION_KEY,RUN_ACCESSION_KEY], inplace=True)

        return metadata