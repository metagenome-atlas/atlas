import os, shutil
from pathlib2 import Path
import pandas as pd

test_directory= Path("test/Test_sra_init")
test_directory.mkdir( exist_ok=True)

Run_info_table= str( test_directory / "Runinfo.tsv" )
id="PRJEB20796"




from atlas.init.get_SRA_runinfo import get_runtable_from_ids
from atlas.init.parse_sra import filter_runinfo, validate_runinfo_table


get_runtable_from_ids(id,output_file=Run_info_table )




RunTable= pd.read_csv(Run_info_table, sep='\t', index_col=0)

validate_runinfo_table(RunTable)

RunTable_filtered= filter_runinfo(RunTable)

RunTable_filtered.BioSample.to_csv( str( test_directory / "Run2Sample.tsv" ) , sep= '\t')


if RunTable.shape[0] != RunTable.BioSample.unique().shape[0] :

    logger.debug("I will automatically merge runs from the same biosample.")







