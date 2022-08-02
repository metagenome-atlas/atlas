
import pandas as pd
import numpy as np


def split_orf_to_index(orf_names):
  """Split the typical prodigal orf name `Sample_contigNr_OrfNR` into three column dataset and conversts the numbers to the smallest unsigned number. 
  """
  result=  pd.Series(orf_series).str.rsplit("_",n=2,expand=True).rename(columns={0:"Sample",1:"ContigNr",2:"OrfNr"})
  result = result.apply(pd.to_numeric,errors="ignore",downcast='unsigned')

  are_numeric= result.dtypes.map(lambda x: np.issubdtype(x,np.number))
  assert all(are_numeric == np.array([False,True,True]) ) ,f"datatypes are not as expected {result.dtypes}"

  return result


def geneNr_to_string(GeneNrs, Ngenes=None):
  """Convert the array of gene number to the corresponding string, e.g.: 5 -> Gene0005
  The leading zeros depends on the number of genes
  """

  assert np.issubdtype(GeneNrs.dtype ,np.number)

  if Ngenes is None:
    Ngenes = GeneNrs.max()

  n_leading_zeros = len(str(Ngenes))

  number_format = f"Gene{{:0{n_leading_zeros}d}}"

  return GeneNrs.apply(number_format.format)
