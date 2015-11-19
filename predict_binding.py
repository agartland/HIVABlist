import pandas as pd

from HLAPredCache import hlaPredCache

from HIVABlist import loadEpitopes

adf, bdf = loadEpitopes()

mers = list(set(adf.Epitope.tolist() + bdf.Epitope.tolist()))
hlas = list(set(adf.HLA.tolist() + bdf.HLA.tolist()))

ba = HLAPredCache.hlaPredCache()
ba.addPredictions(method='netmhcpan', alleles=hlas, peptides=mers)
ba.dumpToFile('data/predictions.csv')