import pandas as pd
import os.path as op

from HLAPredCache import hlaPredCache

"""TODO:
    (1) These paths could be updated so that they are relative to the package."""

def loadEpitopes(dataPath=None):
    """Load epitope lists and existing predictions."""
    dataPath = GIT_PATH + 'HIVABlist/data/'
    adf = pd.read_csv(op.join(dataPath, 'optimal_ctl_summary.csv'))
    print 'A-list version: "%s"' % adf.iloc[0,0]
    adf = adf.drop(0, axis=0)

    bdf = pd.read_csv(op.join(dataPath, 'ctl_summary.csv'))
    print 'B-list version: "%s"' % bdf.iloc[0,0]
    bdf = bdf.drop(0, axis=0)

    def _addRows(df):
        """Add a row for each individual HLA allele"""
        outDf = df.copy()
        for i,row in df.iterrows():
            if ',' in row['HLA']:
                hlas = row['HLA'].split(',')
                outDf.loc[i,'HLA'] = hlas[0]
                for h in hlas[1:]:
                    tmpRow = row.copy()
                    tmpRow['HLA'] = h
                    outDf = outDf.append(tmpRow, ignore_index=True)
        return outDf

    adf = _addRows(adf)
    bdf = _addRows(bdf)

    def _fixAllele(h):
        return h

    return adf, bdf

def parseAllele(h):
    """Parse HLA alleles into standard 6 digits"""
    pass


def loadPredictions(dataPath=None, predFile=None):
    if dataPath is None:
        dataPath = GIT_PATH + 'HIVABlist/data/'
    if predFile is None:
        predFile = 'predictions.csv'
    
    adf,bdf = loadEpitopes()

    ba = hlaPredCache()
    ba.addFromFile(op.join(dataPath, predFile))

    """Check that all predictions are present."""
    missingA = [pep for pep,hla in zip(adf.Epitope, adf.HLA)]
    missingB = [pep for pep,hla in zip(bdf.Epitope, bdf.HLA)]

    return adf, bdf, ba