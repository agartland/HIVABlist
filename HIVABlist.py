import pandas as pd
import numpy as np
import os.path as op

from HLAPredCache import hlaPredCache

"""TODO:
    (1) These paths could be updated so that they are relative to the package."""

__all__ = ['loadEpitopes', 'loadPredictions']

def loadEpitopes(dataPath='data/'):
    """Load epitope lists and existing predictions."""
    ST = {'B70':'B*1503',
          'B63':'B*1516',
          'A8': 'B*0801',
          'B5': 'B*5101',
          'B62': 'B*1501',
          'A19': 'A*2901',
          'A28': 'A*6802',
          'B61': 'B*4002',
          'B17': 'B*5701',
          'B60': 'B*4001'}
    rep1to2 = ['A*24',
               'B*44',
               'B*07',
               'C*12',
               'C*15',
               'C*01',
               'C*06',
               'C*14',
               'C*07',
               'C*02']
    repD = {'A*0801':'B*0801',
            'A*02.01':'A*0201',
            'A2.1':'A*0201',
            'A11.1':'A*1101',
            'B*7001':'B*1503',
            'B*6001':'B*4001'}

    dropHLAs = []

    def _loadOneSet(fn):
        df = pd.read_csv(op.join(dataPath, fn))
        print '%s version: "%s"' % (fn, df.iloc[0,0])
        df = df.drop(0, axis=0)
        df = df.loc[(df.Epitope.str.len() >= 8) & (df.Epitope.str.len() <= 15)]
        df = _addRows(df)

        df['oldHLA'] = df['HLA'].copy()
        df['HLA'] = df['HLA'].map(_fixAllele)
        df = df.dropna(subset=['HLA'], axis=0)
        df = df.drop_duplicates()
        return df

    def _addRows(df):
        """Add a row for each individual HLA allele"""
        outDf = df.copy()
        addRows = []
        for i,row in df.iterrows():
            if isinstance(row['HLA'], basestring) and ',' in row['HLA']:
                hlas = row['HLA'].split(',')
                outDf.loc[i,'HLA'] = hlas[0]
                for h in hlas[1:]:
                    tmpRow = row.copy()
                    tmpRow['HLA'] = h
                    addRows.append(tmpRow)
        outDf = outDf.append(addRows, ignore_index=False)
        outDf.index = np.arange(outDf.shape[0])
        return outDf

    def _fixAllele(h):
        if isinstance(h, basestring):
            h = h.replace('supertype', '')
            h = h.replace('W', 'w')
            h = h.replace('Cw*', 'C*')
            h = h.replace('Cw', 'C')
            h = h.replace('?', '')
            h = h.strip()
            
            if not h[0] in 'ABC':
                return np.nan

            if h[1] == '*':
                if len(h) == 4:
                    h += '01'
            else:
                h = ST.get(h,h)
                if len(h)==2:
                    h = '%s*0%s01' % (h[0], h[1])
                elif len(h)==3:
                    h = '%s*%s01' % (h[0], h[1:])
            
            if h[-2:] == '01' and h[:4] in rep1to2:
                h = h[:4] + '02'
            h = repD.get(h, h)
            if h in dropHLAs:
                h = np.nan
        return h

    adf = _loadOneSet(fn='optimal_ctl_summary.csv')
    bdf = _loadOneSet(fn='ctl_summary.csv')
    return adf, bdf

def loadPredictions(adf, bdf, predFile='data/predictions.csv'):
    ba = hlaPredCache()
    ba.addFromFile(predFile)

    """Check that all predictions are present."""
    missingA = [pep for pep,hla in zip(adf.Epitope, adf.HLA) if np.isnan(ba[(hla,pep)])]
    missingB = [pep for pep,hla in zip(bdf.Epitope, bdf.HLA) if np.isnan(ba[(hla,pep)])]

    if len(missingA) + len(missingB) == 0:
        print "Loaded predictions for all epitopes."
    else:
        print "Missing %d predictions." % (len(missingA) + len(missingB))

    adf['log_IC50'] = [ba[(row['HLA'], row['Epitope'])] for i,row in adf.iterrows()]
    bdf['log_IC50'] = [ba[(row['HLA'], row['Epitope'])] for i,row in bdf.iterrows()]
    return adf, bdf, ba