import pandas as pd
import numpy as np
import os.path as op
import operator
from functools import partial
import matplotlib.pyplot as plt

try:
    from HLAPredCache import hlaPredCache
    import statsmodels.api as sm
    import palettable
except ImportError:
    pass

DATA_FOLDER = op.join(op.dirname(op.abspath(__file__)), 'data')

__all__ = ['loadEpitopes',
           'loadPredictions',
           'plotEpitopeCDF',
           'isLANLEpitope']

def loadEpitopes(dataPath=DATA_FOLDER):
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
        print('%s version: "%s"' % (fn, df.iloc[0, 0]))
        df = df.drop(0, axis=0)
        df = df.loc[(df.Epitope.str.len() >= 8) & (df.Epitope.str.len() <= 15)]
        df = _addRows(df)

        df['oldHLA'] = df['HLA'].copy()
        df['HLA'] = df['HLA'].map(_fixAllele)
        df = df.dropna(subset=['HLA'], axis=0)
        df = df.drop_duplicates()
        df['ba_key'] = [(h, p) for h, p in zip(df.HLA, df.Epitope)]
        return df

    def _addRows(df):
        """Add a row for each individual HLA allele"""
        outDf = df.copy()
        addRows = []
        for i, row in df.iterrows():
            if isinstance(row['HLA'], str) and ',' in row['HLA']:
                hlas = row['HLA'].split(',')
                outDf.loc[i, 'HLA'] = hlas[0]
                for h in hlas[1:]:
                    tmpRow = row.copy()
                    tmpRow['HLA'] = h
                    addRows.append(tmpRow)
        outDf = outDf.append(addRows, ignore_index=False)
        outDf.index = np.arange(outDf.shape[0])
        return outDf

    def _fixAllele(h):
        if isinstance(h, str):
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
                h = ST.get(h, h)
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

def loadPredictions(adf, bdf, predFile=op.join(DATA_FOLDER, 'predictions.csv')):
    ba = hlaPredCache()
    ba.addFromFile(predFile)

    """Check that all predictions are present."""
    missingA = [pep for pep, hla in zip(adf.Epitope, adf.HLA) if np.isnan(ba[(hla, pep)])]
    missingB = [pep for pep, hla in zip(bdf.Epitope, bdf.HLA) if np.isnan(ba[(hla, pep)])]

    if len(missingA) + len(missingB) == 0:
        print("Loaded predictions for all epitopes.")
    else:
        print("Missing %d predictions." % (len(missingA) + len(missingB)))

    adf['log_IC50'] = [ba[(row['HLA'], row['Epitope'])] for i, row in adf.iterrows()]
    bdf['log_IC50'] = [ba[(row['HLA'], row['Epitope'])] for i, row in bdf.iterrows()]
    return adf, bdf, ba

def plotEpitopeCDF(df, ba):
    plt.clf()
    ic50Ticks = [20, 50, 200, 500, 1000, 2000, 5000, 10000, 20000]
    ECDF = sm.distributions.empirical_distribution.ECDF
    colors = palettable.colorbrewer.qualitative.Set1_3.mpl_colors
    for x, locus in enumerate('ABC'):
        obsInd = df.HLA.str.slice(stop=1) == locus
        ecdf = ECDF(df['log_IC50'].loc[obsInd])
        plt.step(ecdf.x, ecdf.y, '-', color=colors[x], lw=3)

        ref = [v for k, v in list(ba.items()) if k[0][0]==locus and not k in df.ba_key]
        ecdf = ECDF(ref)
        plt.step(ecdf.x, ecdf.y, '--', color=colors[x], alpha=1)

    plt.ylim((0, 1))
    plt.legend([plt.Line2D([0,0], [1,1], color=c) for c in colors], ['HLA-A', 'HLA-B', 'HLA-C'])
    plt.xlim((0, 15))
    plt.xticks(np.log(ic50Ticks), ic50Ticks)
    plt.xlabel('Predicted HLA binding $IC_{50}$ (nM)')
    plt.ylabel('Fraction of epitopes')


def _hamming_distance(str1, str2):
    """Hamming distance between str1 and str2."""
    assert len(str1) == len(str2), "Inputs must have the same length."
    return np.sum([i for i in map(operator.__ne__, str1, str2)])

def _optimalOverlap(pep1, pep2, minOverlap):
    #print()
    assert len(pep1) >= minOverlap
    assert len(pep2) >= minOverlap
    mn = None
    for i in range(len(pep2) - minOverlap, 0, -1):
        tmp2 = pep2[i:]
        L = min(len(tmp2), len(pep1))
        tmp1 = pep1[:L]
        tmp2 = tmp2[:L]
        for length in range(minOverlap, L + 1):
            d = _hamming_distance(tmp1[:length], tmp2[:length])
            if mn is None or d < mn:
                #print((tmp1[:length]))
                #print((tmp2[:length]))
                #print(d)
                mn = d

    for i in range(0, len(pep1) - minOverlap):
        tmp1 = pep1[i:]
        L = min(len(tmp1), len(pep2))
        tmp1 = tmp1[:L]
        tmp2 = pep2[:L]
        for length in range(minOverlap, L + 1):
            d = _hamming_distance(tmp1[:length], tmp2[:length])
            if mn is None or d < mn:
                #print((tmp1[:length]))
                #print((tmp2[:length]))
                #print(d)
                mn = d
    return mn

def _HLAissimilar(h1, h2):
    locus1, allele1 = h1.split('*')
    locus2, allele2 = h2.split('*')
    if locus1 == locus2 and allele1[:2] == allele2[:2]:
        return True
    else:
        return False

def _HLAissimilarIC50(h1, h2):
    pass


def isLANLEpitope(h, peptide, lanlDf=None, minOverlap=8, maxMM=4, bindingThresh=6.2, deltaBinding=1, twoDigitMatching=True):
    if lanlDf is None:
        lanlDf = loadEpitopes[1]
    dist = lanlDf.Epitope.map(partial(_optimalOverlap, pep2=peptide, minOverlap=minOverlap))
    tmpDf = lanlDf.loc[dist <= maxMM]
    hlaMatch = tmpDf.HLA.map(partial(_HLAissimilar, h2=h))
    return tmpDf.loc[hlaMatch]