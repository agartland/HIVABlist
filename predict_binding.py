import pandas as pd
import argparse

import HLAPredCache
from HIVABlist import loadEpitopes

def generatePredictionsFile(method='netmhcpan', outfile='data/predictions.csv'):
    adf, bdf = loadEpitopes()

    mers = list(set(adf.Epitope.tolist() + bdf.Epitope.tolist()))
    hlas = list(set(adf.HLA.tolist() + bdf.HLA.tolist()))

    ba = HLAPredCache.hlaPredCache()
    ba.addPredictions(method=method, hlas=hlas, peptides=mers, cpus=6, verbose=True)
            #except HLAPredCache.iedb_src.util.UnexpectedInputError:
                #print 'UnexpectedInput: Could not predict with method %s, allele "%s" (%dmers)' % (method, h, k)
            #except HLAPredCache.iedb_src.util.PredictorError:
                #print 'PredictorError: Could not predict with method %s, allele "%s" (%dmers)' % (method, h, k)

    ba.dumpToFile(outfile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate HLA binding predictions file for HIV A/B-list epitopes.')
    parser.add_argument('--method', metavar='METHOD', type = str, default = 'netmhcpan',
                       help='prediction method')
    parser.add_argument('--out', metavar='OUT_FILE', type = str, default = 'data/predictions.csv',
                       help='output filename')
    args = parser.parse_args()
    
    generatePredictionsFile(method=args.method, outfile=args.out)