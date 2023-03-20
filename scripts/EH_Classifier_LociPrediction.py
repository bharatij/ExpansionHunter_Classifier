# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 14:31:14 2021
@author    : Bharati Jadhav (github: bharatij)
Goal	   : Run classifer on Eexpansion Hunter Genoype Features 
"""
import argparse
import sys
import pickle
import pandas as pd
import numpy as np

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = argparse.ArgumentParser(description='Classify list of Expansion Hunter Loci as True/False expansions ')
    parser.add_argument(
        '--model', type=str, required = True,
        help='Trained model file')
    parser.add_argument(
        '--seqInfo', type=str, required = True,
        help='Sequence composition metrics for the list of TRs to classify')
    parser.add_argument(
        '--features', type=str, required = True,
        help='EH genotype features extracted from EH vcf')
    parser.add_argument(
        '--out', type=str, default = 'EH_Loci',
        help='Prefix for output file (default: %(default)s)')
    return parser.parse_args()

def main():
    args = parse_args()
    modelFile = args.model
    seqInfoFile = args.seqInfo
    featureFile = args.features
    dname = args.out

    model = pickle.load(open(modelFile, 'rb'))
    seqInfo = pd.read_csv(seqInfoFile, sep='\t')
    seqInfo = seqInfo[['chrom','start','end','CG','PerMotifMatch']]
    seqInfo = seqInfo.drop_duplicates()

    ### data predictions
    #PredictLabels(outprefix)
    print("Processing dataset ...." + dname)
    testdf = pd.read_csv(featureFile, sep='\t')
    testdf = testdf.merge(seqInfo[['chrom','start','end','CG','PerMotifMatch']], how='left')
    testdf['CG'] = testdf['CG'].replace(np.nan, 0)
    testdf['PerMotifMatch'] = testdf['PerMotifMatch'].replace(np.nan, 0)
    testdf['Ratio_IRR'] = (testdf['IRR_A1'] + 1)/(testdf['IRR_A2']+1)
    testdf['Ratio_SPR'] = (testdf['SPR_A1'] + 1)/(testdf['SPR_A2']+1)
    testdf['Ratio_FR'] = (testdf['FR_A1'] + 1)/(testdf['FR_A2']+1)
    testdf['MotifSize'] = testdf['RefUnit'].apply(lambda x: len(x))
    testdf['RefUnitInBP'] = testdf['RefCopy'] *  testdf['MotifSize']
    testdf['LongAllleleInBP'] = testdf['MotifSize'] *  testdf['LongAllele']
    conditions = [
        (testdf['LongAllele_Readtype'] == 'SPANNING'),
        (testdf['LongAllele_Readtype'] == 'INREPEAT')
        ]
    values = [testdf['LongAllele_SPR'], testdf['LongAllele_IRR']]
    testdf['LongAllele_IRR_SPR'] = np.select(conditions, values)
    testdf['RatioRefLongAllele'] = (testdf['RefUnitInBP'] )/(testdf['LongAllleleInBP'])
    testdf = testdf.reset_index()
    testSet = testdf[['Ratio_FR','LongAllele_IRR_SPR','LongAllleleInBP','RefUnitInBP','PerMotifMatch','CG']]
    testPreds = model.predict(testSet)
    res = testdf
    res['PredictedLabels'] = testPreds
    print("Writing predictions in output file ...." )

    fout= dname + '_ClassifierPredictions.tsv'
    res.to_csv(fout, sep= '\t', index = False, na_rep='NaN')


if __name__ == '__main__':
    main()
