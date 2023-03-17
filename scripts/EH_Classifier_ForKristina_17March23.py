# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 14:31:14 2021

@author: jadhab01
"""

import sys
import pickle
import pandas as pd
import numpy as np

modelFile = sys.argv[1]
seqInfoFile = sys.argv[2]
featureFile = sys.argv[3]
cohort = sys.argv[4]

def PredictLabels(dname):
    print("Processing dataset ...." + dname)
    testdf = pd.read_csv(featureFile, sep='\t')
    testdf = testdf.merge(seqInfo[['VARID','CG','PerMotifMatch']], how='left')  
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
    fout= dname + '_ClassifierPredictions.tsv'
    res.to_csv(fout, sep= '\t', index = False, na_rep='NaN')
    return(print("Done Processing dataset ...." + dname))

########## Input processing ##########
model = pickle.load(open(modelFile, 'rb'))
seqInfo = pd.read_csv(seqInfoFile, sep='\t')
seqInfo = seqInfo[['Locus','CG','PerMotifMatch']]
seqInfo = seqInfo.drop_duplicates()
seqInfo=seqInfo.rename(columns={'Locus':'VARID'})

### data predictions
PredictLabels(cohort)

