# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 20:45:29 2021
Calculate the A, T, C and G, seq purity and GC content
@@author    : Bharati Jadhav (github: bharatij)
@author: jadhab01
"""
import argparse
import sys
import pandas as pd
import numpy as np
import sys
import subprocess

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = argparse.ArgumentParser(description='Combine multiple TR lists from different sources for EH Catalog')
    parser.add_argument(
        '--TRLists', type=str, nargs='+', required = True,
        help='.List of files with TRs to check with classifier (tab delimited with columns chrom, start, end, motif)
    parser.add_argument(
        '--fasta', type=str, required = True,
        help='.Reference genome fasta file')
    parser.add_argument(
        '--out', type=str, default = 'SequenceComposition',
        help='.Prefix for output file (prefix will be TRList) (default: %(default)s)')
    return parser.parse_args()


def parse_TRLists(filename):
    """Parse all TR list"""
    try:
        trdf = pd.read_csv(filename, sep='\t' ,low_memory=False )
        trdf.columns = ['chrom', 'start', 'end', 'motif']
    except pd.io.common.EmptyDataError:
        sys.exit('ERROR: file {0} was empty.\n'.format(filename))
    return(trdf)

def getFastaSeq(fasta, inBed, outBed):
    cmd = 'ml bedtools \n '
    subprocess.call(cmd, shell=True)
    cmd= cmd + 'bedtools getfasta -fi '+ fasta +  ' -bed ' +  inBed + ' -tab > ' + outBed
    subprocess.call(cmd, shell=True)
    seqDf = pd.read_csv(outBed, sep='\t', header=None)
    seqDf.columns =['Locus', 'Seq']
    seqDf[['chrom','tmp']] = seqDf.Locus.str.split(":",expand=True)
    seqDf['start'] = seqDf.tmp.str.split("-").str[0].astype(int)
    seqDf['end'] = seqDf.tmp.str.split("-").str[1].astype(int)
    seqDf=seqDf.drop(columns=['tmp','Locus'])
    cmd = 'rm ' + outBed + ' '+ inBed
    subprocess.call(cmd, shell=True)
    return(seqDf)

def getRefPerMatch(Seq, motif):
    Seq = Seq.lower()
    motif = motif.lower()
    nMatch = Seq.count(motif)
    perM = round(((nMatch * len(motif)) / len(Seq)) * 100,2)
    return perM

def main():
    args = parse_args()
    trfile = args.TRLists
    base_filename = args.out
    fasta = args.fasta
    inBed= base_filename + 'tempIn.bed'
    outBed = base_filename + 'tempOut.bed'
    print('Reading list of TR files.....')
    allStr = parse_TRLists(trfile)
    allStr = allStr.sort_values(by=['chrom', 'start'],ignore_index=True)
  
    trBed = allStr[['chrom','start','end']]
    trBed = trBed.drop_duplicates()
    trBed.to_csv(inBed, sep= '\t', index = False, na_rep='NaN',header=False )
        
    seqDf = getFastaSeq(fasta, inBed, outBed)
    allStrAnn = allStr.merge(seqDf, on=['chrom', 'start','end'])
    allStrAnn[['A']] = allStrAnn['Seq'].str.count('a|A')
    allStrAnn[['T']] = allStrAnn['Seq'].str.count('t|T')
    allStrAnn[['C']] = allStrAnn['Seq'].str.count('c|C')
    allStrAnn[['G']] = allStrAnn['Seq'].str.count('g|G')
    allStrAnn[['CG']] = (allStrAnn['C'] + allStrAnn['G'])/allStrAnn['Seq'].str.len()

    allStrAnn = allStrAnn.assign( PerTRRefMotifMatch = allStrAnn.apply(lambda x: getRefPerMatch(x['Seq'],x['motif']), axis=1))
    
    fout =   base_filename + '.tsv'
    allStrAnn.to_csv(fout, sep= '\t', index = False, na_rep='NaN')

if __name__ == '__main__':
    main()
