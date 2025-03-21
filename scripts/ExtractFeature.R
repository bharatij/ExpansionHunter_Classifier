#R code snippet to  extract details to create feature table
#df is the datafram with EH genotypes for samples/loci to check with following columns extracted from the EH vcf file
#chrom,start,end,RefCopy,RefUnit,SampleId,VARID,ReadType,Alleles,SpanningRead,FlankingRead,InRepeatRead,LocusCoverage,CI

library(tidyverse)
library(data.table)

df =  fread('EH_GenotypeTableToClassify.tsv', sep = "\t", header=T, check.names = F,nThread=8)
###Alleles
df [,A1 := sapply(strsplit(as.character(Alleles), "\\/"), '[', 1)]
df [,A2 := sapply(strsplit(as.character(Alleles), "\\/"), '[', 2)]
df[ ,c("A1","A2") :=  lapply(.SD, as.integer),.SDcols = c("A1","A2")]
###LongAllele
df[,LongAllele :=  pmax(A1,A2, na.rm = T )]

###IRR
  df [,IRR_A1 := sapply(strsplit(as.character(InRepeatRead), "\\/"), '[', 1)]
  df [,IRR_A2 := sapply(strsplit(as.character(InRepeatRead), "\\/"), '[', 2)]
  df[ ,c("IRR_A1","IRR_A2") :=  lapply(.SD, as.integer),.SDcols = c("IRR_A1","IRR_A2")]
 
  ###Spanning
  df [,SPR_A1 := sapply(strsplit(as.character(SpanningRead), "\\/"), '[', 1)]
  df [,SPR_A2 := sapply(strsplit(as.character(SpanningRead), "\\/"), '[', 2)]
  df[ ,c("SPR_A1","SPR_A2") :=  lapply(.SD, as.integer),.SDcols = c("SPR_A1","SPR_A2")]
 
  ###Flanking
  df [,FR_A1 := sapply(strsplit(as.character(FlankingRead), "\\/"), '[', 1)]
  df [,FR_A2 := sapply(strsplit(as.character(FlankingRead), "\\/"), '[', 2)]
  df[ ,c("FR_A1","FR_A2") :=  lapply(.SD, as.integer),.SDcols = c("FR_A1","FR_A2")]
 
  ###Read Type
  df [,Readtype_A1 := sapply(strsplit(as.character(ReadType), "\\/"), '[', 1)]
  df [,Readtype_A2 := sapply(strsplit(as.character(ReadType), "\\/"), '[', 2)]
 
  ###long allele spanning and IRR
  df [,LongAllele_Readtype := ifelse((!is.na(A2) & A1 <= A2), Readtype_A2, Readtype_A1)]
  df[,LongAllele_IRR := ifelse(!is.na(A2) & A1 <= A2, IRR_A2, IRR_A1)]
  df[,LongAllele_SPR := ifelse(!is.na(A2) & A1 <= A2, SPR_A2, SPR_A1)]
  df[,LongAllele_FR := ifelse(!is.na(A2) & A1 <= A2, FR_A2, FR_A1)]

  df[,A1temp := A1]
  df[,A2temp := A2]
  df[is.na(A2), A2 := A1[is.na(A2)]]
  df[is.na(A2temp),A1 := 0]
 
  df[,A1temp := IRR_A1]
  df[,A2temp := IRR_A2]
  df[is.na(IRR_A2), IRR_A2 := IRR_A1[is.na(IRR_A2)]]
  df[is.na(A2temp),IRR_A1 := 0]
 
  df[,A1temp := SPR_A1]
  df[,A2temp := SPR_A2]
  df[is.na(SPR_A2), SPR_A2 := SPR_A1[is.na(A2temp)]]
  df[is.na(A2temp),SPR_A1 := 0]
 
  df[,A1temp := FR_A1]
  df[,A2temp := FR_A2]
  df[is.na(FR_A2),FR_A2 := FR_A1[is.na(A2temp)]]
  df[is.na(A2temp),FR_A1 := 0]
 
  df[ ,':='(A1temp = NULL, A2temp = NULL)]

df= df[,.SD, .SDcols = c( 'chrom', 'start', 'end', 'RefCopy', 'RefUnit','VARID','SampleId','LongAllele', 'IRR_A1', 'IRR_A2', 'SPR_A1', 'SPR_A2', 'FR_A1', 'FR_A2', 'LongAllele_Readtype', 'LongAllele_IRR', 'LongAllele_SPR', 'LongAllele_FR')]
fwrite(df, file='EH_GenotypeFeaturesToUseWithClassifier.tsv', sep = "\t", quote = F, row.names = F,nThread=8)
####after running this code your table shold have following columns
SampleId, VARID,chrom, start, end, RefCopy, RefUnit, LongAllele, IRR_A1, IRR_A2, SPR_A1, SPR_A2, FR_A1, FR_A2, LongAllele_Readtype, LongAllele_IRR, LongAllele_SPR, LongAllele_FR
