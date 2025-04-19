# ExpansionHunter_Classifier
Predict if the EH reported genotypes are True or False based on the read support for the reported alleles.

Developed and tested with Anaconda3. It should work with Python >= 3.7.

## 1. Extract Features from EH VCF file

Use the R script to generate features from the EH genotypes table extracted from the EH VCF file. This table includes the following columns extracted from the EH vcf.
chrom, start, end, Ref, Alt, Filter, RefCopy, RefUnit, SampleId, VARID, REPID, GT, ReadType, Alleles, SpanningRead, FlankingRead, InRepeatRead, LocusCoverage, CI

## 2. Calculate % of GC content and the percentage of motif match with a reference sequence

python ~/CreatSequenceComposition.py -h

usage: CreatSequenceComposition.py [-h] --TRLists TRLISTS --fasta FASTA [--out OUT]

optional arguments:
  -h, --help         show this help message and exit
  --TRLists TRLISTS  .List of TRs to check with classifier (tab delimited with columns chrom, start, end, motif)
  --fasta FASTA      .Reference genome fasta file
  --out OUT          .Prefix for output file (default: SequenceComposition)


## 3. Classify the EH genotype as True/False expansions
### NOTE: Uploaded the Random Forest classifier model trained using Python v3.12. If you are using Python v3.7, use EHClassifier_TrainedModel_pythonV3.7.pkl model and with Python v3.10 and up use EHClassifier_TrainedModel_pythonV3.10andUp.pkl

python ~/EH_Classifier_LociPrediction.py -h

usage: EH_Classifier_LociPrediction.py [-h] --model MODEL --seqInfo SEQINFO --features FEATURES [--out OUT]

Classify the list of Expansion Hunter Loci as True/False expansions

optional arguments:
  -h, --help           show this help message and exit
  --model MODEL        Trained model file 
  --seqInfo SEQINFO    Sequence composition metrics for the list of TRs to classify  (generated in step 2)
  --features FEATURES  EH genotype features extracted from EH genotype table(generated step 1)
  --out OUT            Prefix for output file (default: EH_Loci)


Cite the code: [![DOI](https://zenodo.org/badge/615451022.svg)](https://zenodo.org/doi/10.5281/zenodo.10821642)

