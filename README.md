# ExpansionHunter_Classifier
Predict if the EH reported genotypes are True or False based on the read support for the reported alleles.

Developed and tested with anaconda3. It should work with python >= 3.7

## 1. Extract Features from EH vcf file

Use the R snippet to generate features from EH genotypes.

## 2. Calculate % of GC content and percentage of motif match with reference sequence

python ~/CreatSequenceComposition.py -h
usage: CreatSequenceComposition.py [-h] --TRLists TRLISTS --fasta FASTA [--out OUT]

optional arguments:
  -h, --help         show this help message and exit
  --TRLists TRLISTS  .List of TRs to check with classifier (tab delimited with columns chrom, start, end, motif)
  --fasta FASTA      .Reference genome fasta file
  --out OUT          .Prefix for output file (default: SequenceComposition)


## 3. Classify EH genotype as True/False expansions

python ~/EH_Classifier_LociPrediction.py -h
usage: EH_Classifier_LociPrediction.py [-h] --model MODEL --seqInfo SEQINFO --features FEATURES [--out OUT]

Classify list of Expansion Hunter Loci as True/False expansions

optional arguments:
  -h, --help           show this help message and exit
  --model MODEL        Trained model file 
  --seqInfo SEQINFO    Sequence composition metrics for the list of TRs to classify  (generated in step2)
  --features FEATURES  EH genotype features extracted from EH genotype table(generated step1)
  --out OUT            Prefix for output file (default: EH_Loci)



