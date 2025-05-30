# ExpansionHunter_Classifier
We developed a Random Forest classifier to recognize low-confidence TRE calls based on characteristics of the TR locus and the output metrics supporting each genotype generated by ExpansionHunter. We first trained this tool on a set of ~3,500 putative TREs from ExpansionHunter, for which we performed manual curation of read pileup plots to determine true and false positives. For each putative TRE, we generated a visualization of reads overlapping the TR locus using GraphAlignmentViewer and by manual inspection classified these as either true or false based on the number and quality of the read alignments supporting the putative TRE. In general, TREs considered as true were supported by multiple reads with few mismatches, whereas those considered as false were often supported by one or very few reads with multiple mismatches in the alignment. Please check our manuscript (https://www.nature.com/articles/s41588-024-01917-1) for the details of the performance of this model for reducing false positives.

This classsifier is developed and tested with Anaconda3. It should work with Python >= 3.7 and recently modified to work with Python v3.10 and above.

## Step 1. Extract Features from EH VCF file

Use the R script (scripts/ExtractFeature.R) to generate features for classifier from the EH genotypes table. The EH genotypes table contains columns extracted from the EH VCF file. 
### File Format
#### Input
1) The EH genotypes table is tsv format file which includes the following columns extracted from the EH vcf.

chrom, start, end, Ref, Alt, Filter, RefCopy, RefUnit, SampleId, VARID, REPID, GT, ReadType, Alleles, SpanningRead, FlankingRead, InRepeatRead, LocusCoverage, CI

2) List of sample and TR to check with classifier
   This tsv file includes two columns: SampleId and VARID
   
#### Output
The output tsv file contains the following columns describing features supporting  long allele which are used for predictions in classifier.

SampleId, VARID,chrom, start, end, RefCopy, RefUnit, LongAllele, IRR_A1, IRR_A2, SPR_A1, SPR_A2, FR_A1, FR_A2, LongAllele_Readtype, LongAllele_IRR, LongAllele_SPR, LongAllele_FR

## Step 2. Calculate % of GC content and the percentage of motif match with a reference sequence
Generate additional features such as sequence purity score and GC content in the reeference sequence.

### Quick Start
python ~/CreatSequenceComposition.py -h

usage: CreatSequenceComposition.py [-h] --TRLists TRLISTS --fasta FASTA [--out OUT]

optional arguments:
  -h, --help         show this help message and exit
  --TRLists TRLISTS  .List of TRs to check with classifier (tab delimited with columns chrom, start, end, motif)
  --fasta FASTA      .Reference genome fasta file
  --out OUT          .Prefix for output file (default: SequenceComposition)

### Tutorial
Example using TRs in (https://github.com/Illumina/ExpansionHunter/blob/master/variant_catalog/hg38/variant_catalog.json)  where variant_catalog_hg38.tsv contains columns: chrom, start, end, motif

CreatSequenceComposition.py --TRLists variant_catalog_hg38.tsv \
                            --fasta Homo_sapiens_assembly38.fasta \
                            --out variant_catalog_hg38_SequenceComposition 
### File Format
#### Input
The tsv file conaining the following columns for the list of TRs included in the EH json catalog file.
chrom   start   end     motif

#### Output
The outfile in tsv format containing collowing columns.

chrom   start   end     CG      PerMotifMatch

## 3. Classify the EH genotype as True/False expansions
This step uses the EH genotype features and Sequence Composition files to perform predictions using Random Forest classifier with the provided trained model in pickle file format.

##### NOTE: Added Random Forest classifier model trained using Python v3.12. If you are using Python v3.7, use EHClassifier_TrainedModel_pythonV3.7.pkl model and with Python v3.10 and up use EHClassifier_TrainedModel_pythonV3.10andUp.pkl

### Quick Start
python ~/EH_Classifier_LociPrediction.py -h

usage: EH_Classifier_LociPrediction.py [-h] --model MODEL --seqInfo SEQINFO --features FEATURES [--out OUT]

Classify the list of Expansion Hunter Loci as True/False expansions

optional arguments:
  -h, --help           show this help message and exit
  --model MODEL        Trained model file 
  --seqInfo SEQINFO    Sequence composition metrics for the list of TRs to classify  (generated in step 2)
  --features FEATURES  EH genotype features extracted from EH genotype table(generated step 1)
  --out OUT            Prefix for output file (default: EH_Loci)

### Tutorial
Example using Random Forest model trained with Python v3.12

python EH_Classifier_LociPrediction.py --model EHClassifier_TrainedModel_pythonV3.10andUp.pkl \
--seqInfo variant_catalog_hg38_SequenceComposition.tsv \
--features AllSamples_AllTRs_GT_SampleTRs_FeaturesToUseWithClassifier.tsv \
--out AllSamples_AllTRs_GT_SampleTRs

### File Format
#### Input
1) Trained model in the pickle file format provided in this repository
2) Feature generated from EH genotype using R script in Step 1
3) Sequence composition file generated for the list of TRs included in the EH catalog in Step2

#### Output
The output file with ClassifierPredictions.tsv suffix containing the column "PredictedLabels" with prediction labeled as T (True) and F (False)

Cite the code: [![DOI](https://zenodo.org/badge/615451022.svg)](https://zenodo.org/doi/10.5281/zenodo.10821642)

