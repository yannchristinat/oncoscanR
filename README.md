# oncoscanR

## Description
OncoscanR is an R package to handle Copy Number Variation analyses originating from the Oncoscan assay (Affymetrix). It allows computation of two homologous recombination default (HRD) scores, LST and HR-LOH as defined by Telli et al. [Clin Cancer Res 2016], along with the tandem duplication plus score (TDplus) to identify CDK12-mutated tumors [Popova et al., Cancer Res 2016].
The package also allows for identification of arm-level alterations (i.e. gain of chromosome arm 1p). An arm is declared globally altered if more than 80% of its basees are altered with a similar CNV type (amplifications [3 extra copies or more], gains [1-2 extra copies], losses or copy-neutral losses of heterozygozity [LOH]).

**IMPORTANT**: The package expects as input the text exported file from ChAS (Chromosome Analysis Suite; the Affymetrix software to identify CNV segments from the Oncoscan Assay). The package assumes that all segments given in the file are correct and true. The ChAS text file has to contain the columns `Type`, `CN State` and `Full Location` (to setup in ChAS).

## Usage
The main workflow can be launched either in R via the `workflow_oncoscan.run(chas.fn, gender)` function or via the script "bin/run_oncoscan_workflow.R":

Usage: `Rscript path_to_oncoscanR_package/bin/oncoscan-workflow.R CHAS_FILE GENDER`
- `CHAS_FILE`: Path to the text export file from ChAS.
- `GENDER`: Gender of the sample (used to handle sex chromosomes). Has to be M (male) or F (female).

The script will output a JSON string into the terminal with all the computed information. :

`{
  "armlevel": {
    "AMP": [],
    "LOSS": ["17p", "2q", "4p"],
    "LOH": ["14q", "5q", "8p", "8q"],
    "GAIN": [19p", "19q", "1q", "20p", "20q", "3q", "5p", "6p", "9p", "9q", "Xp", "Xq"]
  },
  "scores": {
    "LST": 12,
    "LOH": 10,
    "TDplus": 66
  },
  "gender": "F",
  "file": "H19001012_gene_list_full_location.txt"
}`

Please read the manual for a description of all available R functions.

## References
1. "Homologous Recombination Deficiency (HRD) Score Predicts Response to Platinum-Containing Neoadjuvant Chemotherapy in Patients with Triple-Negative Breast Cancer.", M. Telli et al., Clin Cancer Res volume 22(15), august 2016.
2. "Ovarian Cancers Harboring Inactivating Mutations in CDK12 Display a Distinct Genomic Instability Pattern Characterized by Large Tandem Duplications.", T. Popova et al., Cancer Res volume 76(7), april 2016.
