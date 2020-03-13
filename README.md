# oncoscanR

## Description
OncoscanR is an R package to handle Copy Number Variation analyses originating from the Oncoscan assay (Affymetrix). It allows computation of two homologous recombination default (HRD) scores, LST and HR-LOH as defined by Telli et al. [Clin Cancer Res 2016], along with the tandem duplication plus score (TDplus) to identify CDK12-mutated tumors [Popova et al., Cancer Res 2016].
The package also allows for identification of arm-level alterations (i.e. gain of chromosome arm 1p). An arm is declared globally altered if more than 80% of its bases are altered with a similar CNV type (amplifications [3 extra copies or more], gains [1-2 extra copies], losses or copy-neutral losses of heterozygozity [LOH]). For instance, "gain of 3p" indicates that there is more than 80% of arm with 3 copies but less than 80% with 5 (otherwise it would be an amplification). Prior to computation segments of same copy number and at a distance <300Kbp (Oncoscan resolution genome-wide) are merged. The remaining segments are filtered to a minimum size of 300Kbp.

Note that the Oncoscan does not cover the p arms of chromosome 15, 16 and 21.

**IMPORTANT**: The package expects as input the text exported file from ChAS (Chromosome Analysis Suite; the Affymetrix software to identify CNV segments from the Oncoscan Assay). The package assumes that all segments given in the file are correct and true. The ChAS text file has to contain the columns `Type`, `CN State` and `Full Location` (to setup in ChAS).

## Installation
There are two options to install the package: 
1. Download the `oncoscanR_0.1.0.tar.gz` file (stable version). Then in R, set the working directory to where the compressed package is and run `install.packages('oncoscanR_0.1.0.tar.gz', repos=NULL, type='source')`.
2. In R, install the devtools package (`install.packages('devtools')`), load it (`library(devtools)`), then run `install_github('yannchristinat/oncoscanR')`.

The package requires the prior installation of the packages `GenomicRanges` (bioconductor), `magrittr`, `jsonlite` and `readr`.

## Testing the installation
Open R and type the following commands:
`library(oncoscanR)`
`segs.filename <- system.file("extdata", "chas_example.txt", package = "oncoscanR")`
`workflow_oncoscan.run(segs.filename, "M")`

If everything is setup fine, it should return a list with arm-level alterations and HRD scores.

## Usage
The main workflow can be launched either in R via the `workflow_oncoscan.run(chas.fn, gender)` function or via the script `bin/run_oncoscan_workflow.R`:

Usage: `Rscript path_to_oncoscanR_package/bin/oncoscan-workflow.R CHAS_FILE GENDER`
- `CHAS_FILE`: Path to the text export file from ChAS or a compatible text file.
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
