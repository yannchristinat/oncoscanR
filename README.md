# oncoscanR
author: Yann Christinat

date: 13 june 2022

version: 1.0.1

## Description
OncoscanR is an R package to handle Copy Number Variation analyses originating from the Oncoscan assay (Affymetrix). It
allows computation of different homologous recombination deficiency (HRD) scores and the tandem duplication plus score
(TDplus) to identify CDK12-mutated tumors [Popova et al., Cancer Res 2016]. The package also allows for identification
of arm-level alterations (e.g. gain of chromosome arm 1p).

**IMPORTANT**: The package expects as input the text exported file from ChAS (Chromosome Analysis Suite; the Affymetrix
software to identify CNV segments from the Oncoscan Assay). The package assumes that all segments given in the file are
correct and true. The ChAS text file has to contain the columns `Type`, `CN State` and `Full Location` (to setup in
ChAS). Any text file that complies with this structure should work equally well.

Note that the Oncoscan does not cover the p arms of chromosome 13, 14, 15 and 22. The coverage on the p arm of
chromosome 21 is only partial and is not included in the standard Oncoscan workflow (function `workflow_oncoscan.run` or 
script `bin/oncoscan-workflow.R`).

### Computation of arm-level alteration
An arm is declared globally altered if more than 90% of its bases are altered with a similar CNV type (amplifications
[3 extra copies or more], gains [1-2 extra copies], losses or copy-neutral losses of heterozygozity [LOH])[Christinat, 
J Mol Diagn 2021; PMID: 34454110]. For
instance, "gain of 3p" indicates that there is more than 90% of arm with 3 copies but less than 90% with 5 (otherwise
it would be an amplification). Prior to computation, segments of same copy number and at a distance <300Kbp (Oncoscan
resolution genome-wide) are merged. The remaining segments are filtered to a minimum size of 300Kbp.


### HRD scores
#### Score nLST
HRD score developed at HUG and based on the LST score by Popova et al. but normalized by an estimation of the number of
whole-genome doubling events.Of note, copy-neutral LOH segments are removed before computation.

`nLST = LST - 7*W/2` where `W` is the number of whole-genome doubling events.

The score is positive if there are at least 15 nLST.

#### LST
Procedure based on the paper from Popova et al, Can. Res. 2012 (PMID: 22933060). First segments
smaller than 3Mb are removed, then segments are smoothed with respect to copy number at a distance of 3Mb.
The number of LSTs is the number of breakpoints (breakpoints closer than 3Mb are merged) that have a segment
larger or equal to 10Mb on each side. This score was linked to BRCA1/2-deficient tumors.

#### LOH
Procedure based on the paper from Abkevich et al., Br J Cancer 2012 (PMID: 23047548). 
Number of LOH segments larger than 15Mb but excluding segments on chromosomes with a global LOH alteration. 
This score was linked to BRCA1/2-deficient tumors.

#### gLOH
The percentage genomic LOH score is computed as described in the FoundationFocus CDx BRCA
LOH assay; i.e. the percentage of bases covered by the Oncoscan that display a loss of heterozygosity
independently of the number of copies, excluding chromosomal arms that have a global LOH
(>=90 arm length). To compute with the armlevel_alt function on LOH segments only). This
score was linked to BRCA1/2-deficient tumors.

### Score TDplus
Procedure based on the paper from Popova et al., Cancer Res 2016 (PMID: 26787835). The TDplus
score is defined as the number of regions larger than 1Mb but smaller or equal to 10Mb with a gain of one
or two copies. This score was linked to CDK12-deficient tumors. 
They also identified as second category of tandem duplication whose size is smaller or equal than 1Mb and around 
300Kb but could not link it to a phenotype. Note that due to its resolution the Oncoscan assaywill most likely miss t
his second category. Nonetheless it is reported by the function but not by the standard workflow.

## Installation
There are two options to install the package: 
1. Download the `oncoscanR_1.0.0.tar.gz` file (stable version). Then in R, set the working directory to where the
compressed package is and run `install.packages('oncoscanR_1.0.0.tar.gz', repos=NULL, type='source')`.
2. In R, install the devtools package (`install.packages('devtools')`), load it (`library(devtools)`), then run
`install_github('yannchristinat/oncoscanR')`.

The package requires the prior installation of the packages `GenomicRanges` (bioconductor), `magrittr`, `jsonlite` and
`readr`.

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicRanges")

install.packages(c("magrittr", "jsonlite", "readr"))
```

## Testing the installation
Open R and type the following commands:
```{r}
library(oncoscanR)
segs.filename <- system.file("extdata", "chas_example.txt", package = "oncoscanR")
workflow_oncoscan.run(segs.filename, "M")
```

If everything is setup fine, it should return a list with no arm-level alterations and a negative HRD score (nLST=1).


## Usage
The main workflow can be launched either in R via the `workflow_oncoscan.run(chas.fn, gender)` function or via the
script `bin/run_oncoscan_workflow.R`:

Usage: `Rscript path_to_oncoscanR_package/bin/oncoscan-workflow.R CHAS_FILE GENDER`
- `CHAS_FILE`: Path to the text export file from ChAS or a compatible text file.
- `GENDER`: Gender of the sample (used to handle sex chromosomes). Has to be M (male) or F (female).

The script will output a JSON string into the terminal with all the computed information. :

```{json}
{
  "armlevel": {
    "AMP": [],
    "LOSS": ["17p", "2q", "4p"],
    "LOH": ["14q", "5q", "8p", "8q"],
    "GAIN": [19p", "19q", "1q", "20p", "20q", "3q", "5p", "6p", "9p", "9q", "Xp", "Xq"]
  },
  "scores": {
    "HRD": "Negative, nLST=12",
    "TDplus": 22,
    "avgCN": "2.43"
  },
  "gender": "F",
  "file": "H19001012_gene_list_full_location.txt"
}
```

Please read the manual for a description of all available R functions.

## References
1. "Homologous Recombination Deficiency (HRD) Score Predicts Response to Platinum-Containing Neoadjuvant Chemotherapy
in Patients with Triple-Negative Breast Cancer.", M. Telli et al., Clin Cancer Res volume 22(15), august 2016.
2. "Ovarian Cancers Harboring Inactivating Mutations in CDK12 Display a Distinct Genomic Instability Pattern
Characterized by Large Tandem Duplications.", T. Popova et al., Cancer Res volume 76(7), april 2016.
3. "Ploidy and large-scale genomic instability consistently identify basal-like breast carcinomas with BRCA1/2
inactivation.", T. Popova et al., Cancer Res volume 72(21), november 2012.
4. "Patterns of genomic loss of heterozygosity predict homologous recombination repair defects in epithelial ovarian
cancer.", V. Abkevich et al., Br J Cancer. 2012 Nov 6;107(10).
5. "Absolute quantification of somatic DNA alterations in human cancer", S. Carter et al., Nat Biotech, 2012 volume
30(5).
6. "Automated Detection of Arm-Level Alterations for Individual Cancer Patients in the Clinical Setting", Y Christinat 
et al., J Mol Diagn 2021, Dec;23(17):1722-1731.

