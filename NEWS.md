## What's new in v1.2.0
- Added two functions to handle ASCAT files instead of ChAS `load_ascat` and
`workflow_oncoscan.ascat`.
- Main workflow function `workflow_oncoscan.run` is renamed 
`workflow_oncoscan.chas` for consistency with the ASCAT function
- Bug fixes in adjust_loh (crashed with segments of length 1 or if no LOH 
segments where present)

## What's new in v1.1.0

- To simplify the workflow, the gender of the patient is not taken into account 
anymore. That implies that in a male sample, a gain of 3 extra copies on the X 
or Y chromosome is considered as a gain and not an amplification anymore. For 
female samples, nothing changes.
- The oncoscan coverage has been corrected to reflect only areas where there are
groups of probes. Isolated probes where causing issues to identify arm-level 
alterations as ChAS segments where never extended to these probes and the 90% 
threshold could never be met (particularly on chromosomal arms 9p and Yq).
- Minor corrections in vignette

## What's new in v1.0.0

- The nLST test has been clinically validated on 384 patients from the PAOLA-1
    trial and the recommended threshold is now >=15.
- The default value for arm-level alterations has been set to 90% as mentioned
    in the publication [Christinat et al., J Mol Diagn 2021; PMID: 34454110].
- The genomic LOH score (percent of LOH bases) has been added; `score_gloh`.
- Adds a flag "no tumor?" if the percentage of altered bases is less than 1%.
- Package to be released on Bioconductor


## What's new in v0.2.0

- Novel HRD score (nLST: number of LSTs, normalized by ploidy): `score_nlst`
- Change in Oncoscan workflow to use the nLST score and thresholds.
- New function to compute the number of Mb altered (with or without LOH):
    `score_mbalt`
- New function to compute the average copy number: `score_avgcn`
- New function to estimate the number of whole-genome duplication events (based
    on the average copy number and the thresholds defined by Carter et al.): 
    `score_estwgd`
