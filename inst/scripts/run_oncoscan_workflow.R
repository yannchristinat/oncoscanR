#!/usr/bin/Rscript
# oncoscan-workflow.R
#
# Run the standard workflow for the detection of arm-level alterations and prints
# the output in JSON.
#
# Author: Yann Christinat
# Date: 1.10.2019

suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(oncoscanR))

#' Retrieve arguments from command line.
#' Expects that the first argument is the ChAS file and the second the gender (F or M).
args <- commandArgs(TRUE)
if(length(args) != 2){
  stop("The first argument has to be the ChAS file and the second the gender (F or M).")
}

chas.fn <- args[1]
gender <- args[2]

dat <- workflow_oncoscan.run(chas.fn, gender)
print(toJSON(dat, auto_unbox=TRUE, pretty=TRUE))

