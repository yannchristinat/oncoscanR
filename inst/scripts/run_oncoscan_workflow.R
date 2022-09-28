#!/usr/bin/Rscript
#
# Run the standard workflow for the detection of arm-level alterations and
# scores calculation and prints the output in JSON.
#
# Author: Yann Christinat
# Date: 9.8.2021

suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(oncoscanR))

#' Retrieve arguments from command line.
#' Expects that the first argument is the ChAS file.
args <- commandArgs(TRUE)
if(length(args) != 1){
  stop("The first argument has to be the ChAS file ")
}

chas.fn <- args[1]

dat <- workflow_oncoscan.chas(chas.fn)
print(toJSON(dat, auto_unbox=TRUE, pretty=TRUE))

