# input.R
#
# Functions to handle loading/cleaning of segments and reference files.
#
# Author: Yann Christinat
# Date: 19.9.2019


#' Load a ChAS text export file.
#'
#' @details The ChAS file is expected to have the following column names: "CN State" (number or empty), "Type" (expected
#' value: "Gain", "Loss" or "LOH") and "Full Location" (in the format "chr:start-end").
#'
#' @param filename Path to the ChAS file.
#' @param kit.coverage A \code{GRanges} object containing the regions covered on each chromosome arm by the kit.
#'
#' @return A \code{GRanges} object containing the segments, their copy number (field \code{cn}), their copy
#' number types (field \code{cn.type}). \code{cn.type} contains either "Gain",
#' "Loss" or "LOH".
#' If the file contains twice the same segment or does not respect the format specifications, then an error is
#' raised. NB. The chromosome name is in the format "1" and not "chr1" and will be transformed if needed.
#'
#' @export
#'
#' @importFrom readr read_tsv cols_only col_number col_character
#' @import GenomicRanges
#' @import IRanges
#'
#' @examples
#' segs.filename <- system.file("extdata", "chas_example.txt", package = "oncoscanR")
#' segs.chas_example <- load_chas(segs.filename, oncoscan_na33.cov)
load_chas <- function(filename, kit.coverage){
  # Reads in the ChAS file
  oncoscan_table <- read_tsv(filename, comment = "#", col_names = TRUE,
                             col_types = cols_only(`CN State` = col_number(),
                                                   Type = col_character(),
                                                   `Full Location` = col_character()))

  if(dim(oncoscan_table)[2]!=3){
    stop('Parsing ChAS file failed.')
  }
  if(dim(oncoscan_table)[1]==0){
    warning('No segments loaded!')
    return(GRanges())
  }


  # Allocate the place for the GRanges segments. At most we will end up with twice as much segments
  # as present in the raw data.
  segments_list <- vector(mode = "list", length = 2*dim(oncoscan_table)[1])

  # Parse throuh all lines of the data
  counter <- 0
  for (i in 1:dim(oncoscan_table)[1]){
    counter <- counter+1
    # Start: Extract chr no, start, end, copy number

    # Full location from oncoscan file
    loc_cord <- oncoscan_table$`Full Location`[i]
    loc_cord_list <- strsplit(loc_cord, split=':', fixed=TRUE)[[1]]

    # seg chr no based on oncoscan file
    seg_chr <- loc_cord_list[1]
    seg_chr <- gsub("chr","",seg_chr) #remove 'chr' if present

    # seg coordinates
    coord_list <- strsplit(loc_cord_list[2],split='-', fixed=TRUE)[[1]]
    seg_start <- as.numeric(coord_list[1])
    seg_end <- as.numeric(coord_list[2])

    seg_cn <- oncoscan_table$`CN State`[i]
    seg_cntype <- oncoscan_table$Type[i]

    # Test if copy number type is correct
    if(!(seg_cntype %in% c(cntype.gain, cntype.loss, cntype.loh))){
      stop(paste('The column "Type" should contain only the following values:',
                 c(cntype.gain, cntype.loss, cntype.loh), 'whereas',seg_cntype, 'was found!'))
    }

    # Get the arms
    parm <- kit.coverage[seqnames(kit.coverage) == paste0(seg_chr, 'p')]
    qarm <- kit.coverage[seqnames(kit.coverage) == paste0(seg_chr, 'q')]

    if(length(parm)>0 || length(qarm)>0){
      if(length(parm)==0 || seg_start > end(parm)){ # Then the segment is only in the q arm
        seg <- GRanges(seqnames = factor(paste0(seg_chr, 'q'),
                                         levels = levels(seqnames(kit.coverage))),
                       ranges = IRanges(start = seg_start, end = seg_end),
                       cn = seg_cn, cn.type = seg_cntype)
        segments_list[[counter]] <- seg
      }
      else if (length(qarm)==0 || seg_end < start(qarm)){ # Then the segment is only in the p arm
        seg <- GRanges(seqnames = factor(paste0(seg_chr, 'p'),
                                         levels = levels(seqnames(kit.coverage))),
                       ranges = IRanges(start = seg_start, end = seg_end),
                       cn = seg_cn, cn.type = seg_cntype)
        segments_list[[counter]] <- seg
      }
      else {
        # Create a segment for each arm
        seg <- GRanges(seqnames = factor(paste0(seg_chr, c('p','q')),
                                         levels = levels(seqnames(kit.coverage))),
                       ranges = IRanges(start = rep(seg_start, 2),
                                        end = rep(seg_end, 2)))

        # Get the segments overlapping with the arms
        o <- IRanges::findOverlapPairs(c(parm, qarm), seg)
        new_segs <- pintersect(o)
        new_segs$cn <- rep(seg_cn, length(new_segs))
        new_segs$cn.type <- rep(seg_cntype, length(new_segs))
        new_segs$hit <- NULL

        segments_list[[counter]] <- new_segs
      }
    }
  }

  segs <- do.call("c", unlist(segments_list))

  # Test for duplicated entries
  sapply(unique(seqnames(segs)), function(arm){
    arm_segs <- sort(segs[seqnames(segs) == arm])
    if(length(arm_segs)>1){
      for(i in 2:length(arm_segs)){
        segA <- arm_segs[i-1]
        segB <- arm_segs[i]
        if(start(segA) == start(segB) & end(segA) == end(segB) & segA$cn.type == segB$cn.type){
          if(segA$cn.type == cntype.loh){
            stop(paste('The file', filename, 'contains duplicated entries.'))
          }
          else if(segA$cn == segB$cn){
            stop(paste('The file', filename, 'contains duplicated entries.'))
          }
        }
      }
    }
  })

  if(length(segs)==0){
    warning('No segments loaded!')
  }

  return(segs)
}


#' Load a CNV file in a BED-like format.
#'
#' @details The file is expected to contain the following columns: chromosome name, segment start, segment end,
#' copy number. An fifth column with the copy number type is optional and if present should contain the values
#' "Gain", "Loss" or "LOH".
#' If the gender is not specified, then the segments on the sexual chromosomes (X and Y) are dropped as the
#' copy number type cannot be determined for those.
#'
#' @param filename Path to the CNV file.
#' @param gender Character to indicate whether the sample is male ("M") or female ("F").
#'
#' @return A \code{GRanges} object containing the segments, their copy number and copy number types. If the
#' file contains twice the same segment or does not respect the format specifications, then an error is raised.
#'
#' @export
#'
#' @examples
#' segs.filename <- system.file("extdata", "cnv_example.bed", package = "oncoscanR")
#' segs.cnv_example <- load_bed(segs.filename, "F")
load_bed <- function(filename, gender){
  warning("Not yet implemented! Please request it if needed.")
  return(NULL)
}


#' Load the oncoscan annotation file and infer the covered regions from it.
#'
#' @details Expects the following columns from the annotation file (as the first line of non-comment):
#'   - "Chromosome": chromosome name
#'   - "Physical Position": genomic position of the SNP
#'   - "Cytoband": cytoband name (e.g. 3q22.1)
#'
#' @param filename Path to the ChAS annotation file (either compressed or not). Defaults to the Oncoscan na33.r1
#' file present in the package.
#'
#' @return A \code{GRanges} object containing the regions covered on each chromosome arm. If the file does
#' not respect the format specifications, then an error is raised.
#'
#' @export
#'
#' @import GenomicRanges
#' @import IRanges
#' @importFrom readr read_csv cols_only col_character col_integer
#'
#' @examples
#' oncoscan_na33.cov <- get_oncoscan_coverage_from_probes()
get_oncoscan_coverage_from_probes <- function(filename = system.file("extdata",
                                                                     "OncoScan.na33.r1.annot.csv.zip",
                                                                     package = "oncoscanR")){
  # Read the annotation file
  dat <- read_csv(filename, comment = "#", col_names = TRUE,
               col_types = cols_only(Chromosome = col_character(),
                                     `Physical Position` = col_integer(),
                                     Cytoband = col_character()))

  # Set a new column containing the chromosomal arm name (e.g. 3p). Remove 'chr' if present
  dat$Arm <- factor(paste0(gsub("chr", "", dat$Chromosome), substr(dat$Cytoband,1,1)))

  # Get the min and max position of probes within each arm
  segs <- lapply(levels(dat$Arm), function(a){
    minpos <- min(dat$`Physical Position`[dat$Arm == a])
    maxpos <- max(dat$`Physical Position`[dat$Arm == a])

    seg <- GRanges(seqnames = factor(a, levels = levels(dat$Arm)), ranges = IRanges(start = minpos, end = maxpos))
    return(seg)
  })

  # Return all arms
  return(do.call("c", segs))
}


#' Get the detailed copy number types.
#'
#' @details The copy number subtypes are defined as follow:
#'   - "Gain": 1-2 extra copies
#'   - "Weak amplification": 3-7 extra copies
#'   - "Strong amplification": 8 or more extra copies
#'   - "Heterozygote loss": Loss of one copy out of two
#'   - "Homozygote loss": Loss of all copies
#'   - "LOH": copy-neutral loss of one parental allele
#'
#' @param segments A \code{GRanges} object containing the segments, their copy number (field \code{cn}) and
#' their copy number types (field \code{cn.type}). The \code{cn.type} is expected to be either "Gain", "Loss"
#' or "LOH".
#' @param gender Character to indicate whether the sample is male ("M") or female ("F"). If \code{NULL}
#' then the field \code{cn.type} is set to \code{NA} for the chromosomes X and Y.
#'
#' @return A list of copy number types (one for each segment). Raises an error if the \code{cn.type},
#' \code{cn} and gender are not compatible.
#'
#' @import GenomicRanges
#'
#' @export
#'
#' @examples
#' subtypes <- get_cn_subtype(segs.chas_example, 'F')
get_cn_subtype <- function(segments, gender){
  is.cn_segment(segments, raise_error = TRUE)

  if(!(gender %in% c('M', 'F'))){
    message("Unspecified gender. Cannot assign subtypes on sexual chromosomes.")
  }
  sexchroms <- c('Xp','Xq','Yp','Yq')

  subtypes <- lapply(seq_along(segments), function(i){
    seg <- segments[i]
    if(length(intersect(seqnames(seg), c(sexchroms, paste0('chr',sexchroms))))>0 & !(gender %in% c('M', 'F'))){
      return(NA)
    }
    if(seg$cn.type == 'LOH'){
      return(cntype.loh)
    }
    s <- NA
    if(length(intersect(seqnames(seg), c(sexchroms, paste0('chr',sexchroms))))>0 & gender == 'M'){
      s <- ifelse(seg$cn < 1,
                  cntype.homloss,
                  ifelse(seg$cn > 1,
                         ifelse(seg$cn > 3,
                                ifelse(seg$cn > 8,
                                       cntype.strongamp,
                                       cntype.weakamp),
                                cntype.gain),
                         NA))
    }
    else {
      s <- ifelse(seg$cn < 2,
                  ifelse(seg$cn < 1,
                         cntype.homloss,
                         cntype.hetloss),
                  ifelse(seg$cn > 2,
                         ifelse(seg$cn > 4,
                                ifelse(seg$cn > 9,
                                       cntype.strongamp,
                                       cntype.weakamp),
                                cntype.gain),
                         NA))
    }
    return(s)
  })
  return(subtypes)
}


#' Trim segments with respect to the kit's coverage.
#'
#' @details All segments that are not entirely contained within the kit coverage will be trimmed to the coverage's
#' limits.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy number and copy number types.
#' @param kit.coverage A \code{GRanges} object containing the regions covered on each chromosome arm.
#'
#' @return A \code{GRanges} object containing the cleaned segments, their copy number and copy number types.
#' @export
#'
#' @import GenomicRanges
#' @import IRanges
#'
#' @examples
#' segs.trimmed <- trim_to_coverage(segs.chas_example, oncoscan_na33.cov)
trim_to_coverage <- function(segments, kit.coverage){
  is.cn_segment(segments, raise_error = TRUE)

  if(length(segments)==0){
    return(segments)
  }

  # Apply on each arm...
  segs.clean <- lapply(unique(seqnames(segments)), function(arm){
    # Trim segments wrt coverage
    arm.cov <- kit.coverage[seqnames(kit.coverage) == arm]
    o <- IRanges::findOverlapPairs(segments[seqnames(segments) == arm], arm.cov)
    new_segs <- sort(pintersect(o))
    if(length(new_segs) == 0){
      return(new_segs)
    }
    new_segs$hit <- NULL
    return(new_segs)
  })
  return(do.call("c", unlist(segs.clean)))
}


#' Merge segments with respect to the kit resolution and the copy number.
#'
#' @details If two segments are at a distance smaller than the resolution, then the segments are merged if the
#' share the same \code{cn} value. Note that the function does not look at the copy number type or subtype but
#' only at the actual copy number to decide whether segments can be merged.
#'
#'
#' @param segments A \code{GRanges} object containing the segments, their copy number and copy number types.
#' @param kit.resolution Number >0 indicating the minimum segment size detectable by the technique (in kilobases).
#' Defaults to the Oncoscan assay resolution outside of cancer genes: 300Kb.
#'
#' @return A \code{GRanges} object containing the cleaned segments, their copy number and copy number types.
#' @export
#'
#' @import GenomicRanges
#'
#' @examples
#' segs.merged <- merge_segments(segs.chas_example)
#' segs.merged_50k <- merge_segments(segs.chas_example, 50)
merge_segments <- function(segments, kit.resolution = 300){
  is.cn_segment(segments, raise_error = TRUE)

  if(kit.resolution<1/1000){
    stop("Kit resolution has to be greater than zero.")
  }

  if(length(segments)==0){
    return(segments)
  }

  # Go through each arm...
  segs.merged <- lapply(unique(seqnames(segments)), function(arm){
    # Go through each copy number
    armsegs <- segments[seqnames(segments) == arm]
    armsegs.merged <- lapply(unique(armsegs$cn), function(cn){
      segs <- NULL
      if(is.na(cn)){
        segs <- armsegs[is.na(armsegs$cn)]
      }
      else {
        segs <- armsegs[!is.na(armsegs$cn) & armsegs$cn == cn]
      }
      new_segs <- reduce(segs, min.gapwidth = kit.resolution*1000)
      # Re-set the copy number, cn type and subtype
      new_segs$cn <- cn
      new_segs$cn.type <- segs[1]$cn.type
      if(!is.null(segs$cn.subtype)){
        new_segs$cn.subtype <- segs[1]$cn.subtype
      }
      return(new_segs)
    })
    do.call("c", unlist(armsegs.merged))
  })
  return(do.call("c", unlist(segs.merged)))
}


#' Trim LOH segments with respect to loss segments.
#'
#' @details LOH segments completely contained within (or equal to)  a copy loss segment are deleted.
#' LOH segments partially overlapping (on one end only) with a copy loss segment are trimmed to
#' remove the overlap. If a copy loss segment is completely contained within (but not equal to) a LOH
#' segment, then nothing is done; the overlap remains.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy number and copy number types.
#'
#' @return A \code{GRanges} object containing the cleaned segments, their copy number and copy number types.
#' @export
#'
#' @import GenomicRanges
#' @importFrom IRanges findOverlapPairs
#' @importFrom S4Vectors first second
#'
#' @examples
#' segs.adj <- adjust_loh(segs.chas_example)
adjust_loh <- function(segments){
  is.cn_segment(segments, raise_error = TRUE)

  if(length(segments)==0){
    return(segments)
  }

  # Apply on each arm...
  loh.adj <- lapply(unique(seqnames(segments)), function(arm){
    # Detect overlaps between LOH and Loss and trim the LOH segment if necessary
    segs.loh <- segments[segments$cn.type == cntype.loh & seqnames(segments) == arm]
    segs.loss <- segments[segments$cn.type == cntype.loss & seqnames(segments) == arm]
    op <- IRanges::findOverlapPairs(segs.loh, segs.loss)

    loh.todelete <- GRanges()
    loh.toadd <- GRanges()

    for(i in seq_along(op)){
      loh <- S4Vectors::first(op[i])
      loss <- S4Vectors::second(op[i])
      if(loh %within% loss){ # If the LOH is completely within the loss, then delete LOH segment
        loh.todelete <- append(loh.todelete, loh)
      }
      else if (!(start(loh)<start(loss) && end(loh)>end(loss))){
        # If the LOH is not larger than the loss on both ends, then
        new_loh <- setdiff(loh, loss)
        #new_loh$cn <- NA
        #new_loh$cn.type <- cntype.loh
        #new_loh$cn.subtype <- cntype.loh


        loh.toadd <- append(loh.toadd, new_loh)
        loh.todelete <- append(loh.todelete, loh)
      }
    }

    adj <- union(setdiff(segs.loh, loh.todelete), loh.toadd)
    if(length(adj)>0){
      adj$cn <- NA
      adj$cn.type <- cntype.loh
      adj$cn.subtype <- cntype.loh
    }

    return(adj)
  })

  return(c(segments[segments$cn.type != cntype.loh], do.call("c", unlist(loh.adj))))
}


#' Remove segments smaller than the kit resolution.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy number and copy number types.
#' @param threshold Number indicating the minimum segment size to be kept (in kilobases).
#' Defaults to the Oncoscan assay resolution outside of cancer genes: 300Kb.
#'
#' @return A \code{GRanges} object containing the cleaned segments, their copy number and copy number types.
#' @export
#'
#' @importFrom IRanges width
#'
#' @examples
#' segs.300k <- prune_by_size(segs.chas_example)
#' segs.50k <- prune_by_size(segs.chas_example, 50)
prune_by_size <- function(segments, threshold = 300){
  is.cn_segment(segments, raise_error = TRUE)

  if(threshold<0){
    stop("Threshold has to be greater than or equal to zero.")
  }
  if(length(segments)==0){
    return(segments)
  }

  return(segments[width(segments) >= threshold*1000])
}
