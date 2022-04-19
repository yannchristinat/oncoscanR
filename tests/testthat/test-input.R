test_that("loading ChAS file works", {
  segs.filename <- system.file("extdata", "chas_example.txt",
                               package = "oncoscanR")
  segs <- load_chas(segs.filename, kit.coverage = oncoscan_na33.cov)

  found <- 0
  for(i in seq_along(segs.chas_example)){
    segA <- segs.chas_example[i]
    for(j in seq_along(segs)){
      segB <- segs[j]
      if(same_segments(segA, segB)){
        found <- found+1
      }
    }
  }

  expect_true(length(segs) == 15 & found == 15) # 14 segments in file but one is covering two arms
})


test_that("loading ChAS file with duplicated rows fails", {
  segs.filename <- system.file("testdata", "chas_example-withDuplicates.txt",
                               package = "oncoscanR")
  expect_error(load_chas(segs.filename, kit.coverage = oncoscan_na33.cov),
              'contains duplicated entries.')
})


test_that("loading large ChAS file", {
  skip_on_cran()

  segs.filename <- system.file("testdata", "LST_gene_list_full_location.txt",
                               package = "oncoscanR")
  segs <- load_chas(segs.filename, kit.coverage = oncoscan_na33.cov)
  expect_true(length(segs)>100) # Mostly testing that it does not fail nor return zero segments
})


test_that("loading ChAS annotation file works", {
  filename <- system.file("extdata", "OncoScan.na33.r1.annot.csv.chr20.zip",
                               package = "oncoscanR")
  cov <- get_oncoscan_coverage_from_probes(filename)

  found <- 0
  chr20 <- oncoscan_na33.cov[as.vector(seqnames(oncoscan_na33.cov)) %in% c('20p','20q')]
  for(i in seq_along(chr20)){
    segA <- chr20[i]
    for(j in seq_along(cov)){
      segB <- cov[j]
      if(same_segments(segA, segB)){
        found <- found+1
      }
    }
  }

  expect_true(length(cov) == 2 & found == 2)
})

test_that("Loading ChAS with missing 'Full Location' column fails",{
  segs.filename <- system.file("testdata", "chas_example-noFullLocation.txt",
                               package = "oncoscanR")
  expect_error(suppressWarnings(load_chas(segs.filename, kit.coverage = oncoscan_na33.cov)),
               'Parsing ChAS file failed')
})

test_that("loading ChAS with 'chr' name scheme works", {
  segs.filename <- system.file("testdata", "chas_example-withChr.txt",
                               package = "oncoscanR")
  segs <- load_chas(segs.filename, kit.coverage = oncoscan_na33.cov)

  found <- 0
  for(i in seq_along(segs.chas_example)){
    segA <- segs.chas_example[i]
    for(j in seq_along(segs)){
      segB <- segs[j]
      if(same_segments(segA, segB)){
        found <- found+1
      }
    }
  }

  expect_true(length(segs) == 15 & found == 15) # 14 segments in file but one is covering two arms
})

test_that("Loading ChAS with empty file works",{
  segs.filename <- system.file("testdata", "chas_example-empty.txt",
                               package = "oncoscanR")
  expect_warning(segs <- load_chas(segs.filename, kit.coverage = oncoscan_na33.cov),
                 'No segments loaded')
  expect_equal(length(segs), 0)
})

