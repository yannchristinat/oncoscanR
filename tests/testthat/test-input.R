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
  segs.filename <- system.file("extdata", "chas_example-withDuplicates.txt",
                               package = "oncoscanR")
  expect_error(load_chas(segs.filename, kit.coverage = oncoscan_na33.cov),
              'contains duplicated entries.')
})


test_that("loading large ChAS file", {
  segs.filename <- system.file("extdata", "LST_gene_list_full_location.txt",
                               package = "oncoscanR")
  segs <- load_chas(segs.filename, kit.coverage = oncoscan_na33.cov)
  expect_true(length(segs)>100) # Mostly testing that it does not fail nor return zero segments
})


test_that("loading ChAS annotation file works", {
  filename <- system.file("extdata", "OncoScan.na33.r1.annot.csv.zip",
                               package = "oncoscanR")
  cov <- get_oncoscan_coverage_from_probes(filename)

  found <- 0
  for(i in seq_along(oncoscan_na33.cov)){
    segA <- oncoscan_na33.cov[i]
    for(j in seq_along(cov)){
      segB <- cov[j]
      if(same_segments(segA, segB)){
        found <- found+1
      }
    }
  }

  expect_true(length(cov) == 44 & found == 44)
})

