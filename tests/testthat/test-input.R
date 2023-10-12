test_that("loading ChAS file works", {
  segs.filename <- system.file("extdata", "chas_example.txt",
                               package = "oncoscanR")
  segs <- load_chas(segs.filename, kit.coverage = oncoscan_na33.cov)

  split.segs <- GRanges(seqnames = c('7p','7q'),
                        ranges=IRanges(start=c(55219200, 61545024),
                                       end = c(58025423, 116328000)))
  split.segs$cn <- c(7, 7)
  split.segs$cn.type <- c('Gain', 'Gain')


  found <- 0
  for(i in seq_along(split.segs)){
    segA <- split.segs[i]
    for(j in seq_along(segs)){
      segB <- segs[j]
      if(same_segments(segA, segB)){
        found <- found+1
      }
    }
  }

  expect_equal(length(segs), 15)  
  # 14 segments in file but one is covering two arms
  expect_equal(found, 2, info = print(found))

})

test_that("loading ChAS file with separators in Full Location column works", {
    segs.filename <- "../testdata/chas_example-withSeparators.txt"
    segs <- load_chas(segs.filename, kit.coverage = oncoscan_na33.cov)
    
    split.segs <- GRanges(seqnames = c('7p','7q'),
                          ranges=IRanges(start=c(55219200, 61545024),
                                         end = c(58025423, 116328000)))
    split.segs$cn <- c(7, 7)
    split.segs$cn.type <- c('Gain', 'Gain')
    
    
    found <- 0
    for(i in seq_along(split.segs)){
        segA <- split.segs[i]
        for(j in seq_along(segs)){
            segB <- segs[j]
            if(same_segments(segA, segB)){
                found <- found+1
            }
        }
    }
    
    expect_equal(length(segs), 15)  
    # 14 segments in file but one is covering two arms
    expect_equal(found, 2, info = print(found))
    
})

test_that("loading ChAS file with duplicated rows fails", {
  segs.filename <- "../testdata/chas_example-withDuplicates.txt"
  expect_error(load_chas(segs.filename, kit.coverage = oncoscan_na33.cov),
              'contains duplicated entries.')
})


test_that("loading large ChAS file", {
  segs.filename <- system.file("extdata", "LST_gene_list_full_location.txt",
                               package = "oncoscanR")
  segs <- load_chas(segs.filename, kit.coverage = oncoscan_na33.cov)
  expect_true(length(segs)>100) 
  # Mostly testing that it does not fail nor return zero segments
})


test_that("Loading ChAS with missing 'Full Location' column fails",{
  segs.filename <- "../testdata/chas_example-noFullLocation.txt"
  expect_error(suppressWarnings(load_chas(segs.filename, 
                                          kit.coverage = oncoscan_na33.cov)),
               'Parsing ChAS file failed')
})

test_that("loading ChAS with 'chr' name scheme works", {
  segs.filename <- "../testdata/chas_example-withChr.txt"
  segs <- load_chas(segs.filename, kit.coverage = oncoscan_na33.cov)

  split.segs <- GRanges(seqnames = c('7p','7q'),
                        ranges=IRanges(start=c(55219200, 61545024),
                                       end = c(58025423, 116328000)))
  split.segs$cn <- c(7, 7)
  split.segs$cn.type <- c('Gain', 'Gain')

  found <- 0
  for(i in seq_along(split.segs)){
    segA <- split.segs[i]
    for(j in seq_along(segs)){
      segB <- segs[j]
      if(same_segments(segA, segB)){
        found <- found+1
      }
    }
  }

  expect_equal(c(length(segs), found), c(15, 2)) 
  # 14 segments in file but one is covering two arms
})

test_that("Loading ChAS with empty file works",{
  segs.filename <- "../testdata/chas_example-empty.txt"
  expect_warning(segs <- load_chas(segs.filename, 
                                   kit.coverage = oncoscan_na33.cov),
                 'No segments loaded')
  expect_equal(length(segs), 0)
})

test_that("Loading ASCAT file works", {
    segs.filename <- system.file("extdata", "ascat_example.txt",
                                 package = "oncoscanR")
    segs <- load_ascat(segs.filename, kit.coverage = oncoscan_na33.cov)
    
    expected.segs <- GRanges(seqnames = c('1q','5p'),
                          ranges=IRanges(start=c(144009053, 38139),
                                         end = c(249212878, 46193462)))
    expected.segs$cn <- c(3, 3)
    expected.segs$cn.type <- c('Gain', 'Gain')
    
    
    found <- 0
    for(i in seq_along(expected.segs)){
        segA <- expected.segs[i]
        for(j in seq_along(segs)){
            segB <- segs[j]
            if(same_segments(segA, segB)){
                found <- found+1
            }
        }
    }
    
    expect_equal(length(segs), 2)  # 2 segments in file 
    expect_equal(found, 2, info = print(found))
})

test_that("Loading ASCAT with empty file works",{
    segs.filename <- "../testdata/ascat_example-empty.txt"
    expect_warning(segs <- load_ascat(segs.filename, 
                                      kit.coverage = oncoscan_na33.cov),
                   'No segments loaded')
    expect_equal(length(segs), 0)
})

