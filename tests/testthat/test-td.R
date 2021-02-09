test_that("TD score works - only gains, border checks", {
  cov <- GRanges(seqnames = factor(c(paste0(1:4, 'p'), paste0(1:4, 'q')),
                                   levels = c(paste0(1:4, 'p'), paste0(1:4, 'q'))),
                 ranges = IRanges(start = rep(1, 8),
                                  end = rep(200*10^6, 8)))
  segs <- GRanges(seqnames = factor(c(rep('1p', 5), rep('2p', 4)), levels = paste0(1:4, 'p')),
                  ranges = IRanges(start = c(1,   2,   3,  5, 20.000001, 50, 90.000001, 95, 100)*10^6 + 1,
                                   end =   c(1.1, 2.3, 4, 15, 30,        80, 91, 96.000001, 110.000001)*10^6),
                  cn = c(3, 4, 3, 4, 3, 4, 3, 4, 3),
                  cn.type = rep(cntype.gain, 9),
                  cn.subtype = rep(cntype.gain, 9))
  n <- score_td(segs)
  expected <- list(TDplus=3, TD=4)
  expect_equal(n, expected)
})

test_that("TD score works - CN checks", {
  cov <- GRanges(seqnames = factor(c(paste0(1:4, 'p'), paste0(1:4, 'q')),
                                   levels = c(paste0(1:4, 'p'), paste0(1:4, 'q'))),
                 ranges = IRanges(start = rep(1, 8),
                                  end = rep(200*10^6, 8)))
  segs <- GRanges(seqnames = factor(c(rep('1p', 7), rep('2p', 7)), levels = paste0(1:4, 'p')),
                  ranges = IRanges(start = rep(c(1,   10),7)*10^6 + 1,
                                   end =   rep(c(1.5, 15),7)*10^6),
                  cn = c(0,0,1,1,NA,NA,3,3,4,4,5,5,10,10),
                  cn.type = c(cntype.loss, cntype.loss,
                              cntype.loss, cntype.loss,
                              cntype.loh, cntype.loh,
                              cntype.gain, cntype.gain,
                              cntype.gain, cntype.gain,
                              cntype.gain, cntype.gain,
                              cntype.gain, cntype.gain),
                  cn.subtype = c(cntype.homloss, cntype.homloss,
                                 cntype.hetloss, cntype.hetloss,
                                 cntype.loh, cntype.loh,
                                 cntype.gain, cntype.gain,
                                 cntype.gain, cntype.gain,
                                 cntype.weakamp, cntype.weakamp,
                                 cntype.strongamp, cntype.strongamp))
  n <- score_td(segs)
  expected <- list(TDplus=2, TD=2)
  expect_equal(n, expected)
})

test_that("TD scores work - real case", {
  chas.fn <- system.file("testdata", "TDplus_gene_list_full_location.txt", package = "oncoscanR")
  segments <- load_chas(chas.fn, oncoscanR::oncoscan_na33.cov)
  segments$cn.subtype <- get_cn_subtype(segments, 'F')
  segs.clean <- trim_to_coverage(segments, oncoscanR::oncoscan_na33.cov) %>%
    prune_by_size()

  n <- score_td(segs.clean)
  expect_true(n$TDplus==102 && n$TD==21) #Verified by hand in Excel and ChAS
})

test_that("TD scores work - empty segments", {
  segs <- GRanges(seqnames = factor(c(), levels = paste0(1:4, 'p')),
                  ranges = IRanges())
  n <- score_td(segs)
  expect_identical(unlist(n), c(TDplus=0,TD=0))
})
