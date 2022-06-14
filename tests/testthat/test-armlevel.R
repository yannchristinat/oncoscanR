#' Test with no segments, 1 segment >80%, 1 segment <80%, 1 segment at 80%, 1 segment at 100%,
#' one arm with 3 segments and >80% in total, one arm with 3 segments and <80% in total.
test_that("Arm-level alteration detection works - basic", {
  cov <- GRanges(seqnames = factor(paste0(1:7, 'p'), levels = paste0(1:7, 'p')),
                 ranges = IRanges(start = rep(1, 7),
                                  end = rep(100, 7)))
  segs <- GRanges(seqnames = factor(c('2p','3p','4p','5p','6p','6p','6p','7p','7p','7p'),
                                    levels = paste0(1:7, 'p')),
                  ranges = IRanges(start = c(6,  21,  1,   1,  1, 21,  61,  1, 31, 91),
                                   end   = c(95, 80, 80, 100, 20, 50, 100, 10, 40, 100)),
                  cn = c(0,1,NA,3,5,10,0,1,NA,3),
                  cn.type = c("Loss", "Loss", "LOH", "Gain",
                              "Gain", "Gain", "Loss", "Loss",
                              "LOH", "Gain"))
  arms <- armlevel_alt(segs, cov, threshold = 0.8)
  expected <- c('2p'=0.9, '4p'=0.8, '5p'=1, '6p'=0.9)
  expect_equal(arms, expected)
})


test_that("Arm-level alteration detection works - 0% threshold", {
  cov <- GRanges(seqnames = factor(paste0(1:7, 'p'), levels = paste0(1:7, 'p')),
                 ranges = IRanges(start = rep(1, 7),
                                  end = rep(100, 7)))
  segs <- GRanges(seqnames = factor(c('2p','3p','4p','5p','6p','6p','6p','7p','7p','7p'),
                                    levels = paste0(1:7, 'p')),
                  ranges = IRanges(start = c(6,  21,  1,   1,  1, 21,  61,  1, 31, 91),
                                   end   = c(95, 80, 80, 100, 20, 50, 100, 10, 40, 100)),
                  cn = c(0,1,NA,3,5,10,0,1,NA,3),
                  cn.type = c("Loss", "Loss", "LOH", "Gain",
                              "Gain", "Gain", "Loss", "Loss",
                              "LOH", "Gain"))
  arms <- armlevel_alt(segs, cov, 0)
  expected <- c('1p'=0, '2p'=0.9, '3p'=0.6, '4p'=0.8, '5p'=1, '6p'=0.9, '7p'=0.3)
  expect_equal(arms, expected)
})


test_that("Arm-level alteration detection works - 100% threshold", {
  cov <- GRanges(seqnames = factor(paste0(1:7, 'p'), levels = paste0(1:7, 'p')),
                 ranges = IRanges(start = rep(1, 7),
                                  end = rep(100, 7)))
  segs <- GRanges(seqnames = factor(c('2p','3p','4p','5p','6p','6p','6p','7p','7p','7p'),
                                    levels = paste0(1:7, 'p')),
                  ranges = IRanges(start = c(6,  21,  1,   1,  1, 21,  61,  1, 31, 91),
                                   end   = c(95, 80, 80, 100, 20, 50, 100, 10, 40, 100)),
                  cn = c(0,1,NA,3,5,10,0,1,NA,3),
                  cn.type = c("Loss", "Loss", "LOH", "Gain",
                              "Gain", "Gain", "Loss", "Loss",
                              "LOH", "Gain"))
  arms <- armlevel_alt(segs, cov, 1)
  expected <- c('5p'=1)
  expect_equal(arms, expected)
})


test_that("Arm-level alteration detection works - overlaps", {
  cov <- GRanges(seqnames = factor(paste0(1:4, 'p'), levels = paste0(1:4, 'p')),
                 ranges = IRanges(start = c(1,1,101,101),
                                  end = c(100,100,200,200)))
  segs <- GRanges(seqnames = factor(c(rep('1p',3), rep('2p',3), rep('3p',3), rep('4p',3)),
                                    levels = paste0(1:4, 'p')),
                  ranges = IRanges(start = c( 1, 16, 11,  1, 16, 11, 51,  99, 181, 51,  99, 171),
                                   end   = c(20, 50, 70, 20, 50, 80, 90, 110, 210, 90, 150, 210)),
                  cn = rep(c(0,6,NA), 4),
                  cn.type = rep(c("Loss", "Gain", "LOH"), 4))
  arms <- armlevel_alt(segs, cov, 0)
  expected <- c('1p'=0.7, '2p'=0.8, '3p'=0.3, '4p'=0.8)
  expect_equal(arms, expected)
})

test_that("Arm-level alteration detection works - empty segments", {
  cov <- GRanges(seqnames = factor(paste0(1:4, 'p'), levels = paste0(1:4, 'p')),
                 ranges = IRanges(start = c(1,1,101,101),
                                  end = c(100,100,200,200)))
  segs <- GRanges(seqnames = factor(c(), levels = paste0(1:4, 'p')),
                  ranges = IRanges())
  arms <- armlevel_alt(segs, cov)
  expect_equal(length(arms), 0)
})

