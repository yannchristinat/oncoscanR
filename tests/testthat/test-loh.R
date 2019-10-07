test_that("HR-LOH score works - simple cases", {
  cov <- GRanges(seqnames = factor(c(paste0(1:4, 'p'), paste0(1:4, 'q')),
                                   levels = c(paste0(1:4, 'p'), paste0(1:4, 'q'))),
                 ranges = IRanges(start = rep(1, 8),
                                  end = rep(200*10^6, 8)))
  segs <- GRanges(seqnames = factor(c(rep('1p', 6), rep('2p', 2)), levels = paste0(1:4, 'p')),
                  ranges = IRanges(start = c( 1, 10, 40, 60.000001, 80, 100,  1, 50)*10^6 + 1,
                                   end =   c(10, 30, 55, 75,        95, 120, 20, 80)*10^6),
                  cn = c(NA, NA, NA, NA, 3, 1, NA, 15),
                  cn.type = c(cntype.loh, cntype.loh, cntype.loh, cntype.loh,
                              cntype.gain, cntype.loss, cntype.loh, cntype.gain),
                  cn.subtype = c(cntype.loh, cntype.loh, cntype.loh, cntype.loh,
                                 cntype.gain, cntype.hetloss, cntype.loh, cntype.strongamp))
  n <- score_loh(segs, cov, c('2p', '2q'))
  expect_equal(n, 3)
})


test_that("HR-LOH score works - overlapping loh", {
  cov <- GRanges(seqnames = factor(c(paste0(1:4, 'p'), paste0(1:4, 'q')),
                                   levels = c(paste0(1:4, 'p'), paste0(1:4, 'q'))),
                 ranges = IRanges(start = rep(1, 8),
                                  end = rep(200*10^6, 8)))
  segs <- GRanges(seqnames = factor(c(rep('1p', 6), rep('2p', 2)), levels = paste0(1:4, 'p')),
                  ranges = IRanges(start = c( 1, 10, 40, 60.000001, 80, 100,  1, 50)*10^6 + 1,
                                   end =   c(10, 20, 55, 75,        95, 120, 20, 80)*10^6),
                  cn = c(NA, NA, NA, NA, 3, 1, NA, 15),
                  cn.type = c(cntype.loh, cntype.loh, cntype.loh, cntype.loh,
                              cntype.gain, cntype.loss, cntype.loh, cntype.gain),
                  cn.subtype = c(cntype.loh, cntype.loh, cntype.loh, cntype.loh,
                                 cntype.gain, cntype.hetloss, cntype.loh, cntype.strongamp))
  n <- score_loh(segs, cov, c())
  expect_equal(n, 4)
})
