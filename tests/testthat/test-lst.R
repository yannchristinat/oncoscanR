test_that("LST works - simple cases", {
  cov <- GRanges(seqnames = factor(paste0(1:4, 'p'), levels = paste0(1:4, 'p')),
                 ranges = IRanges(start = rep(1, 4),
                                  end = rep(100*10^6, 4)))
  segs <- GRanges(seqnames = factor(c(rep('1p', 4), rep('2p', 2)), levels = paste0(1:4, 'p')),
                  ranges = IRanges(start = c(30, 50, 70,  90, 25,  50)*10^6 + 1,
                                   end =   c(50, 70, 90, 100, 75, 100)*10^6),
                  cn = c(NA, 3, 1, 0, NA, 14),
                  cn.type = c(cntype.loh, cntype.gain, cntype.loss, cntype.loss,
                              cntype.loh, cntype.gain),
                  cn.subtype = c(cntype.loh, cntype.gain, cntype.hetloss, cntype.homloss,
                                 cntype.loh, cntype.strongamp))
  n.1p <- score_lst(segs[seqnames(segs) == '1p'], cov)
  n.2p <- score_lst(segs[seqnames(segs) == '2p'], cov)
  n <- score_lst(segs, cov)
  expect_equal(c(n, n.1p, n.2p), c(7, 4, 3))
})

test_that("LST works - case 1", {
  cov <- GRanges(seqnames = factor(paste0(1:4, 'p'), levels = paste0(1:4, 'p')),
                 ranges = IRanges(start = rep(1, 4),
                                  end = rep(100*10^6, 4)))
  segs <- GRanges(seqnames = factor(rep('1p', 2), levels = paste0(1:4, 'p')),
                  ranges = IRanges(start = c(20,  30)*10^6 + 1,
                                   end =   c(22, 100)*10^6),
                  cn = c(3, 3),
                  cn.type = c(cntype.gain, cntype.gain),
                  cn.subtype = c(cntype.gain, cntype.gain))
  n <- score_lst(segs, cov)
  expect_equal(n, 1)
})

test_that("LST works - case 2", {
  cov <- GRanges(seqnames = factor(paste0(1:4, 'p'), levels = paste0(1:4, 'p')),
                 ranges = IRanges(start = rep(1, 4),
                                  end = rep(100*10^6, 4)))
  segs <- GRanges(seqnames = factor(rep('1p', 3), levels = paste0(1:4, 'p')),
                  ranges = IRanges(start = c(0,  20,  24)*10^6 + 1,
                                   end =   c(20, 22, 100)*10^6),
                  cn = c(1, 3, 3),
                  cn.type = c(cntype.loss, cntype.gain, cntype.gain),
                  cn.subtype = c(cntype.hetloss, cntype.gain, cntype.gain))
  n <- score_lst(segs, cov)
  expect_equal(n, 0)
})

test_that("LST works - case 3", {
  cov <- GRanges(seqnames = factor(paste0(1:4, 'p'), levels = paste0(1:4, 'p')),
                 ranges = IRanges(start = rep(1, 4),
                                  end = rep(100*10^6, 4)))
  segs <- GRanges(seqnames = factor(rep('1p', 2), levels = paste0(1:4, 'p')),
                  ranges = IRanges(start = c(20,  31)*10^6 + 1,
                                   end =   c(29, 100)*10^6),
                  cn = c(3, 7),
                  cn.type = c(cntype.gain, cntype.gain),
                  cn.subtype = c(cntype.gain, cntype.weakamp))
  n <- score_lst(segs, cov)
  expect_equal(n, 2)
})

test_that("LST works - case 4", {
  cov <- GRanges(seqnames = factor(paste0(1:4, 'p'), levels = paste0(1:4, 'p')),
                 ranges = IRanges(start = rep(1, 4),
                                  end = rep(100*10^6, 4)))
  segs <- GRanges(seqnames = factor(rep('1p', 2), levels = paste0(1:4, 'p')),
                  ranges = IRanges(start = c(43,  50)*10^6 + 1,
                                   end =   c(57, 100)*10^6),
                  cn = c(NA, 7),
                  cn.type = c(cntype.loh, cntype.gain),
                  cn.subtype = c(cntype.loh, cntype.weakamp))
  n <- score_lst(segs, cov)
  expect_equal(n, 0)
})

test_that("LST works - case 5", {
  cov <- GRanges(seqnames = factor(paste0(1:4, 'p'), levels = paste0(1:4, 'p')),
                 ranges = IRanges(start = rep(1, 4),
                                  end = rep(100*10^6, 4)))
  segs <- GRanges(seqnames = factor(rep('1p', 2), levels = paste0(1:4, 'p')),
                  ranges = IRanges(start = c(31, 61)*10^6 + 1,
                                   end =   c(80, 100)*10^6),
                  cn = c(3, 7),
                  cn.type = c(cntype.gain, cntype.gain),
                  cn.subtype = c(cntype.gain, cntype.weakamp))
  n <- score_lst(segs, cov)
  expect_equal(n, 3)
})

