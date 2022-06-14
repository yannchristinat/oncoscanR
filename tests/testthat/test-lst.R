test_that("LST works - simple cases", {
  cov <- GRanges(seqnames = factor(paste0(1:4, 'p'), levels = paste0(1:4, 'p')),
                 ranges = IRanges(start = rep(1, 4),
                                  end = rep(100*10^6, 4)))
  segs <- GRanges(seqnames = factor(c(rep('1p', 4), rep('2p', 2)), levels = paste0(1:4, 'p')),
                  ranges = IRanges(start = c(30, 50, 70,  90, 25,  50)*10^6 + 1,
                                   end =   c(50, 70, 90, 100, 75, 100)*10^6),
                  cn = c(NA, 3, 1, 0, NA, 14),
                  cn.type = c("LOH", "Gain", "Loss", "Loss",
                              "LOH", "Gain"))
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
                  cn.type = c("Gain", "Gain"))
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
                  cn.type = c("Loss", "Gain", "Gain"))
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
                  cn.type = c("Gain", "Gain"))
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
                  cn.type = c("LOH", "Gain"))
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
                  cn.type = c("Gain", "Gain"))
  n <- score_lst(segs, cov)
  expect_equal(n, 3)
})

test_that("LST works - real case", {
  oncoscan.cov <- oncoscanR::oncoscan_na33.cov[seqnames(oncoscanR::oncoscan_na33.cov) != '21p']

  chas.fn <- system.file("extdata", "LST_gene_list_full_location.txt",
                         package = "oncoscanR")
  segments <- load_chas(chas.fn, oncoscan.cov)

  segs.clean <- trim_to_coverage(segments, oncoscan.cov) %>%
    adjust_loh() %>%
    merge_segments() %>%
    prune_by_size()

  n <- score_lst(segs.clean, oncoscan.cov)
  expect_equal(n, 26) # Verified by hand in ChAS (on the 4 first chromosomes)
})

test_that("LST works - empty segments", {
  cov <- GRanges(seqnames = factor(paste0(1:4, 'p'), levels = paste0(1:4, 'p')),
                 ranges = IRanges(start = c(1,1,101,101),
                                  end = c(100,100,200,200)))
  segs <- GRanges(seqnames = factor(c(), levels = paste0(1:4, 'p')),
                  ranges = IRanges())
  n <- score_lst(segs, cov)
  expect_equal(n, 0)
})

