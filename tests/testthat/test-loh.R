test_that("HR-LOH score works - simple cases", {
  cov <- GRanges(seqnames = factor(c(paste0(1:4, 'p'), paste0(1:4, 'q')),
                                   levels = c(paste0(1:4, 'p'), paste0(1:4, 'q'))),
                 ranges = IRanges(start = rep(1, 8),
                                  end = rep(200*10^6, 8)))
  segs <- GRanges(seqnames = factor(c(rep('1p', 6), rep('2p', 2)), levels = paste0(1:4, 'p')),
                  ranges = IRanges(start = c( 1, 10, 40, 60-0.000001, 80, 100,  1, 50)*10^6 + 1,
                                   end =   c(10, 30, 55, 75,        95, 120, 20, 80)*10^6),
                  cn = c(NA, NA, NA, NA, 3, 1, NA, 15),
                  cn.type = c("LOH", "LOH", "LOH", "LOH",
                              "Gain", "Loss", "LOH", "Gain"))
  n <- score_loh(segs, c('2p', '2q'), c(), cov)
  expect_equal(n, 3)
})

test_that("HR-LOH score works - chrom with LOH and hetloss", {
  cov <- GRanges(seqnames = factor(c(paste0(1:4, 'p'), paste0(1:4, 'q')),
                                   levels = c(paste0(1:4, 'p'), paste0(1:4, 'q'))),
                 ranges = IRanges(start = rep(1, 8),
                                  end = rep(200*10^6, 8)))
  segs <- GRanges(seqnames = factor(c(rep('1p', 6), rep('2p', 2)), levels = paste0(1:4, 'p')),
                  ranges = IRanges(start = c( 1, 10, 40, 60-0.000001, 80, 100,  1, 50)*10^6 + 1,
                                   end =   c(10, 30, 55, 75,        95, 120, 20, 80)*10^6),
                  cn = c(NA, NA, NA, NA, 3, 1, NA, 15),
                  cn.type = c("LOH", "LOH", "LOH", "LOH",
                              "Gain", "Loss", "LOH", "Gain"))
  n <- score_loh(segs, c('2p'), c('2q'), cov)
  expect_equal(n, 4)
})


test_that("HR-LOH score works - overlapping loh", {
  cov <- GRanges(seqnames = factor(c(paste0(1:4, 'p'), paste0(1:4, 'q')),
                                   levels = c(paste0(1:4, 'p'), paste0(1:4, 'q'))),
                 ranges = IRanges(start = rep(1, 8),
                                  end = rep(200*10^6, 8)))
  segs <- GRanges(seqnames = factor(c(rep('1p', 6), rep('2p', 2)), levels = paste0(1:4, 'p')),
                  ranges = IRanges(start = c( 1, 10, 40, 60-0.000001, 80, 100,  1, 50)*10^6 + 1,
                                   end =   c(10, 20, 55, 75,        95, 120, 20, 80)*10^6),
                  cn = c(NA, NA, NA, NA, 3, 1, NA, 15),
                  cn.type = c("LOH", "LOH", "LOH", "LOH",
                              "Gain", "Loss", "LOH", "Gain"))
  n <- score_loh(segs, c(), c(), cov)
  expect_equal(n, 4)
})

test_that("LOH works - real case", {
  oncoscan.cov <- oncoscanR::oncoscan_na33.cov[seqnames(oncoscanR::oncoscan_na33.cov) != '21p']

  chas.fn <- "../testdata/LST_gene_list_full_location.txt"
  segments <- load_chas(chas.fn, oncoscan.cov)

  segs.clean <- trim_to_coverage(segments, oncoscan.cov) %>%
    adjust_loh() %>%
    prune_by_size()

  armlevel.loh <- get_loh_segments(segs.clean) %>%
    armlevel_alt(kit.coverage = oncoscan.cov)
  armlevel.hetloss <- get_hetloss_segments(segs.clean) %>%
    armlevel_alt(kit.coverage = oncoscan.cov)

  n <- score_loh(segs.clean, names(armlevel.loh), names(armlevel.hetloss), oncoscan.cov)
  expect_equal(n, 25) #Verified by hand in ChAS
})

test_that("LOH works - empty segments", {
  cov <- GRanges(seqnames = factor(paste0(1:4, 'p'), levels = paste0(1:4, 'p')),
                 ranges = IRanges(start = c(1,1,101,101),
                                  end = c(100,100,200,200)))
  segs <- GRanges(seqnames = factor(c(), levels = paste0(1:4, 'p')),
                  ranges = IRanges())
  n <- score_loh(segs, c(), c(), cov)
  expect_equal(n, 0)
})

test_that("gLOH works - real case", {
  oncoscan.cov <- oncoscanR::oncoscan_na33.cov[seqnames(oncoscanR::oncoscan_na33.cov) != '21p']

  chas.fn <- "../testdata/LST_gene_list_full_location.txt"
  segments <- load_chas(chas.fn, oncoscan.cov)
  segs.clean <- trim_to_coverage(segments, oncoscan.cov) %>%
    adjust_loh() %>%
    prune_by_size()

  armlevel.loh <- get_loh_segments(segs.clean) %>%
    armlevel_alt(kit.coverage = oncoscan.cov)
  armlevel.hetloss <- get_hetloss_segments(segs.clean) %>%
    armlevel_alt(kit.coverage = oncoscan.cov)

  sel.segs <- segs.clean$cn.type=="LOH" | (segs.clean$cn.type=="Loss" & segs.clean$cn==1)
  sel.arms <- !as.vector(seqnames(segs.clean)) %in% c(names(armlevel.loh), names(armlevel.hetloss))
  expected_p <- sum( width(segs.clean[sel.segs & sel.arms]))/sum(width(oncoscan.cov))

  p <- score_gloh(segs.clean, names(armlevel.loh), names(armlevel.hetloss), oncoscan.cov)
  expect_equal(p, expected_p)
})

