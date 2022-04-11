test_that("Trimming to coverage works", {
  cov <- GRanges(seqnames = factor(c('1p', '1q', 'Xp'), levels = c('1p', '1q', 'Xp')),
                 ranges = IRanges(start = c(100, 1000, 100),
                                  end = c(200, 2000, 200)))
  segs <- GRanges(seqnames = factor(rep('1p',6), levels = c('1p', '1q', 'Xp')),
                  ranges = IRanges(start = c(1, 50, 130, 170, 250, 50),
                                   end = c(50, 130, 170, 250, 300, 250)),
                  cn = c(0,1,NA,3,5,10),
                  cn.type = c(cntype.loss, cntype.loss, cntype.loh, cntype.gain,
                              cntype.gain, cntype.gain),
                  cn.subtype = c(cntype.homloss, cntype.hetloss, cntype.loh, cntype.gain,
                                 cntype.weakamp, cntype.strongamp))

  segs.clean <- trim_to_coverage(segs, cov)

  expected_segs <- GRanges(seqnames = factor(rep('1p',4), levels = c('1p', '1q', 'Xp')),
                           ranges = IRanges(start = c(100, 130, 170, 100),
                                            end = c(130, 170, 200, 200)),
                           cn = c(1,NA,3,10),
                           cn.type = c(cntype.loss, cntype.loh, cntype.gain, cntype.gain),
                           cn.subtype = c(cntype.hetloss, cntype.loh, cntype.gain, cntype.strongamp))

  expect_true(same_segmentsets(segs.clean, expected_segs))
})

test_that("Trimming with empty coverage works", {
  cov <- GRanges(seqnames = factor(c(), levels = c('1p', '1q', 'Xp')),
                 ranges = IRanges())
  segs <- GRanges(seqnames = factor(rep('1p',6), levels = c('1p', '1q', 'Xp')),
                  ranges = IRanges(start = c(1, 50, 130, 170, 250, 50),
                                   end = c(50, 130, 170, 250, 300, 250)),
                  cn = c(0,1,NA,3,5,10),
                  cn.type = c(cntype.loss, cntype.loss, cntype.loh, cntype.gain,
                              cntype.gain, cntype.gain),
                  cn.subtype = c(cntype.homloss, cntype.hetloss, cntype.loh, cntype.gain,
                                 cntype.weakamp, cntype.strongamp))

  segs.clean <- trim_to_coverage(segs, cov)
  expected_segs <- GRanges()

  expect_true(same_segmentsets(segs.clean, expected_segs))
})



test_that("Merging works with different cn", {
  segs <- GRanges(seqnames = factor(rep('1p',7), levels = c('1p', '1q', 'Xp')),
                  ranges = IRanges(start = c (1, 11, 29, 51, 69, 71, 94),
                                   end =   c(10, 20, 40, 60, 80, 90, 110)),
                  cn = c(0,0,0,0,5,NA,NA),
                  cn.type = c(cntype.loss, cntype.loss, cntype.loss, cntype.loss,
                              cntype.gain, cntype.loh, cntype.loh),
                  cn.subtype = c(cntype.homloss, cntype.homloss, cntype.homloss, cntype.homloss,
                                 cntype.weakamp, cntype.loh, cntype.loh))

  segs.clean <- merge_segments(segs, 10/1000)

  expected_segs <- GRanges(seqnames = factor(rep('1p',4), levels = c('1p', '1q', 'Xp')),
                           ranges = IRanges(start = c (1, 51, 69,  71),
                                            end =   c(40, 60, 80, 110)),
                           cn = c(0,0,5,NA),
                           cn.type = c(cntype.loss, cntype.loss, cntype.gain, cntype.loh),
                           cn.subtype = c(cntype.homloss, cntype.homloss, cntype.weakamp, cntype.loh))

  expect_true(same_segmentsets(segs.clean, expected_segs))
})

test_that("Merging works with overlaps", {
  segs <- GRanges(seqnames = factor(c(rep('1p',5), rep('2q',5)), levels = c('1p', '2q')),
                  ranges = IRanges(start = c( 1,  6, 51, 41, 71,  1,  6, 51, 41, 71),
                                   end =   c(10, 20, 60, 70, 80, 10, 20, 60, 70, 80)),
                  cn = c(0,0,0,0,0,3,NA,NA,1,1),
                  cn.type = c(rep(cntype.loss, 5),
                              cntype.gain, cntype.loh, cntype.loh, cntype.loss, cntype.loss),
                  cn.subtype = c(rep(cntype.homloss, 5),
                                 cntype.gain, cntype.loh, cntype.loh, cntype.hetloss, cntype.hetloss))

  segs.clean <- merge_segments(segs, 10/1000)

  expected_segs <- GRanges(seqnames = factor(c(rep('1p',2), rep('2q',4)), levels = c('1p', '2q')),
                           ranges = IRanges(start = c( 1, 41,  1,  6, 51, 41),
                                            end =   c(20, 80, 10, 20, 60, 80)),
                           cn = c(0,0,3,NA,NA,1),
                           cn.type = c(cntype.loss, cntype.loss, cntype.gain,
                                       cntype.loh, cntype.loh, cntype.loss),
                           cn.subtype = c(cntype.homloss, cntype.homloss, cntype.gain,
                                          cntype.loh, cntype.loh, cntype.hetloss))

  expect_true(same_segmentsets(segs.clean, expected_segs))
})

test_that("Merging works with large threshold", {
  segs <- GRanges(seqnames = factor(c(rep('1p',5), rep('2q',6)), levels = c('1p', '2q')),
                  ranges = IRanges(start = c( 1,  6, 51, 41, 71,  1,  6, 51, 41, 71, 91),
                                   end =   c(10, 20, 60, 70, 80, 10, 20, 60, 70, 80, 100)),
                  cn = c(0,0,0,0,0,3,NA,NA,1,1, 3),
                  cn.type = c(rep(cntype.loss, 5),
                              cntype.gain, cntype.loh, cntype.loh, cntype.loss, cntype.loss, cntype.gain),
                  cn.subtype = c(rep(cntype.homloss, 5),
                                 cntype.gain, cntype.loh, cntype.loh, cntype.hetloss, cntype.hetloss, cntype.gain))

  segs.clean <- merge_segments(segs, 1)

  expected_segs <- GRanges(seqnames = factor(c('1p', rep('2q',3)), levels = c('1p', '2q')),
                           ranges = IRanges(start = c( 1,   1,  6, 41),
                                            end =   c(80, 100, 60, 80)),
                           cn = c(0,3,NA,1),
                           cn.type = c(cntype.loss, cntype.gain, cntype.loh, cntype.loss),
                           cn.subtype = c(cntype.homloss, cntype.gain, cntype.loh, cntype.hetloss))

  expect_true(same_segmentsets(segs.clean, expected_segs))
})

test_that("Merging return an error with threshold=0", {
  segs <- GRanges(seqnames = factor(rep('1p',7), levels = c('1p', '1q', 'Xp')),
                  ranges = IRanges(start = c (1, 11, 29, 51, 69, 71, 94),
                                   end =   c(10, 20, 40, 60, 80, 90, 110)),
                  cn = c(0,0,0,0,5,NA,NA),
                  cn.type = c(cntype.loss, cntype.loss, cntype.loss, cntype.loss,
                              cntype.gain, cntype.loh, cntype.loh),
                  cn.subtype = c(cntype.homloss, cntype.homloss, cntype.homloss, cntype.homloss,
                                 cntype.weakamp, cntype.loh, cntype.loh))

  expect_error(merge_segments(segs, 0), 'greater than zero')
})

test_that("Merging works with min threshold (1)", {
  segs <- GRanges(seqnames = factor(c(rep('1p',5), rep('2q',5)), levels = c('1p', '2q')),
                  ranges = IRanges(start = c( 1,  6, 51, 41, 71,  1,  6, 22, 41, 71),
                                   end =   c(10, 20, 60, 70, 80, 10, 20, 60, 70, 80)),
                  cn = c(0,0,0,0,0,3,NA,NA,1,1),
                  cn.type = c(rep(cntype.loss, 5),
                              cntype.gain, cntype.loh, cntype.loh, cntype.loss, cntype.loss),
                  cn.subtype = c(rep(cntype.homloss, 5),
                                 cntype.gain, cntype.loh, cntype.loh, cntype.hetloss, cntype.hetloss))

  segs.clean <- merge_segments(segs, 1/1000)

  expected_segs <- GRanges(seqnames = factor(c(rep('1p',2), rep('2q',4)), levels = c('1p', '2q')),
                           ranges = IRanges(start = c( 1, 41,  1,  6, 22, 41),
                                            end =   c(20, 80, 10, 20, 60, 80)),
                           cn = c(0,0,3,NA,NA,1),
                           cn.type = c(cntype.loss, cntype.loss, cntype.gain,
                                       cntype.loh, cntype.loh, cntype.loss),
                           cn.subtype = c(cntype.homloss, cntype.homloss, cntype.gain,
                                          cntype.loh, cntype.loh, cntype.hetloss))

  expect_true(same_segmentsets(segs.clean, expected_segs))
})

#' Tests the following scenario (symbol + indicate that the LOH region should be trimmed)
#' LOH  ++++----
#' Loss ----
#' 
#' LOH  ----++++
#' Loss     ----
#' 
#' LOH  ----++++----
#' Loss     ----
#' 
#' LOH      ++++
#' Loss ------------
#' 
#' LOH      ----
#' Loss ----
#' 
#' LOH  ----
#' Loss     ----
#' 
#' LOH  ++++
#' Loss --------
#' 
#' LOH      ++++
#' Loss --------
#' 
#' LOH  ++++    ++++
#' Loss ------------
#' 
#' LOH  ++++----++++
#' Loss ----    ----
#' 
test_that("LOH are trimmed with respect to Loss", {
  segs_loh <- GRanges(seqnames = factor(rep('1p',11), levels = c('1p', '1q', 'Xp')),
                  ranges = IRanges(start = c( 1, 51,  81, 131, 171, 191, 221, 261, 281, 301, 321),
                                   end =   c(20, 70, 110, 140, 180, 200, 230, 270, 290, 310, 350)),
                  cn = rep(NA, 11),
                  cn.type = rep(cntype.loh, 11),
                  cn.subtype = rep(cntype.loh, 11))
  segs_loss <- GRanges(seqnames = factor(rep('1p',11), levels = c('1p', '1q', 'Xp')),
                      ranges = IRanges(start = c( 1, 61,  91, 121, 161, 201, 221, 251, 281, 321, 341),
                                       end =   c(10, 70, 100, 150, 170, 210, 240, 270, 310, 330, 350)),
                      cn = c(rep(0, 3), rep(1, 3), rep(1.5, 5)),
                      cn.type = rep(cntype.loss, 11),
                      cn.subtype = c(rep(cntype.homloss, 3), rep(cntype.hetloss, 8)))
  segs <- c(segs_loh, segs_loss)
  segs.clean <- adjust_loh(segs)

  expected_lohsegs <- GRanges(seqnames = factor(rep('1p',7), levels = c('1p', '1q', 'Xp')),
                           ranges = IRanges(start = c(11, 51, 81, 101, 171, 191, 331),
                                            end =   c(20, 60, 90, 110, 180, 200, 340)),
                           cn = rep(NA, 7),
                           cn.type = rep(cntype.loh, 7),
                           cn.subtype = rep(cntype.loh, 7))

  expect_true(same_segmentsets(segs.clean, c(expected_lohsegs, segs_loss)))
})


test_that("Pruning by size works", {
  segs <- GRanges(seqnames = factor(c('1p','2p','3p'), levels = c('1p','2p','3p')),
                  ranges = IRanges(start = c(1, 1, 1),
                                   end = c(49, 50, 51)),
                  cn = c(0,NA,5),
                  cn.type = c(cntype.loss, cntype.loh, cntype.gain),
                  cn.subtype = c(cntype.homloss, cntype.loh, cntype.weakamp))

  segs.clean <- prune_by_size(segs, 50/1000)

  expected_segs <- GRanges(seqnames = factor(c('2p','3p'), levels = c('1p','2p','3p')),
                  ranges = IRanges(start = c(1, 1),
                                   end = c(50, 51)),
                  cn = c(NA,5),
                  cn.type = c(cntype.loh, cntype.gain),
                  cn.subtype = c(cntype.loh, cntype.weakamp))

  expect_true(same_segmentsets(segs.clean, expected_segs))
})

test_that("Pruning by size works with zero threshold (no pruning)", {
  segs <- GRanges(seqnames = factor(c('1p','2p','3p'), levels = c('1p','2p','3p')),
                  ranges = IRanges(start = c(1, 1, 1),
                                   end = c(49, 50, 1)),
                  cn = c(0,NA,5),
                  cn.type = c(cntype.loss, cntype.loh, cntype.gain),
                  cn.subtype = c(cntype.homloss, cntype.loh, cntype.weakamp))

  segs.clean <- prune_by_size(segs, 0)

  expect_true(same_segmentsets(segs.clean, segs))
})

test_that("Cleaning with empty segments works", {
  cov <- GRanges(seqnames = factor(c('1p', '1q', 'Xp'), levels = c('1p', '1q', 'Xp')),
                 ranges = IRanges(start = c(100, 1000, 100),
                                  end = c(200, 2000, 200)))
  segs <- GRanges(seqnames = factor(c(), levels = c('1p', '1q', 'Xp')),
                  ranges = IRanges())
  segs.trimmed <- trim_to_coverage(segs, cov)
  segs.merged <- merge_segments(segs)
  segs.adj <- adjust_loh(segs)
  segs.pruned <- prune_by_size(segs)

  expect_true(same_segmentsets(segs, segs.trimmed))
  expect_true(same_segmentsets(segs, segs.merged))
  expect_true(same_segmentsets(segs, segs.adj))
  expect_true(same_segmentsets(segs, segs.pruned))
})
