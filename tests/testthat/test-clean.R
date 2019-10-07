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


test_that("LOH are trimmed with respect to Loss", {
  segs <- GRanges(seqnames = factor(rep('1p',12), levels = c('1p', '1q', 'Xp')),
                  ranges = IRanges(start = c (1,  1, 61, 51,  91,  81, 121, 131, 151, 161, 191, 181),
                                   end =   c(10, 20, 70, 70, 100, 110, 150, 140, 160, 170, 200, 190)),
                  cn = c(0,NA,1,NA,1,NA,1,NA,1,NA,1,NA),
                  cn.type = c(cntype.loss, cntype.loh,
                              cntype.loss, cntype.loh,
                              cntype.loss, cntype.loh,
                              cntype.loss, cntype.loh,
                              cntype.loss, cntype.loh,
                              cntype.loss, cntype.loh),
                  cn.subtype = c(cntype.homloss, cntype.loh,
                                 cntype.hetloss, cntype.loh,
                                 cntype.hetloss, cntype.loh,
                                 cntype.hetloss, cntype.loh,
                                 cntype.hetloss, cntype.loh,
                                 cntype.hetloss, cntype.loh))

  segs.clean <- adjust_loh(segs)

  expected_segs <- GRanges(seqnames = factor(rep('1p',11), levels = c('1p', '1q', 'Xp')),
                           ranges = IRanges(start = c (1,  11, 61, 51,  91,  81, 121, 151, 161, 191, 181),
                                            end =   c(10, 20, 70, 60, 100, 110, 150, 160, 170, 200, 190)),
                           cn = c(0,NA,1,NA,1,NA,1,1,NA,1,NA),
                           cn.type = c(cntype.loss, cntype.loh,
                                       cntype.loss, cntype.loh,
                                       cntype.loss, cntype.loh,
                                       cntype.loss,
                                       cntype.loss, cntype.loh,
                                       cntype.loss, cntype.loh),
                           cn.subtype = c(cntype.homloss, cntype.loh,
                                          cntype.hetloss, cntype.loh,
                                          cntype.hetloss, cntype.loh,
                                          cntype.hetloss,
                                          cntype.hetloss, cntype.loh,
                                          cntype.hetloss, cntype.loh))

  expect_true(same_segmentsets(segs.clean, expected_segs))
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
