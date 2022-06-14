test_that("Compare segments works - same segment, no cn info", {
  s1 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10))
  s2 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10))
  expect_true(same_segments(s1,s2))
})

test_that("Compare segments works - diff arm, no cn info", {
  s1 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10))
  s2 <- GRanges(seqnames = '1q', ranges = IRanges(start = 1, end = 10))
  expect_false(same_segments(s1,s2))
})

test_that("Compare segments works - diff end, no cn info", {
  s1 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10))
  s2 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 11))
  expect_false(same_segments(s1,s2))
})

test_that("Compare segments works - diff start, no cn info", {
  s1 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10))
  s2 <- GRanges(seqnames = '1p', ranges = IRanges(start = 2, end = 10))
  expect_false(same_segments(s1,s2))
})

test_that("Compare segments works - same segment, same cn", {
  s1 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10),
                cn = 5, cn.type = "Gain")
  s2 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10),
                cn = 5, cn.type = "Gain")
  expect_true(same_segments(s1,s2))
})

test_that("Compare segments works - same segment, diff cn", {
  s1 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10),
                cn = 5, cn.type = "Gain")
  s2 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10),
                cn = 6, cn.type = "Gain")
  expect_false(same_segments(s1,s2))
})

test_that("Compare segments works - same segment, NA cn", {
  s1 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10),
                cn = 5, cn.type = "Gain")
  s2 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10),
                cn = NA, cn.type = "Gain")
  expect_false(same_segments(s1,s2))
})

test_that("Compare segments works - same segment, diff cn.type", {
  s1 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10),
                cn = 5, cn.type = "Gain")
  s2 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10),
                cn = 5, cn.type = "Amplification")
  expect_false(same_segments(s1,s2))
})

test_that("Compare segments works - same segment, null cn", {
  s1 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10),
                cn = 5, cn.type = "Gain")
  s2 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10),
                cn.type = "Gain")
  expect_false(same_segments(s1,s2))
})

test_that("Compare segments works - same segment, null cn.type", {
  s1 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10),
                cn = 5, cn.type = "Gain")
  s2 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10),
                cn = 5)
  expect_false(same_segments(s1,s2))
})

test_that("Same set works - same segments/cn/cntype", {
  segs <- GRanges(seqnames = factor(c(rep('1p',3), rep('1q',3)),
                                    levels = c('1p', '1q', 'Xp')),
                  ranges = IRanges(start = c(1, 50, 130, 170, 250, 50),
                                   end = c(50, 130, 170, 250, 300, 250)),
                  cn = c(0,1,NA,3,5,10),
                  cn.type = c("Loss", "Loss", "LOH", "Gain",
                              "Gain", "Gain"))
  expect_true(same_segmentsets(segs, segs))
})

test_that("Same set works - different arm", {
  segs.a <- GRanges(seqnames = factor(c(rep('1p',3), rep('1q',3)),
                                      levels = c('1p', '1q', 'Xp')),
                    ranges = IRanges(start = c(1, 50, 130, 170, 250, 50),
                                     end = c(50, 130, 170, 250, 300, 250)),
                    cn = c(0,1,NA,3,5,10),
                    cn.type = c("Loss", "Loss", "LOH", "Gain",
                                "Gain", "Gain"))
  segs.b <- GRanges(seqnames = factor(c(rep('1q',3), rep('1q',3)),
                                      levels = c('1p', '1q', 'Xp')),
                    ranges = IRanges(start = c(1, 50, 130, 170, 250, 50),
                                     end = c(50, 130, 170, 250, 300, 250)),
                    cn = c(0,1,NA,3,5,10),
                    cn.type = c("Loss", "Loss", "LOH", "Gain",
                                "Gain", "Gain"))
  expect_false(same_segmentsets(segs.a, segs.b))
})

test_that("Same set works - different position", {
  segs.a <- GRanges(seqnames = factor(c(rep('1p',3), rep('1q',3)),
                                      levels = c('1p', '1q', 'Xp')),
                    ranges = IRanges(start = c(1, 50, 130, 170, 250, 50),
                                     end = c(50, 130, 170, 250, 300, 250)),
                    cn = c(0,1,NA,3,5,10),
                    cn.type = c("Loss", "Loss", "LOH", "Gain",
                                "Gain", "Gain"))
  segs.b <- GRanges(seqnames = factor(c(rep('1p',3), rep('1q',3)),
                                      levels = c('1p', '1q', 'Xp')),
                    ranges = IRanges(start = c(1, 50, 130, 170, 250, 50),
                                     end = c(51, 130, 170, 250, 300, 250)),
                    cn = c(0,1,NA,3,5,10),
                    cn.type = c("Loss", "Loss", "LOH", "Gain",
                                "Gain", "Gain"))
  expect_false(same_segmentsets(segs.a, segs.b))
})

test_that("Same set works - one segment missing", {
  segs.a <- GRanges(seqnames = factor(c(rep('1p',3), rep('1q',3)),
                                      levels = c('1p', '1q', 'Xp')),
                    ranges = IRanges(start = c(1, 50, 130, 170, 250, 50),
                                     end = c(50, 130, 170, 250, 300, 250)),
                    cn = c(0,1,NA,3,5,10),
                    cn.type = c("Loss", "Loss", "LOH", "Gain",
                                "Gain", "Gain"))
  segs.b <- GRanges(seqnames = factor(c(rep('1p',2), rep('1q',3)),
                                      levels = c('1p', '1q', 'Xp')),
                    ranges = IRanges(start = c(50, 130, 170, 250, 50),
                                     end = c(130, 170, 250, 300, 250)),
                    cn = c(1,NA,3,5,10),
                    cn.type = c("Loss", "LOH", "Gain",
                                "Gain", "Gain"))
  expect_false(same_segmentsets(segs.a, segs.b))
})

test_that("Same set works - different copy number", {
  segs.a <- GRanges(seqnames = factor(c(rep('1p',3), rep('1q',3)),
                                      levels = c('1p', '1q', 'Xp')),
                    ranges = IRanges(start = c(1, 50, 130, 170, 250, 50),
                                     end = c(50, 130, 170, 250, 300, 250)),
                    cn = c(0,1,NA,3,5,10),
                    cn.type = c("Loss", "Loss", "LOH", "Gain",
                                "Gain", "Gain"))
  segs.b <- GRanges(seqnames = factor(c(rep('1p',3), rep('1q',3)),
                                      levels = c('1p', '1q', 'Xp')),
                    ranges = IRanges(start = c(1, 50, 130, 170, 250, 50),
                                     end = c(50, 130, 170, 250, 300, 250)),
                    cn = c(NA,1,NA,3,5,10),
                    cn.type = c("LOH", "Loss", "LOH", "Gain",
                                "Gain", "Gain"))
  expect_false(same_segmentsets(segs.a, segs.b))
})

test_that("Same set works - different cntype", {
  segs.a <- GRanges(seqnames = factor(c(rep('1p',3), rep('1q',3)),
                                      levels = c('1p', '1q', 'Xp')),
                    ranges = IRanges(start = c(1, 50, 130, 170, 250, 50),
                                     end = c(50, 130, 170, 250, 300, 250)),
                    cn = c(0,1,NA,3,5,10),
                    cn.type = c("Loss", "Loss", "LOH", "Gain",
                                "Gain", "Gain"))
  segs.b <- GRanges(seqnames = factor(c(rep('1p',3), rep('1q',3)),
                                      levels = c('1p', '1q', 'Xp')),
                    ranges = IRanges(start = c(1, 50, 130, 170, 250, 50),
                                     end = c(50, 130, 170, 250, 300, 250)),
                    cn = c(0,1,NA,3,5,10),
                    cn.type = c("LOH", "Loss", "LOH", "Gain",
                                "Gain", "Gain"))
  expect_false(same_segmentsets(segs.a, segs.b))
})
