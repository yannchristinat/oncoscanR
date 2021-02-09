test_that("CN subtype 'LOH' works for male", {
  segs.loh <- GRanges(seqnames = c('1p','Xq','chrYq'),
                      ranges = IRanges(start = c(1,1,1), end = c(10, 10, 10)))
  segs.loh$cn <- NA
  segs.loh$cn.type <- cntype.loh
  subtypes <- get_cn_subtype(segs.loh, 'M')
  expect_equal(sum(subtypes == cntype.loh), 3)
})

test_that("CN subtype 'Het loss' works for male", {
  segs.hetloss <- GRanges(seqnames = c('1p','3q','chr22q'),
                          ranges = IRanges(start = c(1,1,1), end = c(10, 10, 10)))
  segs.hetloss$cn <- 1
  segs.hetloss$cn.type <- cntype.loss
  subtypes <- get_cn_subtype(segs.hetloss, 'M')
  expect_equal(sum(subtypes == cntype.hetloss), 3)
})

test_that("CN subtype 'Hom loss' works for male", {
  segs.homloss <- GRanges(seqnames = c('1p','Xq','chrYq'),
                          ranges = IRanges(start = c(1,1,1), end = c(10, 10, 10)))
  segs.homloss$cn <- 0
  segs.homloss$cn.type <- cntype.loss
  subtypes <- get_cn_subtype(segs.homloss, 'M')
  expect_equal(sum(subtypes == cntype.homloss), 3)
})

test_that("CN subtype 'Gain' works for male", {
  segs.gain <- GRanges(seqnames = c('1p','Xq','chrYq'),
                       ranges = IRanges(start = c(1,1,1), end = c(10, 10, 10)))
  segs.gain$cn <- c(3,3,2)
  segs.gain$cn.type <- cntype.gain
  subtypes <- get_cn_subtype(segs.gain, 'M')
  expect_equal(sum(subtypes == cntype.gain), 3)
})

test_that("CN subtype 'Weak amp' works for male", {
  segs.weakamp <- GRanges(seqnames = c('1p','Xq','chrYq'),
                          ranges = IRanges(start = c(1,1,1), end = c(10, 10, 10)))
  segs.weakamp$cn <- c(5,4,8)
  segs.weakamp$cn.type <- cntype.gain
  subtypes <- get_cn_subtype(segs.weakamp, 'M')
  expect_equal(sum(subtypes == cntype.weakamp), 3)
})

test_that("CN subtype 'Strong amp' works for male", {
  segs.strongamp <- GRanges(seqnames = c('1p','Xq','chrYq'),
                            ranges = IRanges(start = c(1,1,1), end = c(10, 10, 10)))
  segs.strongamp$cn <- c(10,15,9)
  segs.strongamp$cn.type <- cntype.gain
  subtypes <- get_cn_subtype(segs.strongamp, 'M')
  expect_equal(sum(subtypes == cntype.strongamp), 3)
})


#
# Same as above but for female
#
test_that("CN subtype 'LOH' works for female", {
  segs.loh <- GRanges(seqnames = c('1p','Xq','chrXp'),
                      ranges = IRanges(start = c(1,1,1), end = c(10, 10, 10)))
  segs.loh$cn <- NA
  segs.loh$cn.type <- cntype.loh
  subtypes <- get_cn_subtype(segs.loh, 'F')
  expect_equal(sum(subtypes == cntype.loh), 3)
})

test_that("CN subtype 'Het loss' works for female", {
  segs.hetloss <- GRanges(seqnames = c('1p','Xq','chrXp'),
                          ranges = IRanges(start = c(1,1,1), end = c(10, 10, 10)))
  segs.hetloss$cn <- 1
  segs.hetloss$cn.type <- cntype.loss
  subtypes <- get_cn_subtype(segs.hetloss, 'F')
  expect_equal(sum(subtypes == cntype.hetloss), 3)
})

test_that("CN subtype 'Hom loss' works for female", {
  segs.homloss <- GRanges(seqnames = c('1p','Xq','chrXp'),
                          ranges = IRanges(start = c(1,1,1), end = c(10, 10, 10)))
  segs.homloss$cn <- 0
  segs.homloss$cn.type <- cntype.loss
  subtypes <- get_cn_subtype(segs.homloss, 'F')
  expect_equal(sum(subtypes == cntype.homloss), 3)
})

test_that("CN subtype 'Gain' works for female", {
  segs.gain <- GRanges(seqnames = c('1p','Xq','chrXp'),
                       ranges = IRanges(start = c(1,1,1), end = c(10, 10, 10)))
  segs.gain$cn <- c(3,3,4)
  segs.gain$cn.type <- cntype.gain
  subtypes <- get_cn_subtype(segs.gain, 'F')
  expect_equal(sum(subtypes == cntype.gain), 3)
})

test_that("CN subtype 'Weak amp' works for female", {
  segs.weakamp <- GRanges(seqnames = c('1p','Xq','chrXp'),
                          ranges = IRanges(start = c(1,1,1), end = c(10, 10, 10)))
  segs.weakamp$cn <- c(5,5,9)
  segs.weakamp$cn.type <- cntype.gain
  subtypes <- get_cn_subtype(segs.weakamp, 'F')
  expect_equal(sum(subtypes == cntype.weakamp), 3)
})

test_that("CN subtype 'Strong amp' works for female", {
  segs.strongamp <- GRanges(seqnames = c('1p','Xq','chrXp'),
                            ranges = IRanges(start = c(1,1,1), end = c(10, 10, 10)))
  segs.strongamp$cn <- c(10,15,10)
  segs.strongamp$cn.type <- cntype.gain
  subtypes <- get_cn_subtype(segs.strongamp, 'F')
  expect_equal(sum(subtypes == cntype.strongamp), 3)
})


#
# Test if no gender
#
test_that("CN subtype 'LOH' works for no gender", {
  segs.loh <- GRanges(seqnames = c('1p','Xq','chrYq'),
                      ranges = IRanges(start = c(1,1,1), end = c(10, 10, 10)))
  segs.loh$cn <- NA
  segs.loh$cn.type <- cntype.loh
  subtypes <- get_cn_subtype(segs.loh, NA)
  expect_equal(sum(subtypes == cntype.loh), 1)
})

test_that("CN subtype 'Het loss' works for no gender", {
  segs.hetloss <- GRanges(seqnames = c('1p','Xq','chrYq'),
                          ranges = IRanges(start = c(1,1,1), end = c(10, 10, 10)))
  segs.hetloss$cn <- 1
  segs.hetloss$cn.type <- cntype.loss
  subtypes <- get_cn_subtype(segs.hetloss, NA)
  expect_equal(sum(subtypes == cntype.hetloss), 1)
})

test_that("CN subtype 'Hom loss' works for no gender", {
  segs.homloss <- GRanges(seqnames = c('1p','Xq','chrYq'),
                          ranges = IRanges(start = c(1,1,1), end = c(10, 10, 10)))
  segs.homloss$cn <- 0
  segs.homloss$cn.type <- cntype.loss
  subtypes <- get_cn_subtype(segs.homloss, NA)
  expect_equal(sum(subtypes == cntype.homloss), 1)
})

test_that("CN subtype 'Gain' works for no gender", {
  segs.gain <- GRanges(seqnames = c('1p','Xq','chrYq'),
                       ranges = IRanges(start = c(1,1,1), end = c(10, 10, 10)))
  segs.gain$cn <- c(3,3,2)
  segs.gain$cn.type <- cntype.gain
  subtypes <- get_cn_subtype(segs.gain, NA)
  expect_equal(sum(subtypes == cntype.gain), 1)
})

test_that("CN subtype 'Weak amp' works for no gender", {
  segs.weakamp <- GRanges(seqnames = c('1p','Xq','chrYq'),
                          ranges = IRanges(start = c(1,1,1), end = c(10, 10, 10)))
  segs.weakamp$cn <- c(5,4,8)
  segs.weakamp$cn.type <- cntype.gain
  subtypes <- get_cn_subtype(segs.weakamp, NA)
  expect_equal(sum(subtypes == cntype.weakamp), 1)
})

test_that("CN subtype 'Strong amp' works for no gender", {
  segs.strongamp <- GRanges(seqnames = c('1p','Xq','chrYq'),
                            ranges = IRanges(start = c(1,1,1), end = c(10, 10, 10)))
  segs.strongamp$cn <- c(10,15,9)
  segs.strongamp$cn.type <- cntype.gain
  subtypes <- get_cn_subtype(segs.strongamp, NA)
  expect_equal(sum(subtypes == cntype.strongamp), 1)
})
