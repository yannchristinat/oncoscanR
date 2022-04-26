#' On a simple case with 4 chromosomes and 6 alterations (including an overlapping LOH segment)
#' TOT bp = 4*100 = 400
#' ALT bp weigthed by CN = 3*20 + 1*20 + 0*10 + 14*50 = 780
#' NON-ALT bp weigthed by CN = 2*(TOT bp - (20+20+10+50)) = 2*300 = 600
#' avg CN = (780+600)/400 = 3.45
test_that("Avg CN computes correctly", {
    cov <- GRanges(seqnames = factor(paste0(1:4, 'p'), levels = paste0(1:4, 'p')),
                   ranges = IRanges(start = rep(1, 4),
                                    end = rep(100, 4)))
    segs <- GRanges(seqnames = factor(c(rep('1p', 4), rep('2p', 2)), levels = paste0(1:4, 'p')),
                    ranges = IRanges(start = c(30, 50, 70,  90, 25,  50) + 1,
                                     end =   c(50, 70, 90, 100, 75, 100)),
                    cn = c(NA, 3, 1, 0, NA, 14),
                    cn.type = c(cntype.loh, cntype.gain, cntype.loss, cntype.loss,
                                cntype.loh, cntype.gain),
                    cn.subtype = c(cntype.loh, cntype.gain, cntype.hetloss, cntype.homloss,
                                   cntype.loh, cntype.strongamp))
    avgcn <- score_avgcn(segs, cov)
    expect_equal(avgcn, 3.45)
})

#' Test the limits at 0, 2, 2.19, 2.2, 2.21, 3.38, 3.4, 3.42 and 10
test_that("Est WGD works on limit cases", {
    cov <- GRanges(seqnames = factor('1p', levels = c('1p')),
                   ranges = IRanges(start = 1, end = 100))
    # Test at avgCN=0
    seg <- GRanges(seqnames = factor('1p', levels = c('1p')),
                ranges = IRanges(start = 1, end = 100), cn=0, cn.type = cntype.loss)
    est <- score_estwgd(seg, cov)
    expect_equal(est, c(WGD=0,avgCN=0))

    # Test at avgCN=2 (no segments)
    est <- score_estwgd(GRanges(), cov)
    expect_equal(est, c(WGD=0,avgCN=2))

    # Test at avgCN=2.19
    seg <- GRanges(seqnames = factor('1p', levels = c('1p')),
                   ranges = IRanges(start = 1, end = 19), cn=3, cn.type = cntype.gain)
    est <- score_estwgd(seg, cov)
    expect_equal(est, c(WGD=0,avgCN=2.19))

    # Test at avgCN=2.2
    seg <- GRanges(seqnames = factor('1p', levels = c('1p')),
                   ranges = IRanges(start = 1, end = 20), cn=3, cn.type = cntype.gain)
    est <- score_estwgd(seg, cov)
    expect_equal(est, c(WGD=0,avgCN=2.2))

    # Test at avgCN=2.21
    seg <- GRanges(seqnames = factor('1p', levels = c('1p')),
                   ranges = IRanges(start = 1, end = 21), cn=3, cn.type = cntype.gain)
    est <- score_estwgd(seg, cov)
    expect_equal(est, c(WGD=1,avgCN=2.21))

    # Test at avgCN=3.38
    seg <- GRanges(seqnames = factor('1p', levels = c('1p')),
                   ranges = IRanges(start = 1, end = 69), cn=4, cn.type = cntype.gain)
    est <- score_estwgd(seg, cov)
    expect_equal(est, c(WGD=1,avgCN=3.38))

    # Test at avgCN=3.4
    seg <- GRanges(seqnames = factor('1p', levels = c('1p')),
                   ranges = IRanges(start = 1, end = 70), cn=4, cn.type = cntype.gain)
    est <- score_estwgd(seg, cov)
    expect_equal(est, c(WGD=1,avgCN=3.4))

    # Test at avgCN=3.42
    seg <- GRanges(seqnames = factor('1p', levels = c('1p')),
                   ranges = IRanges(start = 1, end = 71), cn=4, cn.type = cntype.gain)
    est <- score_estwgd(seg, cov)
    expect_equal(est, c(WGD=2,avgCN=3.42))

    # Test at avgCN=10
    seg <- GRanges(seqnames = factor('1p', levels = c('1p')),
                   ranges = IRanges(start = 1, end = 100), cn=10, cn.type = cntype.gain)
    est <- score_estwgd(seg, cov)
    expect_equal(est, c(WGD=2,avgCN=10))
})


test_that("nLST works - simple case and limit cases", {
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

    # Test with WGD=0, threshold=15 (default)
    n <- score_nlst(segs, 0, cov)
    expect_equal(n, c(nLST=4, HRD='Negative'))

    # Test with WGD=0, threshold=4
    n <- score_nlst(segs, 0, cov, 4)
    expect_equal(n, c(nLST=4, HRD='Positive'))

    # Test with WGD=0, threshold=3.5
    n <- score_nlst(segs, 0, cov, 3.5)
    expect_equal(n, c(nLST=4, HRD='Positive'))

    # Test with WGD=1, threshold=15 (default)
    n <- score_nlst(segs, 1, cov)
    expect_equal(n, c(nLST=0.5, HRD='Negative'))

    # Test with WGD=2, threshold=15 (default)
    n <- score_nlst(segs, 2, cov)
    expect_equal(n, c(nLST=0, HRD='Negative'))

    # Test with WGD=3, threshold=15 (default)
    n <- score_nlst(segs, 3, cov)
    expect_equal(n, c(nLST=0, HRD='Negative'))

})

