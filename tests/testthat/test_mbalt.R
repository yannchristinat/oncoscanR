test_that("Mbp alt works", {
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

    mbalt.wLOH <- score_mbalt(segs, cov, loh.rm = FALSE)
    expect_equal(mbalt.wLOH, c(sample=145, kit=400))

    mbalt.noLOH <- score_mbalt(segs, cov, loh.rm = TRUE)
    expect_equal(mbalt.noLOH, c(sample=100, kit=400))

})
