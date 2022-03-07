armlevel.test <- function(dat.armlevel, armlevel.fn){
  # Get ground truth
  dat.true <- read.csv(armlevel.fn, header = TRUE, row.names = 1, stringsAsFactors = FALSE)

  #Check AMP, GAIN, LOSS and LOH
  amp.test <- identical(sort(dat.armlevel[['AMP']]), sort(unlist(strsplit(dat.true['AMPL','Arms'], ','))))
  gain.test <- identical(sort(dat.armlevel[['GAIN']]), sort(unlist(strsplit(dat.true['GAIN','Arms'], ','))))
  loh.test <- identical(sort(dat.armlevel[['LOH']]), sort(unlist(strsplit(dat.true['LOH','Arms'], ','))))
  loss.test <- identical(sort(dat.armlevel[['LOSS']]), sort(unlist(strsplit(dat.true['LOSS','Arms'], ','))))
  return(c(amp.test, gain.test, loh.test, loss.test))
}
