# Get ChAS files to test
samples <- read.table(system.file("testdata", "validation_samples.txt", package = "oncoscanR"),
                      header=TRUE, row.names = 1)

armlevel.test <- function(chas.fn, gender, armlevel.fn){
  # Get ground truth
  dat.true <- read.csv(armlevel.fn, header = TRUE, row.names = 1, stringsAsFactors = FALSE)

  # Run standard workflow from package
  dat <- workflow_oncoscan.run(chas.fn, gender)
  dat.armlevel <- dat[['armlevel']]

  #Check AMP, GAIN, LOSS and LOH
  amp.test <- identical(sort(dat.armlevel[['AMP']]), sort(unlist(strsplit(dat.true['AMPL','Arms'], ','))))
  gain.test <- identical(sort(dat.armlevel[['GAIN']]), sort(unlist(strsplit(dat.true['GAIN','Arms'], ','))))
  loh.test <- identical(sort(dat.armlevel[['LOH']]), sort(unlist(strsplit(dat.true['LOH','Arms'], ','))))
  loss.test <- identical(sort(dat.armlevel[['LOSS']]), sort(unlist(strsplit(dat.true['LOSS','Arms'], ','))))
  return(c(amp.test, gain.test, loh.test, loss.test))
}
