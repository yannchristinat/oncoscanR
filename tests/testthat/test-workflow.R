test_that("Oncoscan workflow works - LST", {
  #skip_on_cran()

  chas.fn <- "../testdata/LST_gene_list_full_location.txt"
  armlevel.fn <- "../testdata/LST_gene_list_full_location.armlevel_scna.csv"

  # Run standard workflow from package
  dat <- workflow_oncoscan.run(chas.fn)

  tests <- armlevel.test(dat[['armlevel']], armlevel.fn)
  expect_true(sum(tests)==4)

  scores <- dat[['scores']]
  expect_equal(as.character(scores['HRD']), "Positive, nLST=22.5")
  expect_equal(as.character(scores['avgCN']), '2.64')
  expect_equal(as.character(scores['TDplus']), '139')
})

test_that("Oncoscan workflow works - triploide", {
  #skip_on_cran()

  chas.fn <- "../testdata/triploide_gene_list_full_location.txt"
  armlevel.fn <- "../testdata/triploide_gene_list_full_location.armlevel_scna.csv"

  # Run standard workflow from package
  dat <- workflow_oncoscan.run(chas.fn)

  tests <- armlevel.test(dat[['armlevel']], armlevel.fn)
  expect_true(sum(tests)==4)

  scores <- dat[['scores']]
  expect_equal(as.character(scores['HRD']), "Negative, nLST=0")
  expect_equal(as.character(scores['avgCN']), '3.24')
  expect_equal(as.character(scores['TDplus']), '18')
})

test_that("Oncoscan workflow works - TDplus", {
  #skip_on_cran()

  chas.fn <- "../testdata/TDplus_gene_list_full_location.txt"
  armlevel.fn <- "../testdata/TDplus_gene_list_full_location.armlevel_scna.csv"

  # Run standard workflow from package
  dat <- workflow_oncoscan.run(chas.fn)

  tests <- armlevel.test(dat[['armlevel']], armlevel.fn)
  expect_true(sum(tests)==4)

  scores <- dat[['scores']]
  expect_equal(as.character(scores['HRD']), "Negative, nLST=10.5")
  expect_equal(as.character(scores['avgCN']), '3.21')
  expect_equal(as.character(scores['TDplus']), '96')
})

test_that("Oncoscan workflow returns the correct structure", {
  segs.filename <- system.file("extdata", "chas_example.txt", package = "oncoscanR")
  dat <- workflow_oncoscan.run(segs.filename)

  armlevel.nametest <- identical(sort(names(dat[['armlevel']])), sort(c('AMP', 'GAIN', 'LOH', 'LOSS')))
  scores.nametest <- identical(sort(names(dat[['scores']])), sort(c("avgCN", "HRD", "TDplus")))
  file.test <- dat[['file']] == basename(segs.filename)
  #print(c(armlevel.nametest, scores.nametest, file.test))

  expect_true(sum(c(armlevel.nametest, scores.nametest, file.test))==3)
})

test_that("Oncoscan workflow works with an empty file", {
  segs.filename <- "../testdata/chas_example-empty.txt"
  armlevel.fn <- "../testdata/chas_example-empty.armlevel_scna.csv"

  expect_warning(dat <- workflow_oncoscan.run(segs.filename))

  # Test arm-level detection
  tests <- armlevel.test(dat[['armlevel']], armlevel.fn)
  expect_true(sum(tests)==4)

  # Test scores
  scores <- unlist(dat[['scores']])
  expect_equal(as.character(scores['HRD']), "Negative (no tumor?), nLST=0")
  expect_equal(as.character(scores['avgCN']), '2')
  expect_equal(as.character(scores['TDplus']), '0')
})
