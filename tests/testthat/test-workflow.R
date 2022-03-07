test_that("Oncoscan workflow works - LST", {
  skip_on_cran()

  chas.fn <- system.file("testdata", "LST_gene_list_full_location.txt", package = "oncoscanR")
  armlevel.fn <- system.file("testdata", "LST_gene_list_full_location.armlevel_scna.csv", package = "oncoscanR")

  # Run standard workflow from package
  dat <- workflow_oncoscan.run(chas.fn, 'F')

  tests <- armlevel.test(dat[['armlevel']], armlevel.fn)
  expect_true(sum(tests)==4)

  scores <- dat[['scores']]
  expect_true(scores['HRD'] == "Positive, nLST=22.5")
  expect_true(scores['avgCN'] == '2.64')
  expect_true(scores['TDplus'] == '131')
})

test_that("Oncoscan workflow works - triploide", {
  skip_on_cran()

  chas.fn <- system.file("testdata", "triploide_gene_list_full_location.txt", package = "oncoscanR")
  armlevel.fn <- system.file("testdata", "triploide_gene_list_full_location.armlevel_scna.csv", package = "oncoscanR")

  # Run standard workflow from package
  dat <- workflow_oncoscan.run(chas.fn, 'F')

  tests <- armlevel.test(dat[['armlevel']], armlevel.fn)
  expect_true(sum(tests)==4)

  scores <- dat[['scores']]
  expect_true(scores['HRD'] == "Negative, nLST=0")
  expect_true(scores['avgCN'] == '3.24')
  expect_true(scores['TDplus'] == '18')
})

test_that("Oncoscan workflow works - TDplus", {
  skip_on_cran()

  chas.fn <- system.file("testdata", "TDplus_gene_list_full_location.txt", package = "oncoscanR")
  armlevel.fn <- system.file("testdata", "TDplus_gene_list_full_location.armlevel_scna.csv", package = "oncoscanR")

  # Run standard workflow from package
  dat <- workflow_oncoscan.run(chas.fn, 'F')

  tests <- armlevel.test(dat[['armlevel']], armlevel.fn)
  expect_true(sum(tests)==4)

  scores <- dat[['scores']]
  expect_true(scores['HRD'] == "Negative, nLST=10.5")
  expect_true(scores['avgCN'] == '3.21')
  expect_true(scores['TDplus'] == '94')
})

test_that("Oncoscan workflow returns the correct structure", {
  segs.filename <- system.file("extdata", "chas_example.txt", package = "oncoscanR")
  dat <- workflow_oncoscan.run(segs.filename, 'M')

  armlevel.nametest <- identical(sort(names(dat[['armlevel']])), sort(c('AMP', 'GAIN', 'LOH', 'LOSS')))
  scores.nametest <- identical(sort(names(dat[['scores']])), sort(c("avgCN", "HRD", "TDplus")))
  gender.test <- dat[['gender']] == 'M'
  file.test <- dat[['file']] == basename(segs.filename)
  #print(c(armlevel.nametest, scores.nametest, gender.test, file.test))

  expect_true(sum(c(armlevel.nametest, scores.nametest, gender.test, file.test))==4)
})
