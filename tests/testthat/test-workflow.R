test_that("Oncoscan workflow works - H18015654", {
  skip_on_cran()

  sample <- 'H18015654'
  tests <- armlevel.test(
    system.file("testdata", samples[sample, 'chas.fn'], package = "oncoscanR"),
    samples[sample, 'gender'],
    system.file("testdata", samples[sample, 'armlevel.fn'], package = "oncoscanR"))
  expect_true(sum(tests)==4)
})

test_that("Oncoscan workflow works - H19005015", {
  skip_on_cran()

  sample <- 'H19005015'
  tests <- armlevel.test(
    system.file("testdata", samples[sample, 'chas.fn'], package = "oncoscanR"),
    samples[sample, 'gender'],
    system.file("testdata", samples[sample, 'armlevel.fn'], package = "oncoscanR"))
  expect_true(sum(tests)==4)
})

test_that("Oncoscan workflow works - H19005036", {
  skip_on_cran()

  sample <- 'H19005036'
  tests <- armlevel.test(
    system.file("testdata", samples[sample, 'chas.fn'], package = "oncoscanR"),
    samples[sample, 'gender'],
    system.file("testdata", samples[sample, 'armlevel.fn'], package = "oncoscanR"))
  expect_true(sum(tests)==4)
})

test_that("Oncoscan workflow works - H19005761", {
  skip_on_cran()

  sample <- 'H19005761'
  tests <- armlevel.test(
    system.file("testdata", samples[sample, 'chas.fn'], package = "oncoscanR"),
    samples[sample, 'gender'],
    system.file("testdata", samples[sample, 'armlevel.fn'], package = "oncoscanR"))
  expect_true(sum(tests)==4)
})

test_that("Oncoscan workflow works - H19005930", {
  skip_on_cran()

  sample <- 'H19005930'
  tests <- armlevel.test(
    system.file("testdata", samples[sample, 'chas.fn'], package = "oncoscanR"),
    samples[sample, 'gender'],
    system.file("testdata", samples[sample, 'armlevel.fn'], package = "oncoscanR"))
  expect_true(sum(tests)==4)
})

test_that("Oncoscan workflow works - H19005931", {
  skip_on_cran()

  sample <- 'H19005931'
  tests <- armlevel.test(
    system.file("testdata", samples[sample, 'chas.fn'], package = "oncoscanR"),
    samples[sample, 'gender'],
    system.file("testdata", samples[sample, 'armlevel.fn'], package = "oncoscanR"))
  expect_true(sum(tests)==4)
})

test_that("Oncoscan workflow works - H19006144", {
  skip_on_cran()

  sample <- 'H19006144'
  tests <- armlevel.test(
    system.file("testdata", samples[sample, 'chas.fn'], package = "oncoscanR"),
    samples[sample, 'gender'],
    system.file("testdata", samples[sample, 'armlevel.fn'], package = "oncoscanR"))
  expect_true(sum(tests)==4)
})

test_that("Oncoscan workflow works - H19006249", {
  skip_on_cran()

  sample <- 'H19006249'
  tests <- armlevel.test(
    system.file("testdata", samples[sample, 'chas.fn'], package = "oncoscanR"),
    samples[sample, 'gender'],
    system.file("testdata", samples[sample, 'armlevel.fn'], package = "oncoscanR"))
  expect_true(sum(tests)==4)
})

test_that("Oncoscan workflow works - H19006250", {
  skip_on_cran()

  sample <- 'H19006250'
  tests <- armlevel.test(
    system.file("testdata", samples[sample, 'chas.fn'], package = "oncoscanR"),
    samples[sample, 'gender'],
    system.file("testdata", samples[sample, 'armlevel.fn'], package = "oncoscanR"))
  expect_true(sum(tests)==4)
})

test_that("Oncoscan workflow works - H19006593", {
  skip_on_cran()

  sample <- 'H19006593'
  tests <- armlevel.test(
    system.file("testdata", samples[sample, 'chas.fn'], package = "oncoscanR"),
    samples[sample, 'gender'],
    system.file("testdata", samples[sample, 'armlevel.fn'], package = "oncoscanR"))
  expect_true(sum(tests)==4)
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
