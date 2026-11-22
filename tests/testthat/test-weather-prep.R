test_that("prepare_weather adds required derived columns", {
  x <- data.frame(date = c("2000-01-01", "2000-01-02"), tmin = c(1, 2), tmax = c(5, 6))
  y <- prepare_weather(x)
  expect_true(all(c("year", "doy", "tmean") %in% names(y)))
  expect_equal(y$tmean, c(3, 4))
})

test_that("prepare_weather rejects missing columns", {
  expect_error(prepare_weather(data.frame(date = "2000-01-01", tmin = 1)))
})

test_that("prepare_weather preserves grouping values and normalises names", {
  input <- data.frame(
    DATE = c("2000-01-02", "2000-01-01"), TMIN = c(2, 3), TMAX = c(6, 7),
    Scenario = c("s1", "s2"), MODEL = c("m1", "m2"),
    Period = c("historical", "future"), Station = c("04709", "00001")
  )
  result <- prepare_weather(input)
  groups <- c("scenario", "model", "period", "station")
  expected <- input[2:1, 4:7]
  names(expected) <- groups
  expect_identical(result[groups], expected)
  input$TMAX[1] <- NA_real_
  expect_error(prepare_weather(input), "cannot contain NA")
  input$TMAX[1] <- "warm"
  expect_error(prepare_weather(input), "must be numeric")
})

test_that("the observed weather example covers every day in its source year", {
  path <- system.file("extdata", "dwd_simmern_wahlbach_2005.csv", package = "ThermoPheno")
  example <- prepare_weather(read.csv(path, colClasses = c(station = "character")))
  expect_equal(example$date, seq(as.Date("2005-01-01"), as.Date("2005-12-31"), by = "day"))
  expect_identical(unique(example$station), "04709")
  expect_true(all(example$tmin <= example$tmax))
})
