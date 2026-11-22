test_that("simple thermal time is correct", {
  expect_equal(calc_daily_tt(tmin = 5, tmax = 15, t_base = 5, mode = "simple"), 5)
  expect_equal(calc_daily_tt(tmin = 0, tmax = 4, t_base = 5, mode = "simple"), 0)
})

test_that("capped thermal time is correct", {
  expect_equal(calc_daily_tt(tmin = 20, tmax = 40, t_base = 5, t_opt = 25, mode = "capped"), 20)
})

test_that("triangular thermal time is correct", {
  expect_equal(calc_daily_tt(tmin = 20, tmax = 40, t_base = 5, t_opt = 25, t_max_cut = 35, mode = "triangular"), 10)
  expect_error(calc_daily_tt(tmin = 20, tmax = 40, t_base = 5, mode = "triangular"))
})

test_that("all methods handle thresholds and subzero weather", {
  temperatures <- c(-4, 5, 25, 35)
  expect_equal(calc_daily_tt(temperatures, temperatures, 5, 25, 35, "simple"), c(0, 0, 20, 30))
  expect_equal(calc_daily_tt(temperatures, temperatures, 5, 25, 35, "capped"), c(0, 0, 20, 20))
  expect_equal(calc_daily_tt(temperatures, temperatures, 5, 25, 35, "triangular"), c(0, 0, 20, 0))
})

test_that("reference calendar thermal time uses the requested baseline years", {
  weather <- data.frame(
    date = as.Date(c("2000-03-01", "2000-03-02", "2000-03-03",
                     "2001-03-01", "2001-03-02", "2001-03-03")),
    tmin = c(0, 4, 8, 10, 10, 10),
    tmax = c(4, 8, 12, 20, 20, 20)
  )
  result <- estimate_required_tt(weather, 2000L, "03-01", 3, t_base = 5)
  expect_equal(result$required_tt, 6)
  expect_equal(result$yearly_required_tt$year, 2000)
})
