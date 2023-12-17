test_that("coordinates", {
  data("visium_data")
  coords <- vrCoordinates(visium_data)
  coords <- vrCoordinates(visium_data, reg = TRUE)
  segments <- vrSegments(visium_data)
  segments <- vrSegments(visium_data, reg = TRUE)
  expect_equal(1,1L)
})
