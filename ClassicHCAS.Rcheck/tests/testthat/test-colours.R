test_that("palettes returns the expected palette variants", {
  expect_equal(palettes(4), palettes(4, "hcas"))
  expect_length(palettes(7, "hcas"), 7)
  expect_length(palettes(7, "ref_density"), 7)
  expect_false(identical(palettes(4, "hcas"), palettes(4, "ref_density")))
})

test_that("palettes rejects unknown palette names", {
  expect_error(palettes(4, "unknown"), "should be one of")
})
