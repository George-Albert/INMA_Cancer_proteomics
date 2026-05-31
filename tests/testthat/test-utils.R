test_that("build_threshold_tag generates safe names", {
  source(file.path("..", "..", "src", "utils.R"))
  expect_equal(build_threshold_tag(0.2), "mean_gt_0.2")
  expect_equal(build_threshold_tag(1), "mean_gt_1")
})
