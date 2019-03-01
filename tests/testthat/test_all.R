
test_that("total pipelines number", {
    expect_equal(nrow(data(package="RcwlPipelines")$results),
                 7)})
