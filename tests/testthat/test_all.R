
## mc3 removed
test_that("total pipelines number", {
    expect_equal(nrow(data(package="RcwlPipelines")$results),
                 6)
})

test_that("check source tool scripts", {
    tools <- cwlTools(tempdir())
    expect_is(tools, "BiocFileCache")
})
