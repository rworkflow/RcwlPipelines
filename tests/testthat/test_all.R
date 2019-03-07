
test_that("total pipelines number", {
    expect_equal(nrow(data(package="RcwlPipelines")$results),
                 7)})

test_that("check source tool scripts", {
    tools <- cwlTools(tempdir())
    fpath <- bfcinfo(tools)$fpath
    expect_equal(ncol(sapply(fpath, source)),
    length(fpath))
})
