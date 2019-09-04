
tools <- cwlTools(tempdir())
## mc3 removed
test_that("check source tool scripts", {
    expect_is(tools, "BiocFileCache")
})

test_that("check tool type", {
    expect_equal(bfcquery(tools, "bcfview")$Type, "tool")
})

test_that("check pipeline type", {
    expect_equal(bfcquery(tools, "neusomatic$")$Type, "pipeline")
})
