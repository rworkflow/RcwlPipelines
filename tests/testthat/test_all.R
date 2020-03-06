
tools <- cwlUpdate(tempdir())
test_that("check source tool scripts", {
    expect_is(tools, "BiocFileCache")
})

test_that("check tool type", {
    expect_equal(cwlSearch("tl_bcfview", tools)$Type, "tool")
})

test_that("check pipeline type", {
    expect_equal(cwlSearch("pl_neusomatic", tools)$Type, "pipeline")
})
