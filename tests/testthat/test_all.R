
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

test_that("check cwlLoad", {
    bwaMRecal <- cwlLoad(cwlSearch("pl_bwaMRecal")$rname)
    expect_true(exists("bwaMRecal"))
    expect_true(exists("bwa"))
    expect_true(exists("BaseRecal"))
})
