
tools <- cwlUpdate(tempdir())
test_that("check source tool scripts", {
    expect_is(tools, "cwlHub")
})

test_that("check tool type", {
    expect_equal(Type(cwlSearch("tl_bcfview", tools)), "tool")
})

test_that("check pipeline type", {
    expect_equal(Type(cwlSearch("pl_neusomatic", tools)), "pipeline")
})

test_that("check cwlLoad", {
    bwaMRecal <- cwlLoad(title(cwlSearch("pl_bwaMRecal", tools)), tools)
    expect_true(exists("bwaMRecal"))
    expect_true(exists("bwa"))
    expect_true(exists("BaseRecal"))
})
