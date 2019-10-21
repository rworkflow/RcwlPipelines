## merge BAMs
p1 <- InputParam(id = "ibam", type = InputArrayParam(items = "File",
                                                    prefix = "I=",
                                                    separate = FALSE))
p2 <- InputParam(id = "obam", type = "string", prefix = "O=", separate = FALSE)
o1 <- OutputParam(id = "oBam", type = "File", glob = "$(inputs.obam)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "quay.io/biocontainers/picard:2.21.1--0")
mergeBam <- cwlParam(baseCommand = c("picard",
                                     "MergeSamFiles"),
                     requirements = list(req1),
                     inputs = InputParamList(p1, p2),
                     outputs = OutputParamList(o1))
