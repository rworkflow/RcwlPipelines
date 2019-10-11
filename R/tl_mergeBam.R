## merge BAMs
p1 <- InputParam(id = "ibam", type = InputArrayParam(items = "File",
                                                    prefix = "I=",
                                                    separate = FALSE))
p2 <- InputParam(id = "obam", type = "string", prefix = "O=", separate = FALSE)
o1 <- OutputParam(id = "oBam", type = "File", glob = "$(inputs.obam)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "broadinstitute/picard")
mergeBam <- cwlParam(baseCommand = c("java", "-jar", "/usr/picard/picard.jar",
                                     "MergeSamFiles"),
                     requirements = list(req1),
                     inputs = InputParamList(p1, p2),
                     outputs = OutputParamList(o1))
