## picard markdup
p1 <- InputParam(id = "ibam", type = "File", prefix = "I=", separate = FALSE)
p2 <- InputParam(id = "obam", type = "string", prefix = "O=", separate = FALSE)
p3 <- InputParam(id = "matrix", type = "string", prefix = "M=", separate = FALSE)
o1 <- OutputParam(id = "mBam", type = "File", glob = "$(inputs.obam)")
o2 <- OutputParam(id = "Mat", type = "File", glob = "$(inputs.matrix)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "broadinstitute/picard")
markdup <- cwlParam(baseCommand = c("java", "-jar", "/usr/picard/picard.jar",
                                    "MarkDuplicates"),
                    requirements = list(req1),
                    inputs = InputParamList(p1, p2, p3),
                    outputs = OutputParamList(o1, o2))
