
p1 <- InputParam(id = "bam", type = "File", prefix = "I=", separate = FALSE)
p2 <- InputParam(id = "fq1", type = "string", prefix = "F=", separate = FALSE)
p3 <- InputParam(id = "fq2", type = "string", prefix = "F2=", separate = FALSE)
o1 <- OutputParam(id = "FQ1", type = "File", glob = "$(inputs.fq1)")
o2 <- OutputParam(id = "FQ2", type = "File", glob = "$(inputs.fq2)")

req1 <- list(class = "DockerRequirement",
             dockerPull = "quay.io/biocontainers/picard:2.21.1--0")
SamToFastq <- cwlParam(baseCommand = c("picard",
                                       "SamToFastq"),
                       requirements = list(req1),
                       inputs = InputParamList(p1, p2, p3),
                       outputs = OutputParamList(o1, o2))
