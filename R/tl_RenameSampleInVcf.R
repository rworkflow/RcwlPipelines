
p1 <- InputParam(id = "vcf", type = "File", prefix = "I=", separate = FALSE)
p2 <- InputParam(id = "ovcf", type = "string", prefix = "O=", separate = FALSE)
p3 <- InputParam(id = "NewName", type = "string", prefix = "NEW_SAMPLE_NAME=", separate = FALSE)
o1 <- OutputParam(id = "oVcf", type = "File", glob = "$(inputs.ovcf)")

req1 <- list(class = "DockerRequirement",
             dockerPull = "quay.io/biocontainers/picard:2.21.1--0")
RenameSampleInVcf <- cwlParam(baseCommand = c("picard",
                                              "RenameSampleInVcf"),
                              requirements = list(req1),
                              inputs = InputParamList(p1, p2, p3),
                              outputs = OutputParamList(o1))
