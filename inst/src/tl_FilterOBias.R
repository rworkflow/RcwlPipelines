## FilterByOrientationBias
p1 <- InputParam(id = "vcf", type = "File", prefix = "-V")
p2 <- InputParam(id = "art", type = "File", prefix = "-P")
p3 <- InputParam(id = "mode", type = InputArrayParam(items = "string", prefix = "--artifact-modes"))
p4 <- InputParam(id = "avcf", type = "string", prefix = "-O")
o1 <- OutputParam(id = "fout", type = "File", glob = "$(inputs.avcf)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "broadinstitute/gatk:latest")

FilterOBias <- cwlParam(baseCommand = c("gatk",
                                        "FilterByOrientationBias"),
                        requirements = list(req1),
                        inputs = InputParamList(p1, p2, p3, p4),
                        outputs = OutputParamList(o1))
