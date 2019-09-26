## BaseRecalibrator
p1 <- InputParam(id = "bam", type = "File", prefix = "-I")
p2 <- InputParam(id = "ref", prefix = "-R", type = "File", secondaryFiles = c(".fai", "$(self.nameroot).dict"))
p3 <- InputParam(id = "rtable", type = "File", prefix = "--bqsr-recal-file")
p4 <- InputParam(id = "oBam", type = "string", prefix = "-O")
o1 <- OutputParam(id = "Bam", type = "File",
                  glob = "$(inputs.oBam)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "broadinstitute/gatk:latest")

ApplyBQSR <- cwlParam(baseCommand = c("gatk",
                                      "ApplyBQSR"),
                      requirements = list(req1),
                      inputs = InputParamList(p1, p2, p3, p4),
                      outputs = OutputParamList(o1))
