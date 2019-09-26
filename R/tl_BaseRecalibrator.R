## BaseRecalibrator
p1 <- InputParam(id = "bam", type = "File", prefix = "-I")
p2 <- InputParam(id = "ref", prefix = "-R", type = "File", secondaryFiles = c(
                                                               ".fai",
                                                               "$(self.nameroot).dict"))
p3 <- InputParam(id = "knowSites", type = InputArrayParam(items = "File",
                                                          prefix = "--known-sites"),
                 secondaryFiles = ".idx")
p4 <- InputParam(id = "recal", type = "string", prefix = "-O")
o1 <- OutputParam(id = "rtable", type = "File",
                  glob = "$(inputs.recal)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "broadinstitute/gatk:latest")

BaseRecalibrator <- cwlParam(baseCommand = c("gatk",
                                             "BaseRecalibrator"),
                             requirements = list(req1),
                             inputs = InputParamList(p1, p2, p3, p4),
                             outputs = OutputParamList(o1))
