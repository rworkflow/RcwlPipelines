## GetPileupSummaries
p1 <- InputParam(id = "bam", type = "File", prefix = "-I")
p2 <- InputParam(id = "vcf", type = "File", prefix = "-V", secondaryFiles = ".idx")
p3 <- InputParam(id = "interval", type = "File", prefix = "-L", secondaryFiles = ".idx")
p4 <- InputParam(id = "pileup", type = "string", prefix = "-O")
o1 <- OutputParam(id = "pout", type = "File", glob = "$(inputs.pileup)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "broadinstitute/gatk:latest")

GetPileupSummaries <- cwlParam(baseCommand = c("gatk", "GetPileupSummaries"),
                               requirements = list(req1),
                               inputs = InputParamList(p1, p2, p3, p4),
                               outputs = OutputParamList(o1))
