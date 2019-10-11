
p1 <- InputParam(id = "vcf", type = "File", prefix = "--variant")
p2 <- InputParam(id = "bam", type = "File", prefix = "-I", secondaryFiles = ".bai")
p3 <- InputParam(id = "ref", type = "File", prefix = "-R",
                 secondaryFiles = c(".fai", "$(self.nameroot).dict"))
p4 <- InputParam(id = "ovcf", type = "string", prefix = "-o")
p5 <- InputParam(id = "region", type = "File", prefix = "-L")
o1 <- OutputParam(id = "oVcf", type = "File", glob = "$(inputs.ovcf)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "broadinstitute/gatk3:3.8-1")
ReadBackedPhasing <- cwlParam(baseCommand = c("java", "-jar", "/usr/GenomeAnalysisTK.jar",
                                              "-T", "ReadBackedPhasing"),
                              requirements = list(req1),
                              inputs = InputParamList(p1, p2, p3, p4, p5),
                              outputs = OutputParamList(o1))
