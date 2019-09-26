## CollectSequencingArtifactMetrics
p1 <- InputParam(id = "bam", type = "File", prefix = "-I")
p2 <- InputParam(id = "ref", type = "File", prefix = "-R")
p3 <- InputParam(id = "ext", type = "string", prefix = "--FILE_EXTENSION", default = ".txt")
p4 <- InputParam(id = "art", type = "string", prefix = "-O")
o1 <- OutputParam(id = "aout", type = "File",
                  glob = "$(inputs.art).pre_adapter_detail_metrics.txt")
req1 <- list(class = "DockerRequirement",
             dockerPull = "broadinstitute/gatk:latest")

ColSeqArtifact <- cwlParam(baseCommand = c("gatk",
                                           "CollectSequencingArtifactMetrics"),
                           requirements = list(req1),
                           inputs = InputParamList(p1, p2, p3, p4),
                           outputs = OutputParamList(o1))
