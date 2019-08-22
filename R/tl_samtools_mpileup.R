
p1 <- InputParam(id = "bam", type = "File")
p2 <- InputParam(id = "ref", type = "File",
                 secondaryFiles = ".fai", prefix = "-f")
p3 <- InputParam(id = "region", type = "File?", prefix = "-l")
o1 <- OutputParam(id = "pileup", type = "File", glob = "$(inputs.bam.nameroot).pileup")
req1 <- list(class = "DockerRequirement",
             dockerPull = "biocontainers/samtools:v1.7.0_cv3")
samtools_mpileup <- cwlParam(baseCommand = c("samtools", "mpileup"),
                             requirements = list(req1),
                             inputs = InputParamList(p1, p2, p3),
                             outputs = OutputParamList(o1),
                             stdout = "$(inputs.bam.nameroot).pileup")
