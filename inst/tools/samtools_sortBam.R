## samtools sort bam
p1 <- InputParam(id = "bam", type = "File")
o1 <- OutputParam(id = "sbam", type = "File", glob = "$(inputs.bam.nameroot).sorted.bam")
req2 <- list(class = "DockerRequirement",
             dockerPull = "biocontainers/samtools:v1.7.0_cv3")
sortBam <- cwlParam(baseCommand = c("samtools", "sort"),
                    requirements = list(req2),
                    inputs = InputParamList(p1),
                    outputs = OutputParamList(o1),
                    stdout = "$(inputs.bam.nameroot).sorted.bam")
