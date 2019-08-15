## samtools stats
p1 <- InputParam(id = "bam", type = "File")
o1 <- OutputParam(id = "stats", type = "File", glob = "$(inputs.bam.nameroot).stats.txt")
req2 <- list(class = "DockerRequirement",
             dockerPull = "biocontainers/samtools:v1.7.0_cv3")
samtools_stats <- cwlParam(baseCommand = c("samtools", "stats"),
                           requirements = list(req2),
                           inputs = InputParamList(p1),
                           outputs = OutputParamList(o1),
                           stdout = "$(inputs.bam.nameroot).stats.txt")
