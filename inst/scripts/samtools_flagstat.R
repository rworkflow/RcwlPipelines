## bam flagstat
## Note: result will be generated to current dir.
p1 <- InputParam(id = "bam", type = "File")
o1 <- OutputParam(id = "flagstat", type = "File", glob = "$(inputs.bam.nameroot).flagstat.txt")

req1 <- list(class = "DockerRequirement",
             dockerPull = "hubentu/rcwl-rnaseq")
samtools_flagstat <- cwlParam(baseCommand = c("samtools", "flagstat"),
                              requirements = list(req1),
                              inputs = InputParamList(p1),
                              outputs = OutputParamList(o1),
                              stdout = "$(inputs.bam.nameroot).flagstat.txt")
