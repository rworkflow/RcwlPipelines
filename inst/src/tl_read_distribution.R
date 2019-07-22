## read distribution
p1 <- InputParam(id = "bam", type = "File", prefix = "-i")
p2 <- InputParam(id = "bed", type = "File", prefix = "-r")
o1 <- OutputParam(id = "distOut", type = "File", glob = "$(inputs.bam.nameroot).distribution.txt")
req1 <- list(class = "DockerRequirement",
             dockerPull = "hubentu/rcwl-rnaseq")
read_distribution <- cwlParam(baseCommand = c("read_distribution.py"),
                              requirements = list(req1),
                              inputs = InputParamList(p1, p2),
                              outputs = OutputParamList(o1),
                              stdout = "$(inputs.bam.nameroot).distribution.txt")
