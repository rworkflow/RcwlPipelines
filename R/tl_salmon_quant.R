## salmon_quant
## Input: fq or fq.gz

rep1 <- list(class = "DockerRequirement",
             dockerPull = "combinelab/salmon")
rep2 <- list(class = "InlineJavascriptRequirement")

p1 <- InputParam(id = "threadN", type = "int", prefix = "-p", position = 1)
p2 <- InputParam(id = "ref", type = "Directory", prefix = "-i", position = 2)
p3 <- InputParam(id = "fq1", type = "File", prefix = "-1", position = 3)
p4 <- InputParam(id = "fq2", type = "File", prefix = "-2", position = 4)
p5 <- InputParam(id = "outPrefix", type = "string", prefix = "-o", position = 5)
o1 <- OutputParam(id = "out1", type = "Directory", glob = "$(inputs.outPrefix)")

salmon_quant <- cwlParam(baseCommand = c("salmon", "quant"), 
                         requirements = list(rep1, rep2),
                         arguments = list("-l", "A",
                                          "--validateMappings",
                                          "--seqBias",
                                          "--gcBias",
                                          "--posBias"),
                         inputs = InputParamList(p1, p2, p3, p4, p5), 
                         outputs = OutputParamList(o1))
