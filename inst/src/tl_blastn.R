## blast:v2.2.31_cv2

p1 <- InputParam(id = "ThreadN", type = "int", prefix = "-num_threads")
p2 <- InputParam(id = "Ref", type = "File", prefix = "-db", secondaryFiles = c(".nhr", ".nin", ".nsq"))
p3 <- InputParam(id = "Query", type = "File", prefix = "-query")
p4 <- InputParam(id = "IdenPerc", type = "int", prefix = "-perc_identity")
p5 <- InputParam(id = "WordSize", type = "int", prefix = "-word_size")
p6 <- InputParam(id = "Evalue", type = "float", prefix = "-evalue")
p7 <- InputParam(id = "OutFormat", type = "int", prefix = "-outfmt")
p8 <- InputParam(id = "OutFile", type = "string", prefix = "-out")
o1 <- OutputParam(id = "Output", type = "File", glob = "$(inputs.OutFile)")

req1 <- list(class = "DockerRequirement", 
             dockerPull = "biocontainers/blast:v2.2.31_cv2")
blastn <- cwlParam(baseCommand = "blastn",
                   requirements = list(req1),
                   inputs = InputParamList(p1, p2, p3, p4, p5, p6, p7, p8),
                   outputs = OutputParamList(o1))
