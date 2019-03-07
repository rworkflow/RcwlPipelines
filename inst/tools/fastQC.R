## fastqc
f1 <- InputParam(id = "seqfile", type = "File")
o1 <- OutputParam(id = "QCfile", type = "File", glob = "*.zip")

req1 <- list(class = "DockerRequirement",
             dockerPull = "hubentu/rcwl-rnaseq")
fastqc <- cwlParam(baseCommand = "fastqc",
                   requirements = list(req1),
                   arguments = list("--outdir", "./"),
                   inputs = InputParamList(f1),
                   outputs = OutputParamList(o1))

