## mutliqc
p1 <- InputParam(id = "dir", type = "Directory")
o1 <- OutputParam(id = "qc", type = "File", glob = "*.html")
o2 <- OutputParam(id = "qcDat", type = "Directory", glob = "multiqc_data")
req1 <- list(class = "DockerRequirement",
             dockerPull = "hubentu/rcwl-rnaseq")

multiqc <- cwlParam(baseCommand = "multiqc",
                    requirements = list(req1),
                    inputs = InputParamList(p1),
                    outputs = OutputParamList(o1, o2))
