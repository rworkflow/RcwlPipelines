
p1 <- InputParam(id = "bam", type = "File",
                 prefix = "-b", secondaryFiles = ".bai")
p2 <- InputParam(id = "bw", type = "string", prefix = "-o")
p3 <- InputParam(id = "binsize", type = "int",
                 prefix = "-bs", default = 1L)
p4 <- InputParam(id = "processors", type = "int",
                 prefix = "-p", default = 2L)
o1 <- OutputParam(id = "bigwig", type = "File", glob = "$(inputs.bw)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "quay.io/biocontainers/deeptools:3.3.1--py_0")
bamCoverage <- cwlParam(baseCommand = "bamCoverage",
                        requirements = list(req1),
                        arguments = list("--ignoreDuplicates"),
                        inputs = InputParamList(p1, p2, p3, p4),
                        outputs = OutputParamList(o1))
