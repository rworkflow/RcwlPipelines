## https://pvactools.readthedocs.io/en/latest/pvacseq/getting_started.html#running-pvacseq-using-docker

p1 <- InputParam(id = "ivcf", type = "File", position = 1, secondaryFiles = ".tbi")
p2 <- InputParam(id = "sample", type = "string", position = 2)
p3 <- InputParam(id = "allele", type = "string[]",
                 itemSeparator = ",", position = 3)
p4 <- InputParam(id = "algorithms", type = "string[]", position = 4)
p5 <- InputParam(id = "outdir", type = "string", position = 5, default = "pvacseq_out")
p6 <- InputParam(id = "length", type = "string",
                 position = 6, prefix = "-e", default = "8,9,10,11")
p7 <- InputParam(id = "phasedVcf", type = "File?",
                 position = 7, prefix = "-p", secondaryFiles = ".tbi")
o1 <- OutputParam(id = "Out", type = "Directory", glob = "$(inputs.outdir)")

req1 <- list(class = "DockerRequirement",
             dockerPull = "griffithlab/pvactools")
pvacseq <- cwlParam(baseCommand = c("pvacseq", "run"),
                    requirements = list(req1),
                    inputs = InputParamList(p1, p2, p3, p4, p5, p6, p7),
                    outputs = OutputParamList(o1))
