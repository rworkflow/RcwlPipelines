## vatools.org

p1 <- InputParam(id = "ivcf", type = "File", position = 1)
p2 <- InputParam(id = "readcount", type = "File", position = 2)
p3 <- InputParam(id = "ntype", type = "string", position = 3, default = "DNA")
p4 <- InputParam(id = "sample", type = "string?", prefix = "-s")
p5 <- InputParam(id = "vtype", type = "string", prefix = "-t", default = "snv")
p6 <- InputParam(id = "ovcf", type = "string", prefix = "-o")
o1 <- OutputParam(id = "oVcf", type = "File", glob = "$(inputs.ovcf)")

req1 <- list(class = "DockerRequirement",
             dockerPull = "griffithlab/vatools:3.1.0")
vcf_readcount_annotator <- cwlParam(baseCommand = "vcf-readcount-annotator",
                                    requirements = list(req1),
                                    inputs = InputParamList(p1, p2, p3, p4, p5, p6),
                                    outputs = OutputParamList(o1))
