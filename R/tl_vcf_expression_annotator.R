## vatools.org

p1 <- InputParam(id = "ivcf", type = "File", position = 1)
p2 <- InputParam(id = "expression", type = "File", position = 2)
p3 <- InputParam(id = "etype", type = "string", position = 3, default = "kallisto")
p4 <- InputParam(id = "gtype", type = "string", position = 4, default = "transcript")
p5 <- InputParam(id = "idCol", type = "string?", prefix = "-i")
p6 <- InputParam(id = "expCol", type = "string?", prefix = "-e")
p7 <- InputParam(id = "sample", type = "string?", prefix = "-s")
p8 <- InputParam(id = "ovcf", type = "string", prefix = "-o")
o1 <- OutputParam(id = "oVcf", type = "File", glob = "$(inputs.ovcf)")

req1 <- list(class = "DockerRequirement",
             dockerPull = "griffithlab/vatools:3.1.0")
vcf_expression_annotator <- cwlParam(baseCommand = "vcf-expression-annotator",
                                     requirements = list(req1),
                                     arguments = list("--ignore-transcript-version"),
                                     inputs = InputParamList(p1, p2, p3, p4, p5, p6, p7, p8),
                                     outputs = OutputParamList(o1))
