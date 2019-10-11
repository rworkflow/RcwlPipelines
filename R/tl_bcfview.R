
p1 <- InputParam(id = "vcf", type = "File")
p2 <- InputParam(id = "filter", type = "string?", prefix = "-f", default = "PASS")
p3 <- InputParam(id = "fout", type = "string", prefix = "-o")
p4 <- InputParam(id = "otype", type = "string?", prefix = "-O", default = "v")
p5 <- InputParam(id = "sample", type = "string?", prefix = "-s")
p6 <- InputParam(id = "samplefile", type = "File?", prefix = "-S")
p7 <- InputParam(id = "genotype", type = "string?", prefix = "-g")
p8 <- InputParam(id = "include", type = "string?", prefix = "-i")
p9 <- InputParam(id = "exclude", type = "string?", prefix = "-e")

o1 <- OutputParam(id = "Fout", type = "File", glob = "$(inputs.fout)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "biocontainers/bcftools:v1.5_cv3")

bcfview <- cwlParam(baseCommand = c("bcftools", "view"),
                    requirements = list(req1),
                    inputs = InputParamList(p1, p2, p3, p4, p5, p6, p7, p8, p9),
                    outputs = OutputParamList(o1))
