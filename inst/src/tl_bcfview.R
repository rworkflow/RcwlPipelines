
p1 <- InputParam(id = "vcf", type = "File")
p2 <- InputParam(id = "filter", type = "string", prefix = "-f", default = "PASS")
p3 <- InputParam(id = "fout", type = "string", prefix = "-o")
o1 <- OutputParam(id = "Fout", type = "File", glob = "$(inputs.fout)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "biocontainers/bcftools:v1.5_cv3")

bcfview <- cwlParam(baseCommand = c("bcftools", "view"),
                    requirements = list(req1),
                    inputs = InputParamList(p1, p2, p3),
                    outputs = OutputParamList(o1))
