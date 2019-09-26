
p1 <- InputParam(id = "vcf", type = "File", prefix = "I=", separate = FALSE)
p2 <- InputParam(id = "ovcf", type = "string", prefix = "O=", separate = FALSE)
p3 <- InputParam(id = "dict", type = "File?", prefix = "SD=", separate = FALSE)
o1 <- OutputParam(id = "oVcf", type = "File", glob = "$(inputs.ovcf)")

req1 <- list(class = "DockerRequirement",
             dockerPull = "biocontainers/picard:2.3.0")
SortVcf <- cwlParam(baseCommand = c("picard", "SortVcf"),
                    requirements = list(req1),
                    inputs = InputParamList(p1, p2, p3),
                    outputs = OutputParamList(o1))
