## CalculateContamination
p1 <- InputParam(id = "ttable", type = "File", prefix = "-I")
p2 <- InputParam(id = "ntable", type = "File", prefix = "-matched")
p3 <- InputParam(id = "cont", type = "string", prefix = "-O")
o1 <- OutputParam(id = "cout", type = "File", glob = "$(inputs.cont)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "broadinstitute/gatk:latest")

CalculateContamination <- cwlParam(baseCommand = c("gatk", "CalculateContamination"),
                                   requirements = list(req1),
                                   inputs = InputParamList(p1, p2, p3),
                                   outputs = OutputParamList(o1))
