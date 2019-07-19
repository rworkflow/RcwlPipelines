## genePred to bed
p1 <- InputParam(id = "genePred", type = "File", position = 1)
p2 <- InputParam(id = "Bed", type = "string", position = 2)
o1 <- OutputParam(id = "bed", type = "File", glob = "$(inputs.Bed)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "hubentu/rcwl-rnaseq")
genePredToBed <- cwlParam(baseCommand = "genePredToBed",
                          requirements = list(req1),
                          inputs = InputParamList(p1, p2),
                          outputs = OutputParamList(o1))
