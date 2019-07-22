## gtf to genePred
p1 <- InputParam(id = "gtf", type = "File", position = 1)
p2 <- InputParam(id = "gPred", type = "string", position = 2)
o1 <- OutputParam(id = "genePred", type = "File", glob = "$(inputs.gPred)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "hubentu/rcwl-rnaseq")
gtfToGenePred <- cwlParam(baseCommand = "gtfToGenePred",
                          requirements = list(req1),
                          inputs = InputParamList(p1, p2),
                          outputs = OutputParamList(o1))
