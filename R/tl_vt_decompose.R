## https://genome.sph.umich.edu/wiki/Vt

p1 <- InputParam(id = "ivcf", type = "File", position = 1)
p2 <- InputParam(id = "ovcf", type = "string", prefix = "-o", position = 2)
o1 <- OutputParam(id = "oVcf", type = "File", glob = "$(inputs.ovcf)")
    
req1 <- list(class = "DockerRequirement",
             dockerPull = "hubentu/vt")
vt_decompose <- cwlParam(baseCommand = c("vt", "decompose"),
                         requirements = list(req1),
                         arguments = list("-s"),
                         inputs = InputParamList(p1, p2),
                         outputs = OutputParamList(o1))
