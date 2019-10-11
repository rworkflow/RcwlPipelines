## bowtie:v1.2.2dfsg

req1 <- list(class = "DockerRequirement", 
             dockerPull = "biocontainers/bowtie:v1.2.2dfsg-4-deb_cv1")
p1 <- InputParam(id = "ref", type = "File", position = 1)
p2 <- InputParam(id = "outPrefix", type = "string", position = 2)
o1 <- OutputParam(id = "idx", type = "File[]",
                  glob = "$(inputs.outPrefix).*")
bowtie_build <- cwlParam(baseCommand = c("bowtie-build"), 
                         requirements = list(req1),
                         inputs = InputParamList(p1, p2), 
                         outputs = OutputParamList(o1))
