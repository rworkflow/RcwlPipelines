## hisat2:v2.0.5-1-deb_cv1
## The ref does not exceed 4 billion characters
## Basename of "ref" should be exactly same with Name of "outPrefix" 

req1 <- list(class = "DockerRequirement",
             dockerPull = "biocontainers/hisat2:v2.0.5-1-deb_cv1")
p1 <- InputParam(id = "ref", type = "File", position = 1)
p2 <- InputParam(id = "outPrefix", type = "string", position = 2)
o1 <- OutputParam(id = "idx", type = "File[]", 
                  glob = "$(inputs.outPrefix).*")
hisat2_build <- cwlParam(baseCommand = c("hisat2-build"), 
                         requirements = list(req1),
                         inputs = InputParamList(p1, p2), 
                         outputs = OutputParamList(o1))
