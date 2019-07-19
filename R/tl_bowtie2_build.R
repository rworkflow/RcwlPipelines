## bowtie2:v2.2.9_cv2
## Basename of "ref" should be exactly same with Name of "outPrefix"

req1 <- list(class = "DockerRequirement", 
             dockerPull = "biocontainers/bowtie2:v2.2.9_cv2")
p1 <- InputParam(id = "ref", type = "File", position = 1)
p2 <- InputParam(id = "outPrefix", type = "string", position = 2)
o1 <- OutputParam(id = "idx", type = "File[]", 
                  glob = "$(inputs.outPrefix).*")
bowtie2_build <- cwlParam(baseCommand = c("bowtie2-build"), 
                         requirements = list(req1),
                         inputs = InputParamList(p1, p2), 
                         outputs = OutputParamList(o1))
