
p1 <- InputParam(id = "tfile", type = "File", position = 1)
o1 <- OutputParam(id = "idx", type = "File", glob = "$(inputs.tfile.basename)",
                  secondaryFiles = list(".tbi"))
req1 <- list(class = "DockerRequirement",
             dockerPull = "biocontainers/tabix:v1.3.2-2-deb_cv1")
req2 <- list(class = "InitialWorkDirRequirement",
             listing = list("$(inputs.tfile)"))
tabix_index <- cwlParam(baseCommand = "tabix",
                        requirements = list(req1, req2),
                        inputs = InputParamList(p1),
                        outputs = OutputParamList(o1))
