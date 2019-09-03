## bgzip -c file > file.gz
p1 <- InputParam(id = "file", type = "File")
o1 <- OutputParam(id = "zfile", type = "File", glob = "$(inputs.file).gz")
req1 <- list(class = "DockerRequirement",
             dockerPull = "biocontainers/tabix:v1.3.2-2-deb_cv1")
bgzip <- cwlParam(baseCommand = c("bgzip", "-c"),
                  requirements = list(req1),
                  inputs = InputParamList(p1),
                  outputs = OutputParamList(o1),
                  stdout = "$(inputs.file).gz")

