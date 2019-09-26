
p1 <- InputParam(id = "bed", type = "File", prefix = "I=", separate = FALSE)
p2 <- InputParam(id = "SD", type = "File", prefix = "SD=", separate = FALSE)
p3 <- InputParam(id = "out", type = "string", prefix = "O=", separate = FALSE)
o1 <- OutputParam(id = "intval", type = "File", glob = "$(inputs.out)")

req1 <- list(class = "DockerRequirement",
             dockerPull = "biocontainers/picard:2.3.0")
BedToIntervalList <- cwlParam(baseCommand = c("picard", "BedToIntervalList"),
                              requirements = list(req1),
                              inputs = InputParamList(p1, p2, p3),
                              outputs = OutputParamList(o1))
