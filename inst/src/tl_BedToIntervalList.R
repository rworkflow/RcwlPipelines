
p1 <- InputParam(id = "bed", type = "File", prefix = "I=", separate = FALSE)
p2 <- InputParam(id = "SD", type = "File", prefix = "SD=", separate = FALSE)
p3 <- InputParam(id = "out", type = "string", prefix = "O=", separate = FALSE)
o1 <- OutputParam(id = "intval", type = "File", glob = "$(inputs.out)")

req1 <- list(class = "DockerRequirement",
             dockerPull = "broadinstitute/picard")
BedToIntervalList <- cwlParam(baseCommand = c("java", "-jar", "/usr/picard/picard.jar",
                                              "BedToIntervalList"),
                              requirements = list(req1),
                              inputs = InputParamList(p1, p2, p3),
                              outputs = OutputParamList(o1))
