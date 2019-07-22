## RUN WDL
p1 <- InputParam(id = "cromwell", type = "File", position = 1, prefix = "-jar")
p2 <- InputParam(id = "run", type = "string", default = "run", position = 2)
p3 <- InputParam(id = "wdl", type = "File", position = 3)
p4 <- InputParam(id = "json", type = "File", position = 4, prefix = "-i")
o1 <- OutputParam(id = "log", type = "File", glob="$(inputs.wdl.basename).log")
runWDL <- cwlParam(baseCommand = "java",
                   inputs = InputParamList(p1, p2, p3, p4),
                   outputs = OutputParamList(o1),
                   stdout = "$(inputs.wdl.basename).log")
