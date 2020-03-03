## https://github.com/timothyjamesbecker/SVE
p1 <- InputParam(id = "fqs", type = "File[]",
                 itemSeparator = ",", prefix = "-f")
p2 <- InputParam(id = "ref", type = "File",
                 secondaryFile = ".fai", prefix = "-r")
p3 <- InputParam(id = "outdir", type = "string", prefix = "-o")
p4 <- InputParam(id = "threads", type = "int",
                 prefix = "-T", default = 4L)
o1 <- OutputParam(id = "outs", type = "Directory", glob = "$(inputs.outdir)")

req1 <- list(class = "DockerRequirement",
             dockerPull = "timothyjamesbecker/sve")
SVE <- cwlParam(baseCommand = "/software/SVE/scripts/auto.py",
                requirements = list(req1),
                inputs = InputParamList(p1, p2, p3, p4),
                outputs = OutputParamList(o1))
