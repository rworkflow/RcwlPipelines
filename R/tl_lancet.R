## https://github.com/nygenome/lancet

p1 <- InputParam(id = "tbam", type = "File", prefix = "--tumor",
                 secondaryFiles = ".bai")
p2 <- InputParam(id = "nbam", type = "File", prefix = "--normal",
                 secondaryFiles = ".bai")
p3 <- InputParam(id = "ref", type = "File",
                 prefix = "--ref", secondaryFiles = ".fai")
p4 <- InputParam(id = "region", type = "File", prefix = "--bed")
p5 <- InputParam(id = "threads", type = "int", prefix = "--num-threads")
o1 <- OutputParam(id = "vcf", type = "stdout")

req1 <- list(class = "DockerRequirement",
             dockerPull = "kfdrc/lancet:1.0.7")
lancet <- cwlParam(baseCommand = "/lancet-1.0.7/lancet",
                   requirements = list(req1),
                   inputs = InputParamList(p1, p2, p3, p4, p5),
                   outputs = OutputParamList(o1),
                   stdout = "$(inputs.tbam.nameroot)_$(inputs.nbam.nameroot).vcf")
