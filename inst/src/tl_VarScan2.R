## https://github.com/dkoboldt/varscan
p1 <- InputParam(id = "nbam", type = "File",
                 position = 1, secondaryFiles = ".bai")
p2 <- InputParam(id = "tbam", type = "File",
                 position = 2, secondaryFiles = ".bai")
p3 <- InputParam(id = "ref", type = "File", prefix = "-b",
                 secondaryFiles = c(".fai", "$(self.nameroot).dict"))
p4 <- InputParam(id = "allvcf", type = "string", prefix = "-f")
p5 <- InputParam(id = "somvcf", type = "string", prefix = "-y")
p6 <- InputParam(id = "proc", type = "int", prefix = "-p")
o1 <- OutputParam(id = "allVcf", type = "File", glob = "$(inputs.allvcf)")
o2 <- OutputParam(id = "somVcf", type = "File", glob = "$(inputs.somvcf)")

req1 <- list(class = "DockerRequirement",
             dockerPull = "serge2016/varscan:v0.1.1")
VarScan2 <- cwlParam(baseCommand = "",
                     requirements = list(req1),
                     arguments = list("-o", "."),
                     inputs = InputParamList(p1, p2, p3, p4, p5, p6),
                     outputs = OutputParamList(o1, o2))
