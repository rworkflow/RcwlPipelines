## gatk Mutect2

p1a <- InputParam(id = "tbam", type = "File",
                  prefix = "-I", secondaryFiles = ".bai")
p1b <- InputParam(id = "nbam", type = "File?",
                  prefix = "-I", secondaryFiles = ".bai")
p2 <- InputParam(id = "Ref", prefix = "-R", type = "File", secondaryFiles = c(".fai", "$(self.nameroot).dict"))
p3 <- InputParam(id = "normal", type = "string?", prefix = "-normal")
p4 <- InputParam(id = "germline", type = "File?", prefix = "--germline-resource", secondaryFiles = ".idx")
p5 <- InputParam(id = "pon", type = "File?", prefix = "--panel-of-normals", secondaryFiles = ".idx")
p6 <- InputParam(id = "interval", type = "File?", prefix = "-L")
p7 <- InputParam(id = "out", type = "string", prefix = "-O")
o1 <- OutputParam(id = "vout", type = "File", glob = "$(inputs.out)", secondaryFiles = c(".idx", ".stats"))
req1 <- list(class = "DockerRequirement",
             dockerPull = "broadinstitute/gatk:4.1.3.0")

Mutect2 <- cwlParam(baseCommand = c("gatk", "Mutect2"),
                    requirements = list(req1),
                    inputs = InputParamList(p1a, p1b, p2, p3, p4, p5, p6, p7),
                    outputs = OutputParamList(o1))
