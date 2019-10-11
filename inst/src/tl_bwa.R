## bwa mem
p1 <- InputParam(id = "threads", type = "int", prefix = "-t", position = 1)
p2 <- InputParam(id = "RG", type = "string", prefix = "-R", position = 2)
p3 <- InputParam(id = "Ref", type = "File", position = 3,
                 secondaryFiles = c(".amb", ".ann", ".bwt", ".pac", ".sa"))
p4 <- InputParam(id = "FQ1", type = "File", position = 4)
p5 <- InputParam(id = "FQ2", type = "File?", position = 5)
o1 <- OutputParam(id = "sam", type = "File", glob = "*.sam")
req1 <- list(class = "DockerRequirement",
             dockerPull = "biocontainers/bwa:v0.7.17-3-deb_cv1")
bwa <- cwlParam(baseCommand = c("bwa", "mem"),
                requirements = list(req1),
                inputs = InputParamList(p1, p2, p3, p4, p5),
                outputs = OutputParamList(o1),
                stdout = "bwaOutput.sam")
