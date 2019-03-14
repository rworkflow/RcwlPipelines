## bowtie2_2.3.0-2-deb

p1 <- InputParam(id = "threads", type = "int", prefix = "-p", position = 1)
p2 <- InputParam(id = "Ref Index prefix", type = "string", prefix = "-x", position = 2)
p3 <- InputParam(id = "FQ1", type = "File", prefix = "-U", position = 3)
# p4 <- InputParam(id = "FQ2", type = "File?", prefix = "-2", position = 4)
o1 <- OutputParam(id = "sam", type = "File", prefix = "-S", glob = "*.sam")
req1 <- list(class = "DockerRequirement", 
             dockerPull = "biocontainers/bowtie2:2.3.0-2-deb")
bowtie2 <- cwlParam(baseCommand = c("bowtie2"),
                    requirements = lsit(req1), 
                    inputs = InputParamList(p1, p2, p3),
                    output = OutputParamList(o1),
                    stdout = "bowtie2Output.sam")
