## bowtie2_2.3.0-2-deb

p1 <- InputParam(id = "threadN", type = "int", prefix = "-p", position = 1)
# p2 <- InputParam(id = "IndexPrefix", type = "string", prefix = "-x", position = 2) #To be solved 
p3 <- InputParam(id = "fastq", type = "File[]", prefix = "-U", position = 3)
o1 <- OutputParam(id = "sam", type = "File", glob = "*.sam")
req1 <- list(class = "DockerRequirement", 
             dockerPull = "biocontainers/bowtie2:2.3.0-2-deb")
bowtie2 <- cwlParam(baseCommand = "bowtie2",
                    requirements = list(req1),
                    arguments = list("-S", "output.sam"),
                    inputs = InputParamList(p1, p2, p3),
                    outputs = OutputParamList(o1))
