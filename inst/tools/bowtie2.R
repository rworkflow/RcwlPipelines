## bowtie2:v2.2.9_cv2
## Be careful: bowtie2-build genome.fa genome.fa 
## Basename of .bt2 files must be exactly equal to "genome.fa" which is the reference_in

p1 <- InputParam(id = "threadN", type = "int", prefix = "-p", position = 1)
p2 <- InputParam(id = "IndexPrefix", type = "File", prefix = "-x",
                 valueFrom = "$(self.dirname + '/' + self.nameroot)",
                 secondaryFiles = c("$(self.nameroot + '.1.bt2')",
                                    "$(self.nameroot + '.2.bt2')", 
                                    "$(self.nameroot + '.3.bt2')",
                                    "$(self.nameroot + '.4.bt2')",
                                    "$(self.nameroot + '.rev.1.bt2')", 
                                    "$(sefl.nameroot + '.rev.2.bt2')"), 
                 position = 2) 
p3 <- InputParam(id = "fastq", type = "File[]", prefix = "-U", position = 3)
o1 <- OutputParam(id = "sam", type = "File", glob = "*.sam")
req1 <- list(class = "InlineJavascriptRequirement")
req2 <- list(class = "DockerRequirement", 
             dockerPull = "biocontainers/bowtie2:v2.2.9_cv2")
bowtie2 <- cwlParam(baseCommand = "bowtie2",
                    requirements = list(req1, req2),
                    arguments = list("-S", "output.sam"),
                    inputs = InputParamList(p1, p2, p3),
                    outputs = OutputParamList(o1))
