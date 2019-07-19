## hisat2:v2.0.5-1-deb_cv1

req1 <- list(class = "DockerRequirement",
             dockerPull = "biocontainers/hisat2:v2.0.5-1-deb_cv1")
req2 <- list(class = "InlineJavascriptRequirement")

p1 <- InputParam(id = "threadN", type = "int", prefix = "-p")
p2 <- InputParam(id = "IndexPrefix", type = "File", prefix = "-x", 
                  valueFrom = "$(self.dirname + '/' + self.basename)",
                  secondaryFiles = c("$(self.basename + '.1.ht2')", 
                                     "$(self.basename + '.2.ht2')",
                                     "$(self.basename + '.3.ht2')",
                                     "$(self.basename + '.4.ht2')",
                                     "$(self.basename + '.5.ht2')",
                                     "$(self.basename + '.6.ht2')",
                                     "$(self.basename + '.7.ht2')",
                                     "$(self.basename + '.8.ht2')"))
p3 <- InputParam(id = "fq1", type = "File", prefix = "-1")
p4 <- InputParam(id = "fq2", type = "File", prefix = "-2")
o1 <- OutputParam(id = "sam", type = "File", glob = "*.sam")

hisat2_align <- cwlParam(baseCommand = "hisat2",
                         requirements = list(req1, req2),
                         inputs = InputParamList(p1, p2, p3, p4),
                         outputs = OutputParamList(o1),
                         stdout = "output.sam")
