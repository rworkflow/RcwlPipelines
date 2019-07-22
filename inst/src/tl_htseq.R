## HTseq
## Input BAM needs to be sorted by position or name 
## For gene quantification

req1 <- list(class = "DockerRequirement",
             dockerPull = "genomicpariscentre/htseq")

p1 <- InputParam(id = "minaqual", type = "int", prefix = "-a")
p2 <- InputParam(id = "stranded", type = "string", prefix = "-s")
p3 <- InputParam(id = "bam", type = "File")
p4 <- InputParam(id = "gtf", type = "File")
o1 <- OutputParam(id = "out", type = "File", glob = "$(inputs.bam.nameroot).htseq.txt")

htseq <- cwlParam(baseCommand = "htseq-count",
                 requirements = list(req1),
                 arguments = list("--format", "bam",
                                  "--mode", "intersection-strict"),
                 inputs = InputParamList(p1, p2, p3, p4),
                 outputs = OutputParamList(o1),
                 stdout = "$(inputs.bam.nameroot).htseq.txt")
