## Salmon index
## Reference (refMrna.fa) can be downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/

rep1 <- list(class = "DockerRequirement",
             dockerPull = "combinelab/salmon")
rep2 <- list(class = "InlineJavascriptRequirement")

p1 <- InputParam(id = "threadN", type = "int", prefix = "-p", position = 1)
p2 <- InputParam(id = "kmer", type = "int", prefix = "-k", position = 2)
p3 <- InputParam(id = "refFasta", type = "File", prefix = "-t", position = 3)
p4 <- InputParam(id = "outPrefix", type = "string", prefix = "-i", position = 4)
o1 <- OutputParam(id = "out1", type = "Directory", glob = "$(inputs.outPrefix)")

salmon_index <- cwlParam(baseCommand = c("salmon", "index"), 
                             requirements = list(rep1, rep2),
                             arguments = list("--type", "quasi"),
                             inputs = InputParamList(p1, p2, p3, p4), 
                             outputs = OutputParamList(o1))
