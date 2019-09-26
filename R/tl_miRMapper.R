## https://github.com/rajewsky-lab/mirdeep2
## mapper.pl reads.fa -c -j -k TCGTATGCCGTCTTCTGCTTGT  -l 18 -m -p cel_cluster -s reads_collapsed.fa -t reads_collapsed_vs_genome.arf -v -n

p1 <- InputParam(id = "reads", type = "File", position = 1)
p2 <- InputParam(id = "format", type = "string",
                 position = 2, default = "-c")
p3 <- InputParam(id = "adapter", type = "string",
                 prefix = "-k", position = 3)
p4 <- InputParam(id = "len", type = "int",
                 prefix = "-l", default = 18L, position = 4)
p5 <- InputParam(id = "genome", type = "File", prefix = "-p",
                 valueFrom = "$(self.dirname + '/' + self.nameroot)",
                 secondaryFiles = c("$(self.nameroot + '.1.ebwt')",
                                    "$(self.nameroot + '.2.ebwt')",
                                    "$(self.nameroot + '.3.ebwt')",
                                    "$(self.nameroot + '.4.ebwt')",
                                    "$(self.nameroot + '.rev.1.ebwt')",
                                    "$(self.nameroot + '.rev.2.ebwt')"),
                 position = 5)
p6 <- InputParam(id = "preads", type = "string", prefix = "-s", position = 6)
p7 <- InputParam(id = "arf", type = "string", prefix = "-t", position = 7)
o1 <- OutputParam(id = "pReads", type = "File", glob = "$(inputs.preads)")
o2 <- OutputParam(id = "Arf", type = "File", glob = "$(inputs.arf)")

req1 <- list(class = "DockerRequirement",
             dockerPull = "hubentu/mirdeep2")
req2 <- list(class = "InlineJavascriptRequirement")
miRMapper <- cwlParam(baseCommand = "mapper.pl",
                      requirements = list(req1, req2),
                      arguments = list(
                          list(valueFrom = "-j", position = 8L),
                          list(valueFrom = "-m", position = 9L)),
                      inputs = InputParamList(p1, p2, p3, p4, p5, p6, p7),
                      outputs = OutputParamList(o1, o2))
