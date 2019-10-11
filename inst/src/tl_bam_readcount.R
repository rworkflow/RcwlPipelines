## https://github.com/genome/bam-readcount

p1 <- InputParam(id = "vcf", type = "File", position = 1)
p2 <- InputParam(id = "sample", type = "string", position =2)
p3 <- InputParam(id = "ref", type = "File",
                 secondaryFiles = ".fai", position = 3)
p4 <- InputParam(id = "bam", type = "File",
                 secondaryFiles = ".bai", position = 4)
o1 <- OutputParam(id = "snv", type = "File",
                  glob = "$(inputs.sample)_bam_readcount_snv.tsv")
o2 <- OutputParam(id = "indel", type = "File",
                  glob = "$(inputs.sample)_bam_readcount_indel.tsv")

req1 <- list(class = "DockerRequirement",
             dockerPull = "mgibio/bam_readcount_helper-cwl:1.1.1")

bam_readcount <- cwlParam(baseCommand = c("/usr/bin/python",
                                          "/usr/bin/bam_readcount_helper.py"),
                          requirements = list(req1),
                          arguments = list(
                              list(valueFrom = "NOPREFIX", position = 5L),
                              list(valueFrom = "./", position = 6L, shellQuote = FALSE)),
                          inputs = InputParamList(p1, p2, p3, p4),
                          outputs = OutputParamList(o1, o2))
