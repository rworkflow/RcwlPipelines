## STAR
## Note: output to current dir
## readFilesIn must be full path ##string
p1 <- InputParam(id = "prefix", type = "string", prefix = "--outFileNamePrefix")
p2 <- InputParam(id = "readFilesIn", type = "File[]", prefix = "--readFilesIn")
p3 <- InputParam(id = "genomeDir", type = "Directory", prefix = "--genomeDir")
p4 <- InputParam(id = "sjdbGTFfile", type = "File", prefix = "--sjdbGTFfile")
p5 <- InputParam(id = "runThreadN", type = "int", prefix = "--runThreadN", default = 1L)
o1 <- OutputParam(id = "outBAM", type = "File", glob = "*.bam")
o2 <- OutputParam(id = "outLog", type = "File", glob = "*Log.final.out")
o3 <- OutputParam(id = "outCount", type = "File", glob = "*ReadsPerGene.out.tab")

req1 <- list(class = "DockerRequirement",
             dockerPull = "hubentu/rcwl-rnaseq")
STAR <- cwlParam(baseCommand = "STAR",
                 requirements = list(req1),
                 arguments = list("--outFilterMultimapNmax", "3",
                                  "--outSAMunmapped", "Within",
                                  "--outFilterMismatchNmax", "2",
                                  "--outSAMstrandField", "intronMotif",
                                  "--readFilesCommand", "zcat",
                                  "--outSAMtype", "BAM", "SortedByCoordinate",
                                  "--twopassMode", "Basic",
                                  "--quantMode", "GeneCounts"),
                 inputs = InputParamList(p1, p2, p3, p4, p5),
                 outputs = OutputParamList(o1, o2, o3))
