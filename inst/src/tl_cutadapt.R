## Cutadapt
## Work for paired-end reads
## Input files need to be in one of these formats: 
##   - FASTA with extensions: .fasta, .fa or .fna
##   - FASTQ with extensions: .fastq or .fq
##   - Any of the above, but compressed as .gz, .bz2 or .xz

req1 <- list(class = "DockerRequirement", 
             dockerPull = "kfdrc/cutadapt")

p1 <- InputParam(id = "threadN", type = "int", prefix = "-j", position = 1, default = 1)
p2 <- InputParam(id = "adapter", type = "string", prefix = "-b", position = 2)
p3 <- InputParam(id = "out1prefix", type = "string", prefix = "-o", position = 3)
p4 <- InputParam(id = "out2prefix", type = "string", prefix = "-p", position = 4)
p5 <- InputParam(id = "in1", type = "File", position = 5)
p6 <- InputParam(id = "in2", type = "File", position = 6)
o1 <- OutputParam(id = "out1", type = "File", glob = "$(inputs.out1prefix)")
o2 <- OutputParam(id = "out2", type = "File", glob = "$(inputs.out2prefix)")

cutadapt <- cwlParam(baseCommand = "cutadapt", 
                         requirements = list(req1),
                         inputs = InputParamList(p1, p2, p3, p4, p5, p6), 
                         outputs = OutputParamList(o1, o2))
