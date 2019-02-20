## RSeQC Pipeline
## gtf to genePred
p1 <- InputParam(id = "gtf", type = "File", position = 1)
p2 <- InputParam(id = "gPred", type = "string", position = 2)
o1 <- OutputParam(id = "genePred", type = "File", glob = "$(inputs.gPred)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "hubentu/rcwl-rnaseq")
gtfToGenePred <- cwlParam(baseCommand = "gtfToGenePred",
                          requirements = list(req1),
                          inputs = InputParamList(p1, p2),
                          outputs = OutputParamList(o1))

## genePred to bed
p1 <- InputParam(id = "genePred", type = "File", position = 1)
p2 <- InputParam(id = "Bed", type = "string", position = 2)
o1 <- OutputParam(id = "bed", type = "File", glob = "$(inputs.Bed)")
genePredToBed <- cwlParam(baseCommand = "genePredToBed",
                          requirements = list(req1),
                          inputs = InputParamList(p1, p2),
                          outputs = OutputParamList(o1))

## read distribution
p1 <- InputParam(id = "bam", type = "File", prefix = "-i")
p2 <- InputParam(id = "bed", type = "File", prefix = "-r")
o1 <- OutputParam(id = "distOut", type = "File", glob = "$(inputs.bam.nameroot).distribution.txt")
read_distribution <- cwlParam(baseCommand = c("read_distribution.py"),
                              requirements = list(req1),
                              inputs = InputParamList(p1, p2),
                              outputs = OutputParamList(o1),
                              stdout = "$(inputs.bam.nameroot).distribution.txt")

## gene body coverage
p1 <- InputParam(id = "bam", type = "File", prefix = "-i")
p2 <- InputParam(id = "bed", type = "File", prefix = "-r")
p3 <- InputParam(id = "prefix", type = "string", prefix = "-o")
o1 <- OutputParam(id = "gCovPDF", type = "File", glob = "*.geneBodyCoverage.curves.pdf")
o2 <- OutputParam(id = "gCovTXT", type = "File", glob = "*.geneBodyCoverage.txt")
geneBody_coverage <- cwlParam(baseCommand = c("geneBody_coverage.py"),
                              requirements = list(req1),
                              inputs = InputParamList(p1, p2, p3),
                              outputs = OutputParamList(o1, o2))

## Pipeline
p1 <- InputParam(id = "bam", type = "File")
p2 <- InputParam(id = "gtf", type = "File")
o1 <- OutputParam(id = "distribution", type = "File", outputSource = "r_distribution/distOut")
o2 <- OutputParam(id = "gCovP", type = "File", outputSource = "gCoverage/gCovPDF")
o3 <- OutputParam(id = "gCovT", type = "File", outputSource = "gCoverage/gCovTXT")
req1 <- list(class = "StepInputExpressionRequirement")
RSeQC <- cwlStepParam(requirements = list(req1),
                      inputs = InputParamList(p1, p2),
                      outputs = OutputParamList(o1, o2, o3))

s1 <- Step(id = "gtfToGenePred", run = gtfToGenePred,
           In = list(gtf = "gtf",
                     gPred = list(valueFrom = "$(inputs.gtf.nameroot).genePred")))

s2 <- Step(id = "genePredToBed", run = genePredToBed,
           In = list(genePred = "gtfToGenePred/genePred",
                     Bed = list(valueFrom = "$(inputs.genePred.nameroot).bed")))

s3 <- Step(id = "r_distribution", run = read_distribution,
           In = list(bam = "bam",
                     bed = "genePredToBed/bed"))
s4 <- Step(id = "gCoverage", run = geneBody_coverage,
           In = list(bam = "bam",
                     bed = "genePredToBed/bed",
                     prefix = list(valueFrom = "$(inputs.bam.nameroot)")))

RSeQC  <- RSeQC + s1 + s2 + s3 + s4
