source(system.file("demo", "src", "fastQC.R", package = "RcwlPipelines"))
source(system.file("demo", "src", "featureCounts.R", package = "RcwlPipelines"))
source(system.file("demo", "src", "samtools_flagstat.R", package = "RcwlPipelines"))
source(system.file("demo", "src", "samtools_index.R", package = "RcwlPipelines"))
source(system.file("demo", "src", "RSeQC.R", package = "RcwlPipelines"))
source(system.file("demo", "src", "STAR.R", package = "RcwlPipelines"))

## Pipeline: fastQC + STAR + featureCounts
## Note: output to current dir
p1 <- InputParam(id = "in_seqfiles", type = "File[]")
p2 <- InputParam(id = "in_prefix", type = "string")
p3 <- InputParam(id = "in_genomeDir", type = "Directory")
p4 <- InputParam(id = "in_GTFfile", type = "File")
p5 <- InputParam(id = "in_runThreadN", type = "int", default = 1L)
o1 <- OutputParam(id = "out_fastqc", type = "File[]", outputSource = "fastqc/QCfile")
o2a <- OutputParam(id = "out_BAM", type = "File", outputSource = "STAR/outBAM")
o2b <- OutputParam(id = "out_Log", type = "File", outputSource = "STAR/outLog")
o2c <- OutputParam(id = "out_Count", type = "File", outputSource = "STAR/outCount")
o3 <- OutputParam(id = "out_idx",type = "File", outputSource = "samtools_index/idx")
o4 <- OutputParam(id = "out_stat",type = "File", outputSource = "samtools_flagstat/flagstat")
o5 <- OutputParam(id = "out_count", type = "File", outputSource = "featureCounts/Count")
o6 <- OutputParam(id = "out_distribution", type = "File", outputSource = "RSeQC/distribution")
o7 <- OutputParam(id = "out_gCovP", type = "File", outputSource = "RSeQC/gCovP")
o8 <- OutputParam(id = "out_gCovT", type = "File", outputSource = "RSeQC/gCovT")
req1 <- list(class = "ScatterFeatureRequirement")
req2 <- list(class = "SubworkflowFeatureRequirement")
req3 <- list(class = "StepInputExpressionRequirement")
rnaseq_Sf <- cwlStepParam(requirements = list(req1, req2, req3),
                       inputs = InputParamList(p1, p2, p3, p4, p5),
                       outputs = OutputParamList(o1, o2a, o2b, o2c, o3, o4, o5, o6, o7, o8))

## fastqc
s1 <- Step(id = "fastqc", run = fastqc,
           In = list(seqfile = "in_seqfiles"),
           scatter = "seqfile")
## STAR
s2 <- Step(id = "STAR", run = STAR,
           In = list(prefix = "in_prefix",
                     genomeDir = "in_genomeDir",
                     sjdbGTFfile = "in_GTFfile",
                     readFilesIn = "in_seqfiles",
                     runThreadN = "in_runThreadN"))
## samtools
s3 <- Step(id = "samtools_index", run = samtools_index,
           In = list(bam = "STAR/outBAM"))
s4 <- Step(id = "samtools_flagstat", run = samtools_flagstat,
           In = list(bam = "STAR/outBAM"))
## featureCounts
s5 <- Step(id = "featureCounts", run = featureCounts,
           In = list(gtf = "in_GTFfile",
                     bam = "STAR/outBAM",
                     count = list(valueFrom = "$(inputs.bam.nameroot).featureCounts.txt")))
## RSeQC
s6 <- Step(id = "RSeQC", run = RSeQC,
           In = list(bam = "samtools_index/idx",
                     gtf = "in_GTFfile"))
## pipeline
rnaseq_Sf <- rnaseq_Sf + s1 + s2 + s3 + s4 + s5 + s6
