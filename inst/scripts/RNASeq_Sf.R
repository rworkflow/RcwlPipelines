## Pipeline: fastQC + STAR + featureCounts
## Note: output to current dir
p1 <- InputParam(id = "cwl_seqfiles", type = "File[]")
p2 <- InputParam(id = "cwl_prefix", type = "string")
p3 <- InputParam(id = "cwl_genomeDir", type = "Directory")
p4 <- InputParam(id = "cwl_sjdbGTFfile", type = "File")
p5 <- InputParam(id = "cwl_runThreadN", type = "int", default = 1L)
o1 <- OutputParam(id = "cwl_fastqc", type = "File[]", outputSource = "fastqc/QCfile")
o2a <- OutputParam(id = "cwl_outBAM", type = "File", outputSource = "STAR/outBAM")
o2b <- OutputParam(id = "cwl_outLog", type = "File", outputSource = "STAR/outLog")
o2c <- OutputParam(id = "cwl_outCount", type = "File", outputSource = "STAR/outCount")
o3 <- OutputParam(id = "cwl_idx",type = "File", outputSource = "s_index/idx")
o4 <- OutputParam(id = "cwl_stat",type = "File", outputSource = "s_flagstat/flagstat")
o5 <- OutputParam(id = "cwl_count", type = "File", outputSource = "fcount/count")
o6 <- OutputParam(id = "cwl_dist", type = "File", outputSource = "RSeQC/distribution")
    
req1 <- list(class = "ScatterFeatureRequirement")
req2 <- list(class = "SubworkflowFeatureRequirement")
req3 <- list(class = "StepInputExpressionRequirement")
rnaseq_Sf <- cwlStepParam(requirements = list(req1, req2, req3),
                       inputs = InputParamList(p1, p2, p3, p4, p5),
                       outputs = OutputParamList(o1, o2a, o2b, o2c, o3, o4, o5, o6))

## fastqc
s1 <- Step(id = "fastqc", run = fastqc,
           In = list(seqfile = "cwl_seqfiles"),
           scatter = "seqfile")
## STAR
s2 <- Step(id = "STAR", run = STAR,
           In = list(prefix = "cwl_prefix",
                     genomeDir = "cwl_genomeDir",
                     sjdbGTFfile = "cwl_sjdbGTFfile",
                     readFilesIn = "cwl_seqfiles",
                     runThreadN = "cwl_runThreadN"))
## samtools
s3 <- Step(id = "s_index", run = samtools_index,
           In = list(bam = "STAR/outBAM"))
s4 <- Step(id = "s_flagstat", run = samtools_flagstat,
           In = list(bam = "STAR/outBAM"))
## featureCounts
s5 <- Step(id = "fcount", run = featureCounts,
           In = list(gtf = "cwl_sjdbGTFfile",
                     bam = "STAR/outBAM",
                     count = list(valueFrom = "$(inputs.bam.nameroot).featureCounts.txt")))
## RSeQC
s6 <- Step(id = "RSeQC", run = RSeQC,
           In = list(bam = "STAR/outBAM",
                     gtf = "cwl_sjdbGTFfile"))
## pipeline
rnaseq_Sf <- rnaseq_Sf + s1 + s2 + s3 + s4 + s5 + s6
