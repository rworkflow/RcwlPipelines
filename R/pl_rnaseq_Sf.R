##source(system.file("tools", "fastQC.R", package = "RcwlPipelines"))
##source(system.file("tools", "featureCounts.R", package = "RcwlPipelines"))
##source(system.file("tools", "samtools_flagstat.R", package = "RcwlPipelines"))
##source(system.file("tools", "samtools_index.R", package = "RcwlPipelines"))
##source(system.file("pipelines", "RSeQC.R", package = "RcwlPipelines"))
##source(system.file("tools", "STAR.R", package = "RcwlPipelines"))

## Pipeline: fastQC + STAR + featureCounts
## Note: output to current dir
p1 <- InputParam(id = "in_seqfiles", type = "File[]")
p2 <- InputParam(id = "in_prefix", type = "string")
p3 <- InputParam(id = "in_genomeDir", type = "Directory")
p4 <- InputParam(id = "in_GTFfile", type = "File")
p5 <- InputParam(id = "in_runThreadN", type = "int", default = 1L)
o1 <- OutputParam(id = "out_fastqc", type = "File[]", outputSource = "fastqc/QCfile")
o2a <- OutputParam(id = "out_BAM", type = "File",
                   outputSource = "samtools_index/idx")
o2b <- OutputParam(id = "out_Log", type = "File", outputSource = "STAR/outLog")
o2c <- OutputParam(id = "out_Count", type = "File", outputSource = "STAR/outCount")
##o3 <- OutputParam(id = "out_idx",type = "File", outputSource = "samtools_index/idx")
o4 <- OutputParam(id = "out_stat",type = "File", outputSource = "samtools_flagstat/flagstat")
o5 <- OutputParam(id = "out_count", type = "File", outputSource = "featureCounts/Count")
o6 <- OutputParam(id = "out_distribution", type = "File", outputSource = "r_distribution/distOut")
o7 <- OutputParam(id = "out_gCovP", type = "File", outputSource = "gCoverage/gCovPDF")
o8 <- OutputParam(id = "out_gCovT", type = "File", outputSource = "gCoverage/gCovTXT")
req1 <- list(class = "ScatterFeatureRequirement")
req2 <- list(class = "SubworkflowFeatureRequirement")
req3 <- list(class = "StepInputExpressionRequirement")
rnaseq_Sf <- cwlStepParam(requirements = list(req1, req2, req3),
                       inputs = InputParamList(p1, p2, p3, p4, p5),
                       outputs = OutputParamList(o1, o2a, o2b, o2c, o4, o5, o6, o7, o8))

## fastqc
#' @include tl_fastqc.R
s1 <- Step(id = "fastqc", run = fastqc,
           In = list(seqfile = "in_seqfiles"),
           scatter = "seqfile")
## STAR
#' @include tl_STAR.R
## fix STAR sort problem
arguments(STAR)[[13]] <- "Unsorted"
s2 <- Step(id = "STAR", run = STAR,
           In = list(prefix = "in_prefix",
                     genomeDir = "in_genomeDir",
                     sjdbGTFfile = "in_GTFfile",
                     readFilesIn = "in_seqfiles",
                     runThreadN = "in_runThreadN"))
## samtools
#' @include tl_sortBam.R tl_samtools_index.R tl_samtools_flagstat.R
s3a <- Step(id = "sortBam", run = sortBam,
            In = list(bam = "STAR/outBAM"))
s3 <- Step(id = "samtools_index", run = samtools_index,
           In = list(bam = "sortBam/sbam"))
s4 <- Step(id = "samtools_flagstat", run = samtools_flagstat,
           In = list(bam = "sortBam/sbam"))
## featureCounts
#' @include tl_featureCounts.R
s5 <- Step(id = "featureCounts", run = featureCounts,
           In = list(gtf = "in_GTFfile",
                     bam = "samtools_index/idx",
                     count = list(valueFrom = "$(inputs.bam.nameroot).featureCounts.txt")))
## RSeQC
## #' @include pl_RSeQC.R
## s6 <- Step(id = "RSeQC", run = RSeQC,
##            In = list(bam = "samtools_index/idx",
##                      gtf = "in_GTFfile"))
#' @include tl_gtfToGenePred.R tl_genePredToBed.R tl_read_distribution.R tl_geneBody_coverage.R
s6a <- Step(id = "gtfToGenePred", run = gtfToGenePred,
            In = list(gtf = "in_GTFfile",
                      gPred = list(valueFrom = "$(inputs.gtf.nameroot).genePred")))

s6b <- Step(id = "genePredToBed", run = genePredToBed,
            In = list(genePred = "gtfToGenePred/genePred",
                      Bed = list(valueFrom = "$(inputs.genePred.nameroot).bed")))

s6c <- Step(id = "r_distribution", run = read_distribution,
            In = list(bam = "samtools_index/idx",
                      bed = "genePredToBed/bed"))
gCoverage <- geneBody_coverage
gCoverage@inputs$bam@secondaryFiles <- character()
s6d <- Step(id = "gCoverage", run = gCoverage,
            In = list(bam = "samtools_index/idx",
                      bed = "genePredToBed/bed",
                      prefix = list(valueFrom = "$(inputs.bam.nameroot)")))


## pipeline
rnaseq_Sf <- rnaseq_Sf + s1 + s2 + s3a + s3 + s4 + s5 + s6a + s6b + s6c + s6d
