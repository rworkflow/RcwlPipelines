##source(system.file("tools", "gtfToGenePred.R", package = "RcwlPipelines"))
##source(system.file("tools", "genePredToBed.R", package = "RcwlPipelines"))
##source(system.file("tools", "read_distribution.R", package = "RcwlPipelines"))
##source(system.file("tools", "geneBody_coverage.R", package = "RcwlPipelines"))

## RSeQC Pipeline
#' @include tl_gtfToGenePred.R tl_genePredToBed.R tl_read_distribution.R tl_geneBody_coverage.R
p1 <- InputParam(id = "bam", type = "File", secondaryFiles = ".bai")
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
