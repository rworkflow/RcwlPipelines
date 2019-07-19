##source(system.file("pipelines", "DNASeq_bwa.R", package = "RcwlPipelines"))
##source(system.file("pipelines", "DNASeq_merge.R", package = "RcwlPipelines"))

## bwaAlign + mergeBamDup
#' @include pl_bwaAlign.R pl_mergeBamDup.R
p1 <- InputParam(id = "idBam", type = "string")
p2 <- InputParam(id = "RG", type = "string[]")
p3 <- InputParam(id = "threads", type = "int")
p4 <- InputParam(id = "Ref", type = "File",
                 secondaryFiles = c(".amb", ".ann", ".bwt", ".pac", ".sa"))
p5 <- InputParam(id = "FQ1s", type = "File[]")
p6 <- InputParam(id = "FQ2s", type = "File[]")

o1 <- OutputParam(id = "oBam", type = "File", outputSource = "mergeBamDup/oBam")
o2 <- OutputParam(id = "matrix", type = "File", outputSource = "mergeBamDup/matrix")
o3 <- OutputParam(id = "Idx", type = "File", outputSource = "mergeBamDup/Idx")
o4 <- OutputParam(id = "stat", type = "File", outputSource = "mergeBamDup/stat")

req1 <- list(class = "SubworkflowFeatureRequirement")
req2 <- list(class = "ScatterFeatureRequirement")
alignMerge <- cwlStepParam(requirements = list(req1, req2),
                           inputs = InputParamList(p1, p2, p3, p4, p5, p6),
                           outputs = OutputParamList(o1, o2, o3, o4)
                           )

s1 <- Step(id = "bwaAlign", run = bwaAlign,
           In = list(threads = "threads",
                     RG = "RG",
                     Ref = "Ref",
                     FQ1 = "FQ1s",
                     FQ2 = "FQ2s"),
           scatter = list("RG", "FQ1", "FQ2"),
           scatterMethod = "dotproduct")

s2 <- Step(id = "mergeBamDup", run = mergeBamDup,
           In = list(ibam = "bwaAlign/Bam",
                     obam = "idBam"))

alignMerge  <- alignMerge + s1 + s2
