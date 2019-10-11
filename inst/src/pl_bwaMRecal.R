## bwaAlign + mergeBamDup
#' @include pl_bwaAlign.R tl_markdup.R pl_BaseRecal.R
p1 <- InputParam(id = "outBam", type = "string")
p2 <- InputParam(id = "RG", type = "string")
p3 <- InputParam(id = "threads", type = "int")
p4 <- InputParam(id = "Ref", type = "File",
                 secondaryFiles = c(".amb", ".ann", ".bwt", ".pac",
                                    ".sa", ".fai",
                                    "$(self.nameroot).dict"))
p5 <- InputParam(id = "FQ1s", type = "File")
p6 <- InputParam(id = "FQ2s", type = "File")
p7 <- InputParam(id = "knowSites", type = InputArrayParam(items = "File"),
                 secondaryFiles = ".idx")

o1 <- OutputParam(id = "BAM", type = "File", outputSource = "BaseRecal/rcBam")
o2 <- OutputParam(id = "matrix", type = "File", outputSource = "markdup/Mat")
## o3 <- OutputParam(id = "Idx", type = "File", outputSource = "mergeBamDup/Idx")
o3 <- OutputParam(id = "flagstat", type = "File", outputSource = "BaseRecal/flagstat")
o4 <- OutputParam(id = "stats", type = "File", outputSource = "BaseRecal/stats")

req1 <- list(class = "SubworkflowFeatureRequirement")
req2 <- list(class = "StepInputExpressionRequirement")
## req3 <- list(class = "InlineJavascriptRequirement")
bwaMRecal <- cwlStepParam(requirements = list(req1, req2),
                          inputs = InputParamList(p1, p2, p3, p4, p5, p6, p7),
                          outputs = OutputParamList(o1, o2, o3, o4))

s1 <- Step(id = "bwaAlign", run = bwaAlign,
           In = list(threads = "threads",
                     RG = "RG",
                     Ref = "Ref",
                     FQ1 = "FQ1s",
                     FQ2 = "FQ2s"))

s2 <- Step(id = "markdup", run = markdup,
           In = list(ibam = "bwaAlign/Idx",
                     obam = list(valueFrom="$(inputs.ibam.nameroot).mdup.bam"),
                     matrix = list(
                         source = list("outBam"),
                         valueFrom="$(self).markdup.txt")))

s3 <- Step(id = "BaseRecal", run = BaseRecal,
           In = list(bam = "markdup/mBam",
                     ref = "Ref",
                     knowSites = "knowSites",
                     oBam = "outBam"))

bwaMRecal <- bwaMRecal + s1 + s2 + s3
