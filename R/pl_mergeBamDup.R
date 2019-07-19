##source(system.file("tools", "mergeBam.R", package = "RcwlPipelines"))
##source(system.file("tools", "markdup.R", package = "RcwlPipelines"))
##source(system.file("tools", "samtools_index.R", package = "RcwlPipelines"))
##source(system.file("tools", "samtools_flagstat.R", package = "RcwlPipelines"))


#' @include tl_mergeBam.R tl_markdup.R tl_samtools_index.R tl_samtools_flagstat.R
## mergeBam + markdup + index + flagstat
p1 <- InputParam(id = "ibam", type = "File[]")
p2 <- InputParam(id = "obam", type = "string")
o1 <- OutputParam(id = "oBam", type = "File", outputSource = "markdup/mBam")
o2 <- OutputParam(id = "matrix", type = "File", outputSource = "markdup/Mat")
o3 <- OutputParam(id = "Idx", type = "File", outputSource = "samtools_index/idx")
o4 <- OutputParam(id = "stat", type = "File",
                  outputSource = "samtools_flagstat/flagstat")
req1 <- list(class = "StepInputExpressionRequirement")
req2 <- list(class = "InlineJavascriptRequirement")
mergeBamDup <- cwlStepParam(requirements = list(req1, req2),
                            inputs = InputParamList(p1, p2),
                            outputs = OutputParamList(o1, o2, o3, o4))
s1 <- Step(id = "mergeBam", run = mergeBam,
           In = list(ibam = "ibam",
                     obam = "obam"))
s2 <- Step(id = "markdup", run = markdup,
           In = list(ibam = "mergeBam/oBam",
                     obam = "obam",
                     matrix = list(valueFrom="$(inputs.ibam.nameroot).markdup.txt")))
s3 <- Step(id = "samtools_index", run = samtools_index,
           In = list(bam = "markdup/mBam"))
s4 <- Step(id = "samtools_flagstat", run = samtools_flagstat,
           In = list(bam = "markdup/mBam"))
mergeBamDup <- mergeBamDup + s1 + s2 + s3 + s4
