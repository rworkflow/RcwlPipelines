## BAM + BaseRecalibrator + ApplyBQSR + index + stats
#' @include  tl_BaseRecalibrator.R tl_ApplyBQSR.R tl_samtools_index.R tl_samtools_flagstat.R tl_samtools_stats.R
p1 <- InputParam(id = "bam", type = "File")
p2 <- InputParam(id = "ref", type = "File", secondaryFiles = c(".fai",
                                                               "$(self.nameroot).dict"))
p3 <- InputParam(id = "knowSites", type = InputArrayParam(items = "File"),
                 secondaryFiles = ".idx")
p4 <- InputParam(id = "oBam", type = "string")

s1 <- Step(id = "BaseRecalibrator", run = BaseRecalibrator,
           In = list(bam = "bam",
                     ref = "ref",
                     knowSites = "knowSites",
                     recal = list(
                         valueFrom="$(inputs.bam.nameroot).recal.txt")))
s2 <- Step(id = "ApplyBQSR", run = ApplyBQSR,
           In = list(bam = "bam",
                     ref = "ref",
                     rtable = "BaseRecalibrator/rtable",
                     oBam = "oBam"))
s3 <- Step(id = "samtools_index", run = samtools_index,
           In = list(bam = "ApplyBQSR/Bam"))
s4 <- Step(id = "samtools_flagstat", run = samtools_flagstat,
           In = list(bam = "ApplyBQSR/Bam"))
s5 <- Step(id = "samtools_stats", run = samtools_stats,
           In = list(bam = "ApplyBQSR/Bam"))

o1 <- OutputParam(id = "rcBam", type = "File", outputSource = "samtools_index/idx")
o2 <- OutputParam(id = "flagstat", type = "File",
                  outputSource = "samtools_flagstat/flagstat")
o3 <- OutputParam(id = "stats", type = "File",
                  outputSource = "samtools_stats/stats")

req1 <- list(class = "StepInputExpressionRequirement")
req2 <- list(class = "InlineJavascriptRequirement")
BaseRecal <- cwlStepParam(requirements = list(req1, req2),
                            inputs = InputParamList(p1, p2, p3, p4),
                            outputs = OutputParamList(o1, o2, o3))
BaseRecal <- BaseRecal + s1 + s2 + s3 + s4 + s5
