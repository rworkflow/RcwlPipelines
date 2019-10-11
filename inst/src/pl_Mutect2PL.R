## source(system.file("tools", "gatk_CalculateContamination.R", package = "RcwlPipelines"))
## source(system.file("tools", "gatk_ColSeqArtifact.R", package = "RcwlPipelines"))
## source(system.file("tools", "gatk_FilterMutectCalls.R", package = "RcwlPipelines"))
## source(system.file("tools", "gatk_FilterOBias.R", package = "RcwlPipelines"))
## source(system.file("tools", "gatk_Funcotator.R", package = "RcwlPipelines"))
## source(system.file("tools", "gatk_GenomicsDBImport.R", package = "RcwlPipelines"))
## source(system.file("tools", "gatk_GetPileupSummaries.R", package = "RcwlPipelines"))
## source(system.file("tools", "gatk_mutect2.R", package = "RcwlPipelines"))
## source(system.file("tools", "gatk_PoN.R", package = "RcwlPipelines"))
## source(system.file("tools", "bcftools_filter.R", package = "RcwlPipelines"))

##1. Mutect2: Call variant on normal samples, with arguments: --max-mnp-distance 0.

##2. GenomicsDB + PoN, run with cwltool argument: --relax-path-checks

##3. Mutect2 + GetPileupSummaries + CalculateContamination + FilterMutectCalls + CollectSequencingArtifactMetrics + FilterByOrientationBias

p1a <- InputParam(id = "tbam", type = "File", secondaryFiles = ".bai")
p1b <- InputParam(id = "nbam", type = "File", secondaryFiles = ".bai")
p2 <- InputParam(id = "Ref", type = "File",
                 secondaryFiles = c(".fai", "$(self.nameroot).dict"))
p3 <- InputParam(id = "normal", type = "string")
p4 <- InputParam(id = "tumor", type = "string")
p5 <- InputParam(id = "gresource", type = "File", secondaryFiles = ".idx")
p6 <- InputParam(id = "pon", type = "File", secondaryFiles = ".idx")
p7 <- InputParam(id = "interval", type = "File")
p8 <- InputParam(id = "comvcf", type = "File", secondaryFiles = ".idx")
p9 <- InputParam(id = "artMode", type = InputArrayParam(items = "string"),
                 default = list("G/T", "C/T"))
p10 <- InputParam(id = "filter", type = "string", default = "PASS")
#' @include tl_Mutect2.R
s1 <- Step(id = "Mutect2", run = Mutect2,
           In = list(tbam = "tbam",
                     nbam = "nbam",
                     Ref = "Ref",
                     normal = "normal",
                     germline = "gresource",
                     pon = "pon",
                     interval = "interval",
                     out = list(source = list("normal", "tumor"),
                                valueFrom = "$(self[0]).$(self[1])")))
#' @include tl_GetPileupSummaries.R
s2a <- Step(id = "GetPileupSummariesT", run = GetPileupSummaries,
            In = list(bam = "tbam",
                      vcf = "comvcf",
                      interval = "comvcf",
                      pileup = list(valueFrom = "$(inputs.bam.nameroot).ptable")))
s2b <- Step(id = "GetPileupSummariesN", run = GetPileupSummaries,
            In = list(bam = "nbam",
                      vcf = "comvcf",
                      interval = "comvcf",
                      pileup = list(valueFrom = "$(inputs.bam.nameroot).ptable")))
#' @include tl_CalculateContamination.R
s3 <- Step(id = "CalculateContamination", run = CalculateContamination,
           In = list(ttable = "GetPileupSummariesT/pout",
                     ntable = "GetPileupSummariesN/pout",
                     cont = list(source = list("tumor"),
                                 valueFrom = "$(self[0]).contamination.table")))
#' @include tl_FilterMutectCalls.R
s4 <- Step(id = "FilterMutectCalls", run = FilterMutectCalls,
           In = list(vcf = "Mutect2/vout",
                     cont = "CalculateContamination/cout",
                     ref = "Ref",
                     fvcf = list(source = list("normal", "tumor"),
                                 valueFrom = "$(self[0]).$(self[1]).ctfiltered.vcf")))
#' @include tl_ColSeqArtifact.R
s5 <- Step(id = "ColSeqArtifact", run = ColSeqArtifact,
           In = list(bam = "tbam",
                     ref = "Ref",
                     art = list(valueFrom = "$(inputs.bam.nameroot).art")))
#' @include tl_FilterOBias.R
s6 <- Step(id = "FilterOBias", run = FilterOBias,
           In = list(vcf = "FilterMutectCalls/fout",
                     art = "ColSeqArtifact/aout",
                     mode = "artMode",
                     avcf = list(source = list("normal", "tumor"),
                                 valueFrom = "$(self[0]).$(self[1]).ctfiltered.obfiltered.vcf")))
#' @include tl_bcfview.R
s7 <- Step(id = "bcfview", run = bcfview,
           In = list(vcf = "FilterOBias/fout",
                     filter = "filter",
                     fout = list(valueFrom = "$(inputs.vcf.nameroot).PASS.vcf")))

o1 <- OutputParam(id = "filterVCF", type = "File", outputSource = "FilterOBias/fout")
o2 <- OutputParam(id = "passVCF", type = "File", outputSource = "bcfview/Fout")
o3 <- OutputParam(id = "conTable", type = "File", outputSource = "CalculateContamination/cout")
o4 <- OutputParam(id = "artTable", type = "File", outputSource = "ColSeqArtifact/aout")

req1 <- list(class = "InlineJavascriptRequirement")
req2 <- list(class = "StepInputExpressionRequirement")
req3 <- list(class = "MultipleInputFeatureRequirement")
Mutect2PL <- cwlStepParam(requirements = list(req1, req2, req3),
                          inputs = InputParamList(p1a, p1b, p2, p3, p4,
                                                  p5, p6, p7, p8, p9, p10),
                          outputs = OutputParamList(o1, o2, o3, o4))

Mutect2PL <- Mutect2PL + s1 + s2a + s2b + s3 + s4 +s5 + s6 + s7
