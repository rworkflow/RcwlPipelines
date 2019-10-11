## https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep.html

p1 <- InputParam(id = "svcf", type = "File")
p2 <- InputParam(id = "gvcf", type = "File", secondaryFile = ".tbi")
p3 <- InputParam(id = "ref", type = "File",
                 secondaryFiles = c(".fai", "$(self.nameroot).dict"))
p4 <- InputParam(id = "VepDir", type = "Directory")
p5a <- InputParam(id = "tbam", type = "File", secondaryFile = ".bai")
p5b <- InputParam(id = "rbam", type = "File", secondaryFile = ".bai")
p6 <- InputParam(id = "tsample", type = "string")
p7 <- InputParam(id = "nsample", type = "string")
p8 <- InputParam(id = "rnaseqs", type = "File[]")
p9 <- InputParam(id = "kallistoIdx", type = "File")
p10 <- InputParam(id = "threads", type = "int", default = 16L)
#' @include tl_vep.R
s1 <- Step(id = "VCFvep", run = vep,
           In = list(ivcf = "svcf",
                     ref = "ref",
                     cacheDir = "VepDir",
                     ovcf = list(valueFrom = "$(inputs.ivcf.nameroot)_vep.vcf")))
#' @include pl_vcfCoverage.R
s2a <- Step(id = "dVCFcoverage", run = vcfCoverage,
           In = list(vcf = "VCFvep/oVcf",
                     bam = "tbam",
                     sample = list(valueFrom = "SAMPLE"),
                     ref = "ref"))
s2b <- Step(id = "rVCFcoverage", run = vcfCoverage,
           In = list(vcf = "dVCFcoverage/outvcf",
                     bam = "rbam",
                     sample = list(valueFrom = "SAMPLE"),
                     ntype = list(valueFrom = "RNA"),
                     ref = "ref"))
#' @include pl_vcfExpression.R
s3 <- Step(id = "VCFexpression", run = vcfExpression,
           In = list(rnafqs = "rnaseqs",
                     kallistoIdx = "kallistoIdx",
                     svcf = "rVCFcoverage/outvcf",
                     threads = "threads"))
#' @include pl_phaseVcf.R
s4 <- Step(id = "PhaseVcf", run = phaseVcf,
           In = list(gvariant = "gvcf",
                     svariant = "VCFexpression/ExpVcf",
                     bam = "tbam",
                     outvcf = list(valueFrom = "$(inputs.tsample)_phased.vcf"),
                     nsample = "nsample",
                     tsample = "tsample",
                     ref = "ref"))
o1 <- OutputParam(id = "annVcf", type = "File", outputSource = "VCFexpression/ExpVcf")
o2 <- OutputParam(id = "phasedVCF", type = "File", outputSource = "PhaseVcf/pvcf")
req1 <- list(class = "InlineJavascriptRequirement")
req2 <- list(class = "StepInputExpressionRequirement")
req3 <- list(class = "SubworkflowFeatureRequirement")
ht1 <- list("cwltool:LoadListingRequirement"=
                list(loadListing = "no_listing"))
ext <- list("$namespaces" = list(cwltool = "http://commonwl.org/cwltool#"))
AnnPhaseVcf <- cwlStepParam(requirements = list(req1, req2, req3),
                            inputs = InputParamList(p1, p2, p3, p4, p5a, p5b,
                                                    p6, p7, p8, p9, p10),
                            outputs = OutputParamList(o1, o2),
                            hints = ht1,
                            extensions = ext)
AnnPhaseVcf <- AnnPhaseVcf + s1 + s2a + s2b + s3 + s4
