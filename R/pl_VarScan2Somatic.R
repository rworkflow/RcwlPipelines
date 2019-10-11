
#' @include tl_samtools_mpileup.R tl_VarScan2_somatic.R tl_VarScan2_processSomatic.R tl_VarScan2_somaticFilter.R

p1 <- InputParam(id = "tbam", type = "File")
p2 <- InputParam(id = "nbam", type = "File")
p3 <- InputParam(id = "ref", type = "File", secondaryFiles = ".fai")
p4 <- InputParam(id = "region", type = "File")

s1 <- Step(id = "mpileupT", run = samtools_mpileup,
           In = list(bam = "tbam",
                     ref = "ref",
                     region = "region"))
s2 <- Step(id = "mpileupN", run = samtools_mpileup,
           In = list(bam = "nbam",
                     ref = "ref",
                     region = "region"))
s3 <- Step(id = "somatic", run = VarScan2_somatic,
           In = list(npileup = "mpileupN/pileup",
                     tpileup = "mpileupT/pileup",
                     bname = list(valueFrom = "$(inputs.tpileup.nameroot)")))
s4 <- Step(id = "processSomatic", run = VarScan2_processSomatic,
           In = list(vcf = "somatic/snp"))
s5 <- Step(id = "somaticFilter", run = VarScan2_somaticFilter,
           In = list(vcf = "processSomatic/somaticHC",
                     indel = "somatic/indel",
                     outvcf = list(source = list("tbam", "nbam"),
                                   valueFrom = "$(self[0].nameroot).$(self[1].nameroot).somatic.vcf")))
o1 <- OutputParam(id = "sSnp", type = "File", outputSource = "somatic/snp")
o2 <- OutputParam(id = "sIndel", type = "File", outputSource = "somatic/indel")
o3 <- OutputParam(id = "sVcf", type = "File", outputSource = "somaticFilter/outVcf")

req1 <- list(class = "StepInputExpressionRequirement")
req2 <- list(class = "MultipleInputFeatureRequirement")
VarScan2Somatic <- cwlStepParam(requirements = list(req1, req2),
                                inputs = InputParamList(p1, p2, p3, p4),
                                outputs = OutputParamList(o1, o2, o3))
VarScan2Somatic <- VarScan2Somatic + s1 + s2 + s3 + s4 + s5
