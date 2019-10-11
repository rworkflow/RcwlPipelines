
p1 <- InputParam(id = "vcf", type = "File")
p2 <- InputParam(id = "sample", type = "string")
p3 <- InputParam(id = "bam", type = "File", secondaryFiles = ".bai")
p4 <- InputParam(id = "ntype", type = "string", default = "DNA")
p5 <- InputParam(id = "ref", type = "File", secondaryFiles = ".fai")
#' @include tl_vt_decompose.R
s1 <- Step(id = "decompose", run = vt_decompose,
           In = list(ivcf = "vcf",
                     ovcf = list(valueFrom = "$(inputs.ivcf.nameroot)_dc.vcf")))
#' @include tl_bam_readcount.R
s2 <- Step(id = "readcount", run = bam_readcount,
           In = list(vcf = "decompose/oVcf",
                     sample = "sample",
                     ref = "ref",
                     bam = "bam"))
#' @include tl_vcf_readcount_annotator.R
s3 <- Step(id = "readcount_annotator_snv", run = vcf_readcount_annotator,
           In = list(ivcf = "decompose/oVcf",
                     readcount = "readcount/snv",
                     ntype = "ntype",
                     sample = "sample",
                     ovcf = list(valueFrom = "$(inputs.ivcf.nameroot)_snv.vcf")))
s4 <- Step(id = "readcount_annotator_indel", run = vcf_readcount_annotator,
           In = list(ivcf = "readcount_annotator_snv/oVcf",
                     readcount = "readcount/indel",
                     ntype = "ntype",
                     sample = "sample",
                     vtype = list(valueFrom = "indel"),
                     ovcf = list(valueFrom = "$(inputs.ivcf.nameroot)_indel.vcf")))

o1 <- OutputParam(id = "outvcf", type = "File",
                  outputSource = "readcount_annotator_indel/oVcf")
req1 <- list(class = "InlineJavascriptRequirement")
req2 <- list(class = "StepInputExpressionRequirement")
req3 <- list(class = "SubworkflowFeatureRequirement")
vcfCoverage <- cwlStepParam(requirements = list(req1, req2, req3),
                            inputs = InputParamList(p1, p2, p3, p4, p5),
                            outputs = OutputParamList(o1))
vcfCoverage <- vcfCoverage + s1 + s2 + s3 + s4
