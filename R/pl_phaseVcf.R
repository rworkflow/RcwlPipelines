## https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep/proximal_vcf.html

p1a <- InputParam(id = "gvariant", type = "File", secondaryFile = ".tbi")
p1b <- InputParam(id = "svariant", type = "File", secondaryFile = ".tbi")
p2 <- InputParam(id = "ref", type = "File",
                 secondaryFiles = c(".fai", "$(self.nameroot).dict"))
p3 <- InputParam(id = "bam", type = "File", secondaryFiles = ".bai")
p4 <- InputParam(id = "outvcf", type = "string")
p5 <- InputParam(id = "nsample", type = "string")
p6 <- InputParam(id = "tsample", type = "string")
#' @include tl_bcfview.R
s1 <- Step(id = "splitSample", run = bcfview,
           In = list(vcf = "gvariant",
                     sample = "nsample",
                     fout = list(valueFrom = "$(inputs.sample)_germline.vcf"),
                     genotype = list(valueFrom = "^miss"),
                     exclude = list(valueFrom = "GT='0/0'")))
#' @include tl_RenameSampleInVcf.R
s2a <- Step(id = "renameGVcf", run = RenameSampleInVcf,
           In = list(vcf = "splitSample/Fout",
                     ovcf = list(valueFrom = "$(inputs.vcf.nameroot)_g.vcf"),
                     NewName = "tsample"))
s2b <- Step(id = "renameSVcf", run = RenameSampleInVcf,
           In = list(vcf = "svariant",
                     ovcf = list(valueFrom = "$(inputs.vcf.nameroot)_s.vcf"),
                     NewName = "tsample"))
#' @include tl_CombineVariants.R
s3 <- Step(id = "combineVariants", run = CombineVariants,
           In = list(variants = list(source = list("renameGVcf/oVcf", "renameSVcf/oVcf")),
                     ref = "ref",
                     ovcf = list(valueFrom = "combined_somatic_germline.vcf")))
#' @include tl_SortVcf.R
s4 <- Step(id = "sortVcf", run = SortVcf,
           In = list(vcf = "combineVariants/oVcf",
                     ovcf = list(valueFrom = "$(inputs.vcf.nameroot)_sorted.vcf")))
#' @include tl_ReadBackedPhasing.R
s5 <- Step(id = "ReadBackedPhasing", run = ReadBackedPhasing,
           In = list(vcf = "sortVcf/oVcf",
                     bam = "bam",
                     ref = "ref",
                     region = "sortVcf/oVcf",
                     ovcf = "outvcf"))
#' @include tl_bgzip.R
s6 <- Step(id = "bgzip", run = bgzip,
           In = list(ifile = "ReadBackedPhasing/oVcf"))
#' @include tl_tabix_index.R
s7 <- Step(id = "tabixIndex", run = tabix_index,
           In = list(tfile = "bgzip/zfile",
                     type = list(valueFrom = "vcf")))
o1 <- OutputParam(id = "pvcf", type = "File", outputSource = "tabixIndex/idx")
req1 <- list(class = "InlineJavascriptRequirement")
req2 <- list(class = "StepInputExpressionRequirement")
req3 <- list(class = "MultipleInputFeatureRequirement")
phaseVcf <- cwlStepParam(requirements = list(req1, req2, req3),
                         inputs = InputParamList(p1a, p1b, p2, p3, p4, p5, p6),
                         outputs = OutputParamList(o1))
phaseVcf <- phaseVcf + s1 + s2a + s2b + s3 + s4 + s5 + s6 + s7
