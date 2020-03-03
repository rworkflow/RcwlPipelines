
p1 <- InputParam(id = "tbam", type = "File", secondaryFiles = ".bai")
p2 <- InputParam(id = "nbam", type = "File", secondaryFiles = ".bai")
p3 <- InputParam(id = "Ref", type = "File",
                 secondaryFiles = c(".fai", "$(self.nameroot).dict"))
## p4 <- InputParam(id = "normal", type = "string")
## p5 <- InputParam(id = "tumor", type = "string")
p6 <- InputParam(id = "dbsnp", type = "File", secondaryFiles = ".tbi")
p7 <- InputParam(id = "gresource", type = "File", secondaryFiles = ".idx")
p8 <- InputParam(id = "pon", type = "File", secondaryFiles = ".idx")
p9 <- InputParam(id = "interval", type = "File")
p10 <- InputParam(id = "comvcf", type = "File", secondaryFiles = ".idx")
p11 <- InputParam(id = "artMode", type = InputArrayParam(items = "string"),
                 default = list("G/T", "C/T"))
p12 <- InputParam(id = "filter", type = "string", default = "PASS")
p13 <- InputParam(id = "threads", type = "int", default = 8)

#' @include pl_Mutect2PL.R
s1 <- Step(id = "Mutect2PL", run = Mutect2PL,
           In = list(tbam = "tbam",
                     nbam = "nbam",
                     Ref = "Ref",
                     normal = list(valueFrom = "$(inputs.nbam.nameroot.split('_')[0])"),
                     tumor = list(valueFrom = "$(inputs.tbam.nameroot.split('_')[0])"),
                     gresource = "gresource",
                     pon = "pon",
                     interval = "interval",
                     comvcf = "comvcf"))
#' @include tl_MuSE.R
s2 <- Step(id = "MuSE", run = MuSE,
           In = list(tbam = "tbam",
                     nbam = "nbam",
                     ref = "Ref",
                     region = "interval",
                     dbsnp = "dbsnp",
                     vcf = list(valueFrom = "$(inputs.tbam.nameroot.split('_')[0])_MuSE.vcf")))
#' @include pl_mantaStrelka.R
s3a <- Step(id = "bgzip", run = bgzip,
            In = list(ifile = "interval"))
s3b <- Step(id = "tabixIndex", run = tabix_index,
            In = list(tfile = "bgzip/zfile",
                      type = list(valueFrom = "bed")))
s3 <- Step(id = "mantaStrelka", run = mantaStrelka,
           In = list(tbam = "tbam",
                     nbam = "nbam",
                     ref = "Ref",
                     region = "tabixIndex/idx"))
#' @include tl_SomaticSniper.R
s4 <- Step(id = "SomaticSniper", run = SomaticSniper,
           In = list(tbam = "tbam",
                     nbam = "nbam",
                     ref = "Ref",
                     vcf = list(valueFrom = "$(inputs.tbam.nameroot.split('_')[0])_SomaticSniper.vcf")))
#' @include tl_VarDict.R
s5 <- Step(id = "VarDict", run = VarDict,
           In = list(tbam = "tbam",
                     nbam = "nbam",
                     ref = "Ref",
                     region = "interval",
                     vcf = list(valueFrom = "$(inputs.tbam.nameroot.split('_')[0])_VarDict.vcf")))
#' @include tl_LoFreq.R
s6 <- Step(id = "LoFreq", run = LoFreq,
           In = list(tbam = "tbam",
                     nbam = "nbam",
                     ref = "Ref",
                     region = "interval",
                     dbsnp = "dbsnp",
                     threads = "threads",
                     out = list(valueFrom = "$(inputs.tbam.nameroot.split('_')[0])_LoFreq.vcf")))
#' @include pl_VarScan2Somatic.R
s7 <- Step(id = "VarScanPL", run = VarScan2Somatic,
           In = list(tbam = "tbam",
                     nbam = "nbam",
                     ref = "Ref",
                     region = "interval"))

#' @include tl_SomaticSeq_Wrapper.R
s8 <- Step(id = "Wrapper", run = SomaticSeq_Wrapper,
           In = list(tbam = "tbam",
                     nbam = "nbam",
                     ref = "Ref",
                     region = "interval",
                     dbsnp = "dbsnp",
                     mutect2 = "Mutect2PL/filterVCF",
                     varscanSnv = "VarScanPL/sSnp",
                     varscanIndel = "VarScanPL/sIndel",
                     sniper = "SomaticSniper/outVcf",
                     vardict = "VarDict/outVcf",
                     muse = "MuSE/outVcf",
                     strelkaSnv = "mantaStrelka/snvs",
                     strelkaIndel = "mantaStrelka/indels",
                     lofreqSnv = "LoFreq/snp",
                     lofreqIndel = "LoFreq/indel"))

mergeTSV <- function(esnv, eindel){
    snv1 <- read.table(esnv, header = TRUE)
    indel1 <- read.table(eindel, header = TRUE)
    var1 <- rbind(snv1, indel1)
    var1[is.na(var1)] <- 0
    write.table(var1, file = "Ensemble.sVar.tsv",
                row.names = FALSE, sep = "\t", quote = FALSE)
}

m1 <- InputParam(id = "esnv", type = "File", prefix = "esnv=", separate = FALSE)
m2 <- InputParam(id = "eindel", type = "File", prefix = "eindel=", separate = FALSE)
out1 <- OutputParam(id = "tsv", type = "File", glob = "Ensemble.sVar.tsv")
mTSV <- cwlParam(baseCommand = mergeTSV,
                 inputs = InputParamList(m1, m2),
                 outputs = OutputParamList(out1))

s9 <- Step(id = "mergeTSV", run = mTSV,
           In = list(esnv = "Wrapper/EnsSNV",
                     eindel = "Wrapper/EnsINDEL"))
#' @include pl_neusomatic.R
s10 <- Step(id = "neusomaticPL", run = neusomatic,
            In = list(tbam = "tbam",
                      nbam = "nbam",
                      ref = "Ref",
                      region = "interval",
                      ensemble = "mergeTSV/tsv",
                      threads = "threads",
                      ovcf = list(valueFrom = "$(inputs.tbam.nameroot.split('_')[0])_neusomatic.vcf")))

o1a <- OutputParam(id = "mutect2filterVCF", type = "File", outputSource = "Mutect2PL/filterVCF")
o1b <- OutputParam(id = "mutect2passVCF", type = "File", outputSource = "Mutect2PL/passVCF")
o1c <- OutputParam(id = "mutect2conTable", type = "File", outputSource = "Mutect2PL/conTable")
o1d <- OutputParam(id = "mutect2artTable", type = "File", outputSource = "Mutect2PL/artTable")
o2 <- OutputParam(id = "MuSEout", type = "File", outputSource = "MuSE/outVcf")
o3a <- OutputParam(id = "strelka2snv", type = "File", outputSource = "mantaStrelka/snvs")
o3b <- OutputParam(id = "strelka2indel", type = "File", outputSource = "mantaStrelka/indels")
o4 <- OutputParam(id = "SomaticSniperout", type = "File", outputSource = "SomaticSniper/outVcf")
o5 <- OutputParam(id = "VarDictout", type = "File", outputSource = "VarDict/outVcf")
o6a <- OutputParam(id = "LoFreqsnp", type = "File", outputSource = "LoFreq/snp")
o6b <- OutputParam(id = "LoFreqindel", type = "File", outputSource = "LoFreq/indel")
o6c <- OutputParam(id = "LoFreqsnpdb", type = "File", outputSource = "LoFreq/snpdb")
o6d <- OutputParam(id = "LoFreqindeldb", type = "File", outputSource = "LoFreq/indeldb")
o7a <- OutputParam(id = "VarScanSnp", type = "File", outputSource = "VarScanPL/sSnp")
o7b <- OutputParam(id = "VarScanIndel", type = "File", outputSource = "VarScanPL/sIndel")
o7c <- OutputParam(id = "VarScansVcf", type = "File", outputSource = "VarScanPL/sVcf")
o8a <- OutputParam(id = "mergeTSVout", type = "File", outputSource = "mergeTSV/tsv")
o8b <- OutputParam(id = "WrapperSNV", type = "File", outputSource = "Wrapper/conSNV")
o8c <- OutputParam(id = "WrapperINDEL", type = "File", outputSource = "Wrapper/conINDEL")
o8d <- OutputParam(id = "WrapperESNV", type = "File", outputSource = "Wrapper/EnsSNV")
o8e <- OutputParam(id = "WrapperEINDEL", type = "File", outputSource = "Wrapper/EnsINDEL")
o9 <- OutputParam(id = "neusomaticVCF", type = "File", outputSource = "neusomaticPL/outVcf")

req1 <- list(class = "InlineJavascriptRequirement")
req2 <- list(class = "StepInputExpressionRequirement")
req3 <- list(class = "SubworkflowFeatureRequirement")
SomaticCallers <- cwlStepParam(requirements = list(req1, req2, req3),
                               inputs = InputParamList(p1, p2, p3, p6, p7,
                                                       p8, p9, p10, p11, p12, p13),
                               outputs = OutputParamList(o1a, o1b, o1c, o1d, o2,
                                                         o3a, o3b, o4, o5,
                                                         o6a, o6b, o6c, o6d,
                                                         o7a, o7b, o7c,
                                                         o8a, o8b, o8c, o8d, o8e,
                                                         o9))

SomaticCallers <- SomaticCallers + s1 + s2 + s3a + s3b + s3 + s4 +s5 + s6 + s7 + s8 + s9 + s10

