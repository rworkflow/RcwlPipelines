#' @include tl_manta.R tl_strelka.R
p1 <- InputParam(id = "tbam", type = "File", secondaryFiles = ".bai")
p2 <- InputParam(id = "nbam", type = "File", secondaryFiles = ".bai")
p3 <- InputParam(id = "ref", type = "File", secondaryFiles = ".fai")
p4 <- InputParam(id = "region", type = "File", secondaryFiles = ".tbi")

s1 <- Step(id = "manta", run = manta,
           In = list(tbam = "tbam",
                     nbam = "nbam",
                     ref = "ref",
                     callRegions = "region"))
s2 <- Step(id = "strelka", run = strelka,
           In = list(tbam = "tbam",
                     nbam = "nbam",
                     ref = "ref",
                     callRegions = "region",
                     indelCandidates = "manta/candidateSmallIndels"))

o1 <- OutputParam(id = "snvs", type = "File", outputSource = "strelka/snvs")
o2 <- OutputParam(id = "indels", type = "File", outputSource = "strelka/indels")

mantaStrelka <- cwlStepParam(inputs = InputParamList(p1, p2, p3, p4),
                             outputs = OutputParamList(o1, o2))
mantaStrelka <- mantaStrelka + s1 + s2

mantaStrelka$ref <- "/rpcc/bioinformatics/reference/current/human_g1k_v37.fa"
mantaStrelka$tbam <- "/mnt/lustre/users/qhu/TESLA/output/BAM/Patient20_tumor/Patient20_tumor.bam"
mantaStrelka$nbam <- "/mnt/lustre/users/qhu/TESLA/output/BAM/Patient20_normal/Patient20_normal.bam"
mantaStrelka$region <- "/mnt/lustre/users/qhu/TESLA/output/neusomatic/Patient20_tumor/1/1.test.bed.gz"
