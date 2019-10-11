##2. GenomicsDB + PoN, run with cwltool argument: --relax-path-checks
#' @include tl_GenomicsDB.R tl_PoN.R
p1 <- InputParam(id = "nvcf", type = InputArrayParam(items = "File"), secondaryFiles = ".idx")
p2 <- InputParam(id = "Ref", type = "File",
                 secondaryFiles = c(".fai", "$(self.nameroot).dict"))
p3 <- InputParam(id = "interval", type = "File")
p4 <- InputParam(id = "pvcf", type = "string")
p5 <- InputParam(id = "gresource", type = "File?", secondaryFiles = ".idx")

s1 <- Step(id = "GenomicsDB", run = GenomicsDB,
           In = list(vcf = "nvcf",
                     Ref = "Ref",
                     intervals = "interval"))
s2 <- Step(id = "PoN", run = PoN,
           In = list(db = "GenomicsDB/dbout",
                     Ref = "Ref",
                     pon = "pvcf",
                     gresource = "gresource"))

o1 <- OutputParam(id = "Pvcf", type = "File", outputSource = "PoN/pout")

GPoN <- cwlStepParam(inputs = InputParamList(p1, p2, p3, p4, p5),
                     outputs = OutputParamList(o1))
GPoN <- GPoN + s1 + s2
