
source("fastQC.R")
source("featureCounts.R")
source("samtools_flagstat.R")
source("samtools_index.R")
source("RSeQC.R")
source("RNASeq_Sf.R")
source("multiqc.R")
usethis::use_data(fastqc,
                  featureCounts,
                  samtools_flagstat,
                  samtools_index,
                  RSeQC,
                  rnaseq_Sf,
                  multiqc,
                  overwrite = TRUE)
