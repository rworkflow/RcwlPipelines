
mc3 <- readCWL(system.file("mc3/mc3_vcf2maf_full.cwl", package="RcwlPipelines"))
## fix bugs for vcf2maf pipeline and remove markfile step
steps(mc3) <- steps(mc3)[1:2]
mc3@inputs@inputs <- inputs(mc3)[1:10]
mc3@outputs@outputs$outmaf@outputSource <- "convert/outmaf"
outlist <- mc3@outputs@outputs
outlist <- SimpleList(outmaf=outlist[[1]], outvcf=OutputParam(id = "vepvcf", type = "File", outputSource = "convert/vepvcf"))
mc3@outputs@outputs <- outlist

mc3@inputs@inputs$normal@secondaryFiles <- ".bai"
mc3@inputs@inputs$tumor@secondaryFiles <- ".bai"
mc3@inputs@inputs$refFasta@secondaryFiles <- list(".fai", ".gzi")
mc3@inputs@inputs$dbsnp@secondaryFiles <- ".tbi"
mc3@inputs@inputs$cosmic@secondaryFiles <- ".tbi"

mc3@steps@steps$convert@Out <- list("outmaf", "vepvcf")
usethis::use_data(mc3, overwrite = TRUE)

##
mc3_variant <- readCWL(system.file("mc3/mc3_variant.cwl", package="RcwlPipelines"))
usethis::use_data(mc3_variant, overwrite = TRUE)
##
mc3_vcf2maf <- readCWL(system.file("mc3/mc3_vcf2maf.cwl", package="RcwlPipelines"))
usethis::use_data(mc3_vcf2maf)
