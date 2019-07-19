
mc3 <- readCWL(system.file("mc3", "mc3_vcf2maf_full.cwl", package="RcwlPipelines"))
## fix bugs for vcf2maf pipeline and remove markfile step
steps(mc3) <- steps(mc3)[1:2]
mc3@inputs <- inputs(mc3)[1:10]
mc3@outputs$outmaf@outputSource <- "convert/outmaf"
outlist <- mc3@outputs
outlist <- list(outmaf=outlist[[1]], outvcf=OutputParam(id = "vepvcf", type = "File", outputSource = "convert/vepvcf"))
mc3@outputs@listData <- outlist

mc3@inputs$normal@secondaryFiles <- ".bai"
mc3@inputs$tumor@secondaryFiles <- ".bai"
mc3@inputs$refFasta@secondaryFiles <- list(".fai", ".gzi")
mc3@inputs$dbsnp@secondaryFiles <- ".tbi"
mc3@inputs$cosmic@secondaryFiles <- ".tbi"

mc3@steps$convert@Out <- list("outmaf", "vepvcf")

