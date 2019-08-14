
## cwlFiles <- list.files(".", "tl_")
## for(ss in cwlFiles){
##     cwlEnv <- new.env()
##     source(ss, local = cwlEnv)
##     objs <- ls(envir = cwlEnv)
##     idx <- sapply(objs, function(x)is(get(x, cwlEnv), "cwlParam"))
##     cwl <- get(objs[idx], cwlEnv)
##     if(baseCommand(cwl)[1] == "Rscript"){
##         til <- paste("Rscript:", objs[idx])
##     }else if(all(baseCommand(cwl) == "java")){
##         til <- paste("java:", objs[idx])
##     }else{
##         til <- paste(baseCommand(cwl), collapse = " ")
##     }
##     exp <- c(paste("#'", til), "#' @export",
##              paste0("\"", objs[idx], "\""), "\n")
##     cat(exp, file = "cwl_tools.R", sep = "\n", append = TRUE)
## }

#' bcftools view
#' @export
"bcfview"


#' blastn
#' @export
"blastn"


#' bowtie2
#' @export
"bowtie2"


#' bowtie2-build
#' @export
"bowtie2_build"


#' bwa index
#' @export
"bwa_index"


#' bwa mem
#' @export
"bwa"


#' cutadapt
#' @export
"cutadapt"


#' fastqc
#' @export
"fastqc"


#' featureCounts
#' @export
"featureCounts"


#' gatk CalculateContamination
#' @export
"CalculateContamination"


#' gatk CollectSequencingArtifactMetrics
#' @export
"ColSeqArtifact"


#' gatk FilterMutectCalls
#' @export
"FilterMutectCalls"


#' gatk FilterByOrientationBias
#' @export
"FilterOBias"


#' gatk Funcotator
#' @export
"Funcotator"


#' gatk GenomicsDBImport
#' @export
"GenomicsDB"


#' gatk GetPileupSummaries
#' @export
"GetPileupSummaries"


#' gatk Mutect2
#' @export
"Mutect2"


#' gatk CreateSomaticPanelOfNormals
#' @export
"PoN"


#' geneBody_coverage.py
#' @export
"geneBody_coverage"


#' genePredToBed
#' @export
"genePredToBed"


#' gtfToGenePred
#' @export
"gtfToGenePred"


#' hisat2
#' @export
"hisat2_align"


#' hisat2-build
#' @export
"hisat2_build"


#' htseq-count
#' @export
"htseq"


#' makeblastdb
#' @export
"makeblastdb"


#' picard MarkDuplicates
#' @export
"markdup"


#' picard MergeSamFiles
#' @export
"mergeBam"


#' multiqc
#' @export
"multiqc"


#' Rscript: mvOut
#' @export
"mvOut"


#' read_distribution.py
#' @export
"read_distribution"


#' java: runWDL
#' @export
"runWDL"


#' salmon index
#' @export
"salmon_index"


#' salmon quant
#' @export
"salmon_quant"


#' samtools flagstat
#' @export
"samtools_flagstat"


#' samtools index
#' @export
"samtools_index"


#' samtools view
#' @export
"sam2bam"


#' samtools sort
#' @export
"sortBam"


#' STAR
#' @export
"STAR"

#' STAR-Fusion
#' @export
"starFusion"
