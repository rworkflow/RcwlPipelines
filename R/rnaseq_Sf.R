#' RNASeq pipeline with STAR and featureCounts
#'
#' RNASeq pipeline by STAR and featureCounts.
#' @format A `cwlStepParam` object.
#' \describe{
#'  An RNASeq alignment and quantification pipeline built by `Rcwl`, which contains steps:
#'  \item{fastqc}{The reads QC step by fastQC}
#'  \item{STAR}{The alignment step by STAR}
#'  \item{samtools_index}{Index bam file by samtools}
#'  \item{samtools_flagstat}{Flag stat by samtools}
#'  \item{featureCounts}{Gene level quantification by featureCounts}
#'  \item{RSeQC}{QC for RNASeq alignments by RSeQC}
#' }
#' @source \url{https://hubentu.github.io/others/Rcwl_RNASeq.html}
#' @examples
#' ## scripts to build the pipeline
#' demo("RNASeq_Sf")
"rnaseq_Sf"
