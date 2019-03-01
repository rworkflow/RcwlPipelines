#' DNASeq alignment, merge and markduplicates
#'
#' The DNASeq pipeline to run bwa alignment, merge and mark
#' duplicates.
#' 
#' @format A `cwlStepParam` object.
#' \describe{
#'  \item{bwaAlign}{to align fastqs with bwa and sort with samtools}
#'  \item{mergeBamDup}{to merge BAMs from different flowcells
#' and then mark duplicates with picard}
#' }
#' @source \url{https://hubentu.github.io/others/Rcwl_DNASeq_Align.html}
#' @examples
#' ## scripts to build the pipeline
#' ## demo("alignMerge")
"alignMerge"
