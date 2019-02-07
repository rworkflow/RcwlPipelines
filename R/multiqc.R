#' multiqc
#'
#' The mutliQC command tool
#' @format A `cwlParam` object.
#' \describe{
#'  \item{bwaAlign}{to align with bwa and sort with samtools}
#'  \item{mergeBamDup}{to merge BAMs from different flowcells and then mark duplicates with picard}
#' }
#' @source \url{https://hubentu.github.io/others/Rcwl_DNASeq_Align.html}
#' @source \url{https://multiqc.info}
#' @examples
#' ## scripts to build the pipeline
#' dir(system.file("scripts/DNAalign", package="RcwlPipelines"))
"multiqc"
