#' GATK alignment pipeline
#'
#' Workflows for processing high-throughput sequencing data for
#' variant discovery with GATK4 and related tools.
#' 
#' @format A `cwlStepParam` object.
#' \describe{
#'  \item{fq2ubam}{To covert fastq to ubam with read group information}
#'  \item{align}{To run BWA alignment and BAM BaseRecalibration.}
#' }
#' @source \url{https://github.com/gatk-workflows/seq-format-conversion}
#' @source \url{https://github.com/gatk-workflows/gatk4-data-processing}
#' @source \url{https://hubentu.github.io/others/Rcwl_GATK4.html}
#' @examples
#' ## scripts to build the pipeline
#' demo("GAlign")
"GAlign"
