#' GATK joint discovery pipeline
#'
#' The joint discovery and VQSR filtering portion of the GATK Best
#' Practices (June 2016) for germline SNP and Indel discovery in human
#' whole-genome sequencing (WGS) and exome sequencing data.
#'
#' @format A `cwlStepParam` object.
#' \describe{
#'  \item{JD}{variant joint genotyping}
#' }
#' @source \url{https://github.com/gatk-workflows/gatk4-germline-snps-indels}
#' @source \url{https://hubentu.github.io/others/Rcwl_GATK4.html}
#' @examples
#' ## scripts to build the pipeline
#' demo("jdCall")
"jdCall"
