#' GATK haplotypecaller pipeline
#'
#' The workflow runs HaplotypeCaller from GATK4 in GVCF mode on a
#' single sample according to the GATK Best Practices (June 2016),
#' scattered across intervals. The workflow from github,
#' gatk4-germline-snps-indels (last update: 7/23/2018) was cloned to
#' this package.
#' 
#' @format A `cwlStepParam` object.
#' \describe{
#'  \item{HC}{HaplotypeCaller from GATK4}
#' }
#' @source \url{https://github.com/gatk-workflows/gatk4-germline-snps-indels}
#' @source \url{https://hubentu.github.io/others/Rcwl_GATK4.html}
#' @examples
#' ## scripts to build the pipeline
#' ## demo("hapCall")
"hapCall"
