#' TCGA MC3 pipeline
#'
#' The MC3 somatic variant calling pipeline.
#' @format A `cwlStepParam` object.
#' \describe{
#'  \item{call_variants}{To call somatic variants by multiple pipelines}
#'  \item{convert}{To convert and merge VCFs to MAF}
#' }
#' @source \url{https://github.com/OpenGenomics/mc3}
#' @source \url{https://hubentu.github.io/others/Rcwl_MC3.html}
#' @examples
#' ## scripts to build the pipeline
#' demo("mc3")
"mc3"
