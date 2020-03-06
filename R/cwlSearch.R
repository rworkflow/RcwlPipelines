#' cwlSearch
#'
#' Function to search Rcwl tools and pipelines.
#' @param keyword The keyword to search.
#' @param bfc The BiocFileCache object for the recipes.
#' @param ... More options from `bfcquery`.
#' @return A BiocFileCache tibble.
#' @export
#' @examples
#' \dontrun{
#' cwlSearch("bwa")
#' }
cwlSearch <- function(keyword, bfc = NULL, ...){
    if(is.null(bfc)){
        cachePath <- user_cache_dir("Rcwl")
        bfc <- BiocFileCache(cachePath, ask = FALSE)
    }
    bfcquery(bfc, query = keyword,
             field = c("rname", "rpath", "fpath", "Command", "Container"),
             ...)
}
