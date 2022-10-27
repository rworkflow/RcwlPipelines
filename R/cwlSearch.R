#' cwlSearch
#'
#' Function to search Rcwl tools and pipelines.
#' @param keyword A (vector of) character string as keywords to search
#'     for tools or pipelines.  Will be used to match patterns
#'     (case-insensitive) against `rname`, `rpath`, `fpath`, `Command`
#'     and `Container` column in the `bfc` object. 
#' @param bfc The `BiocFileCache` object for the recipes returned from
#'     `cwlUpdate`. The default is NULL which automatically detect the
#'     "Rcwl" cache directory.
#' @param type The `Type` to filter the results, "pipeline" or "tool".
#' @param ... More options from the internal `bfcquery` function.
#' @return A BiocFileCache tibble.
#' @export
#' @examples
#' \dontrun{
#' tls <- cwlSearch(c("bwa", "mem"))
#' data.frame(tls)
#' }

cwlSearch <- function(keyword, bfc = NULL, type = NULL, ...){
    bfcpath <- Sys.getenv("cachePath")
    if(bfcpath != ""){
        cachePath <- file.path(bfcpath, "Rcwl")
        bfc <- BiocFileCache(cachePath, ask = FALSE)
    }else if(is.null(bfc)){
        cachePath <- user_cache_dir("Rcwl")
        bfc <- BiocFileCache(cachePath, ask = FALSE)
    }
    res <- bfcquery(bfc, query = keyword,
                    field = c("rname", "rpath", "fpath", "Command", "Container"),
                    ignore.case = TRUE, ...)
    if(!is.null(type)){
        res <- res[res$Type == type,]
    }
    cwlHub(bfc[res$rid])
}
