#' cwlTools
#'
#' To generate a file cache object for CWL tools in the package.
#'
#' @param cachePath The cache path of the BiocFileCache object.
#' @param ... options from `bfcadd`.
#' @importFrom rappdirs user_cache_dir
#' @import dplyr
#' @import BiocFileCache
#' @return A BiocFileCache object for existing CWL tools.
#' @export
#' @examples
#' tools <- cwlTools()

cwlTools <- function(cachePath = "Rcwl", ...) {
    f1 <- list.files(system.file("tools", package="RcwlPipelines"),
                     pattern = "*.R", recursive = TRUE,
                     full.names = TRUE)
    f2 <- list.files(system.file("pipelines", package="RcwlPipelines"),
                     pattern = "*.R", recursive = TRUE,
                     full.names = TRUE)
    fpath <- c(f1, f2)
    if(!grepl("^/", cachePath)){
        cachePath <- user_cache_dir(cachePath)
    }
    cwlBFC <- BiocFileCache(cachePath, ask = FALSE)
    fpath0 <- bfcinfo(cwlBFC) %>% pull(fpath)
    fpath <- fpath[!fpath %in% fpath0]
    if(length(fpath) > 0){
        rname <- sub(".R", "", basename(fpath))
        for(i in seq_along(fpath)){
            rname <- sub(".R$", "", basename(fpath[i]))
            bfcadd(cwlBFC, rname, fpath = fpath[i],
                   rtype = "local", action = "asis", ...)
        }
    }

    return(cwlBFC)
}
