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
    f1 <- list.files(system.file("src", package="RcwlPipelines"),
                     pattern = "^tl_", full.names = TRUE)
    f2 <- list.files(system.file("src", package="RcwlPipelines"),
                     pattern = "^pl_", full.names = TRUE)
    fpath <- c(f1, f2)
    if(!grepl("^/", cachePath)){
        cachePath <- user_cache_dir(cachePath)
    }
    cwlBFC <- BiocFileCache(cachePath, ask = FALSE)
    fpath0 <- bfcinfo(cwlBFC) %>% pull(fpath)
    fpath <- fpath[!fpath %in% fpath0]
    if(length(fpath) > 0){
        for(i in seq_along(fpath)){
            rname <- sub(".R$", "", basename(fpath[i]))
            if(grepl("^tl_", rname)){
                rname <- sub("^tl_", "", rname)
                Type <- "tool"
                cwl <- get(rname, environment())
                Command <- paste(baseCommand(cwl), collapse = " ")
                Container <- unlist(requirements(cwl))["dockerPull"]
                if(is.null(Container)) Container <- NA
            }else{
                rname <- sub("^pl_", "", rname)
                Type <- "Pipeline"
                Command <- Container <- NA
            }

            add1 <- bfcadd(cwlBFC, rname, fpath = fpath[i],
                           rtype = "local", action = "asis", ...)
            bfcmeta(cwlBFC, "cwlMeta", append = TRUE) <-
                data.frame(rid = names(add1), Type, Command, Container,
                           stringsAsFactors = FALSE)
        }
    }

    return(cwlBFC)
}
