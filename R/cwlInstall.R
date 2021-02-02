.sourceCWL <- function(rscript, env = .GlobalEnv){


    .env <- new.env()
    source(rscript, .env)
    objs <- ls(.env)
    oidx <- sapply(objs,
                   function(x)is(get(x, envir = .env), "cwlProcess"))
    for(i in seq(sum(oidx))){
        assign(objs[oidx][i],
               get(objs[oidx][i], envir = .env),
               envir = env)
    }
}

#' cwlInstall
#' 
#' To source Rcwl scripts
#' @param rname The name or filepath of tool or pipeline to install
#'     (`rname` or `fpath` column from the `bfc` object returned from
#'     `cwlSearch`).
#' @param bfc The `BiocFileCache` object for the recipes returned from
#'     `cwlUpdate`. The default is NULL which automatically detect the
#'     "Rcwl" cache directory.
#' @param env The R enviroment to export to. The default is
#'     `.GlobalEnv`.
#' @details Note to developers that the dependent Rcwl scripts should
#'     be included in the recipe with `@include` tag.
#' @import methods
#' @export
#' @examples
#' \dontrun{
#' tls <- cwlSearch("bwa")
#' tls$rname
#' cwlInstall("tl_bwa")
#' cwlInstall(tls$fpath[tls$rname == "tl_bwa"])  ## equivalent
#' bwa
#' }
cwlInstall <- function(rname, bfc = NULL, env = .GlobalEnv) {
    if(is.null(bfc)){
        cachePath <- user_cache_dir("Rcwl")
        bfc <- BiocFileCache(cachePath, ask = FALSE)
    }
    if (missing(rname))
        stop("Please provide a valid name or filepath for the tool/pipeline.")
    idx <- match(rname, bfcinfo(bfc)$rname)
    if (!is.na(idx)) {
        fpath <- bfcrpath(bfc)[idx]
    } else {
        if (file.exists(rname)){
            fpath <- rname
        } else {
            stop("Please provide a valid name or filepath for the tool/pipeline.")
        }
    }
    scripts <- readLines(fpath)
    iscripts <- grep("@include", scripts, value = TRUE)
    if(length(iscripts) > 0){
        rscripts <- grep(".R$",
                         unlist(strsplit(iscripts, split = " ")),
                         value = TRUE)
        if(length(rscripts) > 0){
            sapply(rscripts, function(x){
                rscript <- file.path(dirname(fpath), x)
                if(any(grepl("cwlWorkflow", readLines(rscript)))){
                    cwlInstall(rscript, bfc = bfc, env = env)
                }else{
                    .sourceCWL(rscript, env)
                }
            })
        }
    }
    .sourceCWL(fpath, env)
}
