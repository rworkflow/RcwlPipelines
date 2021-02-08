#' cwlUpdate
#'
#' Function to sync and get the most updated Rcwl recipes from the
#' RcwlRecipes github 
#' @param cachePath The cache path of the BiocFileCache object to
#'     store the Rcwl tools and pipelines recipes.
#' @param force Whether to clean existing recipes cache.
#' @param branch The branch of github recipes repository. It can be
#'     "master" and "dev". "force = TRUE" is recommended when swithing
#'     branch.
#' @importFrom rappdirs user_cache_dir
#' @import BiocFileCache
#' @import utils
#' @export
#' @examples
#' \dontrun{
#' tools <- cwlUpdate()
#' }
cwlUpdate <- function(cachePath = "Rcwl", force = FALSE, branch = NULL) {
    if(is.null(branch) & grepl("alpha|unstable", version$status)){
        branch <- "dev"
    }else if(is.null(branch)){
        branch <- "master"
    }
    
    if(!file.exists(cachePath) & !grepl("^/", cachePath)){
        cachePath <- user_cache_dir(cachePath)
    }
    cwlBFC <- BiocFileCache(cachePath, ask = FALSE)

    if(force){
        message("Warning: existing caches will be removed")
        bfcremove(cwlBFC, bfcinfo(cwlBFC)$rid)
    }

    message("Update scripts...")
    download.file(paste0("https://github.com/hubentu/RcwlRecipes/archive/", branch, ".zip"),
                  file.path(cachePath, paste0(branch, ".zip")))
    unzip(file.path(cachePath, paste0(branch, ".zip")), exdir = cachePath)
    fpath <- list.files(file.path(cachePath, paste0("RcwlRecipes-", branch, "/Rcwl")),
                        full.names = TRUE)
    
    if(length(fpath) > 0){
        rnames <- sub(".R$", "", basename(fpath))
        ex <- rnames %in% bfcinfo(cwlBFC)$rname
        if(sum(!ex)>0){
            idx <- which(!ex)
            for(i in idx){
                add1 <- bfcadd(cwlBFC, rnames[i], fpath = fpath[i],
                               rtype = "local", action = "asis")
                message(basename(add1), " added")
            }
        }
    }
    meta <- read.csv(file.path(cachePath, paste0("RcwlRecipes-", branch, "/cwlMeta.csv")), row.names = 1)
    BM <- data.frame(rid = bfcinfo(cwlBFC)$rid,
                     meta[match(bfcinfo(cwlBFC)$rname, rownames(meta)), ])
    bfcmeta(cwlBFC, "cwlMeta", overwrite = TRUE) <- BM
    return(cwlHub(cwlBFC))
}

cwlMeta <- function(fpaths){
    BM <- c()
    for(i in seq(length(fpaths))){
        fpath <- fpaths[i]
        rname <- sub(".R$", "", basename(fpath))
        if(grepl("^tl_", rname)){
            Type <- "tool"
            t1 <- cwlLoad(fpath)
            if(is(baseCommand(t1), "function")){
                Command <- "R function"
            }else{
                Command <- paste(baseCommand(t1), collapse = " ")
            }
            Container <- unlist(requirements(t1))["dockerPull"]
            if(is.null(Container)) Container <- NA
        }else{
            Type <- "pipeline"
            p1 <- cwlLoad(fpath)
            Command <- paste(names(runs(p1)), collapse = "+")
            Container <- NA
        }
        mtime <- file.info(fpath)$mtime
        bm <- data.frame(row.names = rname, Type, Command, Container, mtime,
                         stringsAsFactors = FALSE)
        BM <- rbind(BM, bm)
    }
    BM
}
