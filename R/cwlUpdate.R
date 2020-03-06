#' cwlUpdate
#'
#' Function to update the Rcwl recipes from the RcwlRecipes github
#' repository.
#' @param cachePath The cache path of the BiocFileCache object to
#'     store the Rcwl tools and pipelines recipes.
#' @param force Whether to clean existing recipes cache.
#' @importFrom rappdirs user_cache_dir
#' @import BiocFileCache
#' @import utils
#' @export
#' @examples
#' \dontrun{
#' tools <- cwlUpdate()
#' }
cwlUpdate <- function(cachePath = "Rcwl", force = FALSE) {
    if(!file.exists(cachePath) & !grepl("^/", cachePath)){
        cachePath <- user_cache_dir(cachePath)
    }
    cwlBFC <- BiocFileCache(cachePath, ask = FALSE)
    
    ## req <- GET("https://api.github.com/repos/hubentu/RcwlRecipes/git/trees/master?recursive=1")
    ## stop_for_status(req)
    ## filelist <- unlist(lapply(content(req)$tree, "[", "path"), use.names = F)
    ## filelist <- filelist[grep("Rcwl/", filelist)]
    ## fpath <- paste0("https://raw.githubusercontent.com/hubentu/RcwlRecipes/master/", filelist)

    if(force){
        message("Warning: existing caches will be removed")
        bfcremove(cwlBFC, bfcinfo(cwlBFC)$rid)
    }

    message("Update scripts...")
    download.file("https://github.com/hubentu/RcwlRecipes/archive/master.zip",
                  file.path(cachePath, "master.zip"))
    unzip(file.path(cachePath, "master.zip"), exdir = cachePath)
    fpath <- list.files(file.path(cachePath, "RcwlRecipes-master/Rcwl"),
                        full.names = TRUE)
    
    ## fpath <- setdiff(fpath, bfcinfo(cwlBFC)$fpath)
    if("cwlMeta" %in% bfcmetalist(cwlBFC) & !force){
        BM <- bfcmeta(cwlBFC, "cwlMeta")
    }else{
        BM <- data.frame()
    }
    if(length(fpath) > 0){
        for(i in seq_along(fpath)){
            rname <- sub(".R$", "", basename(fpath[i]))
            ## add BFC
            if(!rname %in% bfcinfo(cwlBFC)$rname){
                add1 <- bfcadd(cwlBFC, rname, fpath = fpath[i],
                               rtype = "local", action = "asis")
                rid <- names(add1)
            }else{
                rid <- bfcrid(cwlBFC)[bfcinfo(cwlBFC)$rname == rname]
            }
            ## if(fpath[i] %in% bfcinfo(cwlBFC)$fpath){
            ##     rid  <- unique(bfcinfo(cwlBFC)$rid[bfcinfo(cwlBFC)$fpath == fpath[i]])
            ##     cwlBFC <- bfcupdate(cwlBFC, rid, fpath = fpath[i], ask = F)
            ## }else{
            ##     ## add1 <- bfcadd(cwlBFC, rname, fpath = fpath[i],
            ##     ##                rtype = "web")
            ##     add1 <- bfcadd(cwlBFC, rname, fpath = fpath[i],
            ##                    rtype = "local", action = "asis")
            ##     rid <- names(add1)
            ## }
            ## collect meta
            if(grepl("^tl_", rname)){
                rname <- sub("^tl_", "", rname)
                Type <- "tool"
                source(bfcrpath(cwlBFC)[rid], environment())
                cwl <- get(rname, environment())
                if(is(baseCommand(cwl), "function")){
                    Command <- "R function"
                }else{
                    Command <- paste(baseCommand(cwl), collapse = " ")
                }
                Container <- unlist(requirements(cwl))["dockerPull"]
                if(is.null(Container)) Container <- NA
            }else{
                rname <- sub("^pl_", "", rname)
                Type <- "pipeline"
                ss <- readLines(bfcrpath(cwlBFC)[rid])
                ss <- gsub("\\s*", "", ss[grep("run", ss)])
                Command <- paste(sub(".*run=(.*),", "\\1", ss), collapse = " + ")
                Container <- NA
            }
            bm <- data.frame(rid = rid, Type, Command, Container,
                             stringsAsFactors = FALSE)
            idx <- match(bm$rid, BM$rid)
            if(!is.na(idx)){
                BM[idx,] <- bm
            }else{
                BM <- rbind(BM, bm)
            }
        }
    }
    BM <- BM[match(bfcrid(cwlBFC), BM$rid),]
    bfcmeta(cwlBFC, "cwlMeta", overwrite = TRUE) <- BM
    return(cwlBFC)
}
