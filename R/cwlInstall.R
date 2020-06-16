.sourceCWL <- function(rscript, env = .GlobalEnv){
    .env <- new.env()
    source(rscript, .env)
    objs <- ls(.env)
    oidx <- sapply(objs,
                   function(x)is(get(x, envir = .env), "cwlParam"))
    for(i in seq(sum(oidx))){
        assign(objs[oidx][i],
               get(objs[oidx][i], envir = .env),
               envir = env)
    }
}

#' cwlInstall
#' 
#' To source Rcwl scripts
#' @param file The Rcwl tool or pipeline receipes. The dependent
#'     Rcwl scripts should be included with @include tag.
#' @param rname The `rname` to install from the `bfc` object.
#' @param bfc The BiocFileCache object for the recipes.
#' @param env The enviroment to export.
#' @import methods
#' @export
#' @examples
#' \dontrun{
#' cwlInstall(rname = "tl_bwa")
#' }
cwlInstall <- function(file, rname, bfc = NULL, env = .GlobalEnv) {
    if(is.null(bfc)){
        cachePath <- user_cache_dir("Rcwl")
        bfc <- BiocFileCache(cachePath, ask = FALSE)
    }
    if(missing(file)){
        file <- bfcrpath(bfc)[bfcinfo(bfc)$rname == rname]
    }
    if(!file.exists(file)){
        rname = file
    }
    if(!missing(rname)){
        file <- bfcrpath(bfc)[bfcinfo(bfc)$rname == rname]
    }
    scripts <- readLines(file)
    iscripts <- grep("@include", scripts, value = TRUE)
    if(length(iscripts) > 0){
        rscripts <- grep(".R$",
                         unlist(strsplit(iscripts, split = " ")),
                         value = TRUE)
        if(length(rscripts) > 0){
            sapply(rscripts, function(x){
                rscript <- file.path(dirname(file), x)
                if(any(grepl("cwlStepParam", readLines(rscript)))){
                    cwlInstall(file = rscript, bfc = bfc, env = env)
                }else{
                    .sourceCWL(rscript, env)
                }
            })
        }
    }
    .sourceCWL(file, env)
}
