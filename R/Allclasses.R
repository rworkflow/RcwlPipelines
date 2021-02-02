setClass("cwlHub", contains = "BiocFileCacheReadOnly")
#' cwlHub
#'
#' The `cwlHub` constructor for `BiocFileCache` object.
#' @param BFC A BiocFileCache created for `RcwlRecipes`.
#' @importFrom S4Vectors DataFrame
#' @export
cwlHub <- function(BFC){
    if(class(BFC) == "BiocFileCache"){
        BFC <- new("BiocFileCacheReadOnly",
                   cache = BFC@cache,
                   rid = bfcrid(BFC))
    }
    new("cwlHub", BFC)
}

#' mcols
#'
#' DataFrame information from the `BicFileCache` object.
#' @param x A `cwlHub` object
#' @importFrom S4Vectors mcols get_showHeadLines get_showTailLines
#' @exportMethod mcols
setMethod("mcols", "cwlHub", function(x){
    mc <- bfcinfo(x)[bfcrid(x) %in% x@rid,]
    DataFrame(mc)
})

setMethod("show", "cwlHub", function(object){
    rid <- object@rid
    mc <- bfcinfo(object)[bfcrid(object) %in% rid,]        

    cat("cwlHub with", length(rid), "records\n")
    cat("cache path: ", bfccache(object), "\n")
    mdate <- tail(sort(as.Date(mc$mtime)), 1)
    cat("# last modified date: ", as.character(mdate), "\n")
    cat("# cwlSearch() to query scripts\n")
    cat("# cwlLoad('title') to load the script\n")
    cat("# additional mcols(): rid, rpath, Type, Container, mtime, ...\n")
    ## https://github.com/Bioconductor/AnnotationHub/blob/master/R/Hub-class.R#L602
    .some <-
        function(elt, nh, nt, fill="...", width=getOption("width") - 13L)
    {
        answer <- if (length(elt) < nh + nt + 1L)
                      elt
                  else
                      c(head(elt, nh), fill, tail(elt, nt))
        ifelse(nchar(answer) > width,
               sprintf("%s...", substring(answer, 1L, width-3L)),
               answer)
    }
    nhead <- get_showHeadLines()
    ntail <- get_showTailLines()
    rownames <- paste0("  ", .some(rid, nhead, ntail))
    out <- matrix(c(.some(rep("|", length(rid)), nhead, ntail, fill=""),
                    .some(mc$rname, nhead, ntail),
                    .some(mc$Command, nhead, ntail)),
                  ncol=3L,
                  dimnames=list(rownames, c("", "title", "Command")))
    cat("\n")
    print(out, quote=FALSE, right=FALSE)
})

#' extract
#' 
#' @rdname cwlHub-methods
#' @param x A `cwlHub` object.
#' @param value The "BFC" ID to extract the subset.
#' @exportMethod [
setMethod("[", "cwlHub", function(x, value){
    idx <- match(value, x@rid)
    isNA <- is.na(idx)
    x@rid[idx[!isNA]] <- value[!isNA]
    x
})

#' title
#'
#' @rdname cwlHub-methods
#' @param object A `cwlHub` object.
#' @export
#' @examples
#' \dontrun{
#' tools <- cwlUpdate()
#' title(tools)
#' }
title <- function(object){
    mcols(object)$rname
}

#' Command
#'
#' @rdname cwlHub-methods
#' @export
Command <- function(object){
    mcols(object)$Command
}

#' Container
#'
#' @rdname cwlHub-methods
#' @export
Container <- function(object){
    mcols(object)$Container
}

#' Type
#'
#' @rdname cwlHub-methods
#' @export
Type <- function(object){
    mcols(object)$Type
}
