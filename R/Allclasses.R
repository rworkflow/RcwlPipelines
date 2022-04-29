#' cwlHub
#'
#' `cwlHub` class, constructor, and methods. 
#' @rdname cwlHub-class 
#' @exportClass cwlHub
#' 
setClass("cwlHub", contains = "BiocFileCacheReadOnly")

#' @rdname cwlHub-class
#' @param BFC A BiocFileCache created for `RcwlRecipes`.
#' @return cwlHub: a `cwlHub` object with slots of `rid` and `cache` path. 
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

#' @rdname cwlHub-class
#' @param x A `cwlHub` object
#' @return mcols: a `DataFrame` with information from the `BicFileCache` object.
#' @importFrom S4Vectors mcols get_showHeadLines get_showTailLines
#' @exportMethod mcols
#' 
setMethod("mcols", "cwlHub", function(x){
    mc <- bfcinfo(x)[bfcrid(x) %in% x@rid,]
    DataFrame(mc)
})

#' @rdname cwlHub-class
#' @param object A `cwlHub` object
#' @exportMethod show
#' 
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
    if (length(rid) > 0) {
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
    }
})

#' @rdname cwlHub-class
#' @param x A `cwlHub` object.
#' @param value The "BFC" ID to extract the subset.
#' @return [: a subset of `cwlHub` records.
#' @exportMethod [
#' 
setMethod("[", "cwlHub", function(x, value){
    idx <- match(value, x@rid)
    isNA <- is.na(idx)
    x@rid[idx[!isNA]] <- value[!isNA]
    x
})

#' @rdname cwlHub-class
#' @return title: the `Rcwl` recipe names for tools or pipelines.
#' @export
#' @examples
#' \dontrun{
#' tools <- cwlUpdate()
#' t1 <- tools["BFC178"]
#' title(t1)
#' Command(t1)
#' Container(t1)
#' Type(t1)
#' }
#' 
title <- function(object){
    mcols(object)$rname
}

#' @rdname cwlHub-class
#' @return Command: The name of `Rcwl` wrapped command line tools.
#' @export
Command <- function(object){
    mcols(object)$Command
}

#' @rdname cwlHub-class
#' @return Container: the container name for the `Rcwl` recipe if exist. Otherwise `NA`. 
#' @export
Container <- function(object){
    mcols(object)$Container
}

#' @rdname cwlHub-class
#' @return Type: The type of the `Rcwl` recipe, "pipeline" or "tool".
#' @export
Type <- function(object){
    mcols(object)$Type
}
