#' seawrch containers
#'
#' To search container images for a tool in certain repository from
#' quay.io or dockerhub.
#' @param tool The tool to search.
#' @param repo The repository to lookup.
#' @param source The container server to search, quay.io or dockerhub.
#' @return A DataFrame contains image tag names, updated dates and
#'     image sizes.
#' @importFrom httr GET content
#' @importFrom S4Vectors DataFrame
#' @export
#' @examples
#' searchContainer("samtools")
searchContainer <- function(tool, repo = "biocontainers",
                            source = c("quay", "dockerhub")){
    source <- match.arg(source)
    switch(source,
           quay = searchQuay(tool, repo),
           dockerhub = searchDockerhub(tool, repo))
}

searchQuay <- function(tool, repo = "biocontainers"){
    api <- "https://quay.io/api/v1/repository"
    url <- paste(api, repo, tool, sep = "/")
    res <- GET(url)
    if(res$status_code == 200){
        cont <- content(res)
        tags <- do.call(rbind, cont$tags)
        tags <- apply(tags, 2, unlist)
        tags <- tags[, c("name", "last_modified", "size")]
        res <- cbind(tool, paste0("quay.io/", repo), tags)
        dates <- as.Date(tags[, "last_modified"], "%a, %d %b %Y %T")
        df <- DataFrame(res[order(dates, decreasing = TRUE), ])
        df$container <- paste0(df$V2, "/", df$tool, ":", df$name)
        return(df)
    }else{
        return(DataFrame())
    }
}

searchDockerhub <- function(tool, repo = "biocontainers"){
    api <- "https://registry.hub.docker.com/v2/repositories/"
    url <- paste(api, repo, tool, "tags", sep = "/")
    res <- GET(url)
    if(res$status_code == 200 & content(res)$count > 0){
        res <- content(res)
        res <- do.call(rbind, res$results)
        res <- res[, c("name", "last_updated", "full_size")]
        res <- apply(res, 2, unlist)
        res <- cbind(tool, repo, res)
        dates <- as.Date(res[,"last_updated"], "%Y-%m-%d")
        df <- DataFrame(res[order(dates, decreasing = TRUE),])
        df$container <- paste0(df$repo, "/", df$tool, ":", df$name)
        return(df)
    }else{
        return(DataFrame())
    }
}
