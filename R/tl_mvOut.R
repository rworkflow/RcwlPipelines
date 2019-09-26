## glob outputs
## rscripts <- "args <- commandArgs(TRUE)
## logFile <- args[1]
## log1 <- readLines(logFile)
## startn <- grep('Final Outputs:', log1)+1
## endn <- grep('}$', log1)
## endn <- endn[endn > startn][1]
## logOut <- jsonlite::fromJSON(log1[startn:endn])
## logOut <- logOut[lengths(logOut)>0]
## dir.create('output', showWarnings = FALSE)
## lapply(logOut, function(x)file.rename(x, file.path('output', basename(x))))"
## rscripts <- gsub("\n", "; ", rscripts)

## p1 <- InputParam(id = "logFile", type = "File")
## o1 <- OutputParam(id = "OutDir", type = "Directory", glob = "output")
## mvOut <- cwlParam(baseCommand = c("Rscript", "-e", rscripts),
##                   inputs = InputParamList(p1),
##                   outputs = OutputParamList(o1))
## #' @importFrom jsonlite fromJSON
mvout <- function(logFile){
    log1 <- readLines(logFile)
    startn <- grep('Final Outputs:', log1)+1
    endn <- grep('}$', log1)
    endn <- endn[endn > startn][1]
    logOut <- jsonlite::fromJSON(log1[startn:endn])
    logOut <- logOut[lengths(logOut)>0]
    dir.create('output', showWarnings = FALSE)
    lapply(logOut, function(x)file.rename(x, file.path('output', basename(x))))
}

p1 <- InputParam(id = "logFile", type = "File", prefix = "logFile=", separate = F)
o1 <- OutputParam(id = "OutDir", type = "Directory", glob = "output")
mvOut <- cwlParam(baseCommand = mvout,
                  inputs = InputParamList(p1),
                  outputs = OutputParamList(o1))
