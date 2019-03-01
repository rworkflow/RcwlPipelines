## RUN WDL
p1 <- InputParam(id = "cromwell", type = "File", position = 1, prefix = "-jar")
p2 <- InputParam(id = "run", type = "string", default = "run", position = 2)
p3 <- InputParam(id = "wdl", type = "File", position = 3)
p4 <- InputParam(id = "json", type = "File", position = 4, prefix = "-i")
o1 <- OutputParam(id = "log", type = "File", glob="$(inputs.wdl.basename).log")
runWDL <- cwlParam(baseCommand = "java",
                   inputs = InputParamList(p1, p2, p3, p4),
                   outputs = OutputParamList(o1),
                   stdout = "$(inputs.wdl.basename).log")

## glob outputs
rscripts <- "args <- commandArgs(TRUE)
logFile <- args[1]
log1 <- readLines(logFile)
startn <- grep('Final Outputs:', log1)+1
endn <- grep('}$', log1)
endn <- endn[endn > startn][1]
logOut <- jsonlite::fromJSON(log1[startn:endn])
logOut <- logOut[lengths(logOut)>0]
dir.create('output', showWarnings = FALSE)
lapply(logOut, function(x)file.rename(x, file.path('output', basename(x))))"
rscripts <- gsub("\n", "; ", rscripts)

p1 <- InputParam(id = "logFile", type = "File")
o1 <- OutputParam(id = "OutDir", type = "Directory", glob = "output")
mvOut <- cwlParam(baseCommand = c("Rscript", "-e", rscripts),
                  inputs = InputParamList(p1),
                  outputs = OutputParamList(o1))
