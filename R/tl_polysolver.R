## https://software.broadinstitute.org/cancer/cga/polysolver_run
## run with --no-read-only
p1 <- InputParam(id = "bam", type = "File",
                 position = 1, secondaryFiles = ".bai")
p2 <- InputParam(id = "race", type = "string",
                 position = 2, default = "Unknown")
p3 <- InputParam(id = "includeFreq", type = "int",
                 position = 3, default = 1L)
p4 <- InputParam(id = "build", type = "string",
                 position = 4, default = "hg19")
p5 <- InputParam(id = "format", type = "string",
                 position = 5, default = "STDFQ")
p6 <- InputParam(id = "insertCalc", type = "int",
                 position = 6, default = 0L)
p7 <- InputParam(id = "outDir", type = "Directory",
                 position = 7)
o1 <- OutputParam(id = "hla", type = "File", glob = "*.hla.txt")
req1 <- list(class = "DockerRequirement",
             dockerPull = "sachet/polysolver:v4")

polysolver <- cwlParam(baseCommand = c("bash", "/home/polysolver/scripts/shell_call_hla_type"),
                       requirements = list(req1),
                      arguments = list(
                          list(valueFrom = "$(runtime.outdir)", position = 7L)
                      ),
                      inputs = InputParamList(p1, p2, p3, p4, p5, p6),
                      outputs = OutputParamList(o1))
