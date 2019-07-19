##source(system.file("tools", "runWDL.R", package = "RcwlPipelines"))
#' @include tl_runWDL.R tl_mvOut.R
## joint discovery
rscripts <- "args <- commandArgs(TRUE)
splitList <- function(s)as.list(unlist(strsplit(s, split = ',')))
sampleName <- args[1]
gvcf <- args[2]
callsetName <- args[3]
intervals <- args[4]
unpadded_intervals <- args[5]
tmpl4 <- args[6]
json1 <- jsonlite::fromJSON(tmpl4, simplifyVector=FALSE)
json1$JointGenotyping.sample_names <- splitList(sampleName)
json1$JointGenotyping.input_gvcfs <- splitList(gvcf)
json1$JointGenotyping.input_gvcfs_indices <- splitList(gsub('gz', 'gz.tbi', gvcf))
json1$JointGenotyping.callset_name <- callsetName
json1$JointGenotyping.eval_interval_list <- intervals
json1$JointGenotyping.unpadded_intervals_file <- unpadded_intervals
cat(jsonlite::toJSON(json1, pretty = TRUE, auto_unbox = T))"
rscripts <- gsub("\n", ";", rscripts)

p1 <- InputParam(id = "sampleName", type = "string", position = 1)
p2 <- InputParam(id = "gvcf", type = "string", position = 2)
p3 <- InputParam(id = "callsetName", type = "string", position = 3)
p4 <- InputParam(id = "intervals", type = "string", position = 4)
p5 <- InputParam(id = "unpadded_intervals", type = "string", position =5)
p6 <- InputParam(id = "tmpl", type = "File", position = 6)
o1 <- OutputParam(id = "json", type = "File", glob = "tmpl4.json")
jdJson <- cwlParam(baseCommand = c("Rscript", "-e", rscripts),
                    inputs = InputParamList(p1, p2, p3, p4, p5, p6),
                    outputs = OutputParamList(o1),
                    stdout = "tmpl4.json")

p7 <- InputParam(id = "cromwell", type = "File")
p8 <- InputParam(id = "wdl", type = "File")
o1 <- OutputParam(id = "hclog", type = "File", outputSource = "JD/log")
o2 <- OutputParam(id = "outdir", type = "Directory", outputSource = "mvOut/OutDir")
jdCall <- cwlStepParam(inputs = InputParamList(p1, p2, p3, p4, p5, p6, p7, p8),
                        outputs = OutputParamList(o1, o2))
s1 <- Step(id = "jdJson", run = jdJson,
           In = list(sampleName = "sampleName",
                     gvcf = "gvcf",
                     callsetName = "callsetName",
                     intervals = "intervals",
                     unpadded_intervals = "unpadded_intervals",
                     tmpl = "tmpl"))
s2 <- Step(id = "JD", run = runWDL,
           In = list(cromwell = "cromwell",
                     wdl = "wdl",
                     json = "jdJson/json"))
s3 <- Step(id = "mvOut", run = mvOut,
           In = list(logFile = "JD/log"))

jdCall <- jdCall + s1 + s2 + s3
