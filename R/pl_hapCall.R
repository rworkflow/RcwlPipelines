##source(system.file("tools", "runWDL.R", package = "RcwlPipelines"))
#' @include tl_runWDL.R tl_mvOut.R
## haplotypecaller
rscripts <- "args <- commandArgs(TRUE)
bam <- args[1]
intervals <- args[2]
tmpl3 <- args[3]
json1 <- jsonlite::fromJSON(tmpl3)
json1$HaplotypeCallerGvcf_GATK4.input_bam <- bam
json1$HaplotypeCallerGvcf_GATK4.input_bam_index <- sub('.bam', '.bai', bam)
json1$HaplotypeCallerGvcf_GATK4.scattered_calling_intervals_list <- intervals
cat(jsonlite::toJSON(json1, pretty = TRUE, auto_unbox = T))"
rscripts <- gsub("\n", ";", rscripts)

p1 <- InputParam(id = "bam", type = "string")
p2 <- InputParam(id = "intervals", type = "string")
p3 <- InputParam(id = "tmpl", type = "File")
o1 <- OutputParam(id = "json", type = "File", glob = "tmpl3.json")
hapJson <- cwlParam(baseCommand = c("Rscript", "-e", rscripts),
                    inputs = InputParamList(p1, p2, p3),
                    outputs = OutputParamList(o1),
                    stdout = "tmpl3.json")

p1 <- InputParam(id = "bam", type = "string")
p2 <- InputParam(id = "intervals", type = "string")
p3 <- InputParam(id = "cromwell", type = "File")
p4 <- InputParam(id = "wdl", type = "File")
p5 <- InputParam(id = "tmpl", type = "File")
o1 <- OutputParam(id = "hclog", type = "File", outputSource = "HC/log")
o2 <- OutputParam(id = "outdir", type = "Directory", outputSource = "mvOut/OutDir")
hapCall <- cwlStepParam(inputs = InputParamList(p1, p2, p3, p4, p5),
                        outputs = OutputParamList(o1, o2))
s1 <- Step(id = "hapJson", run = hapJson,
           In = list(bam = "bam",
                     intervals = "intervals",
                     tmpl = "tmpl"))
s2 <- Step(id = "HC", run = runWDL,
           In = list(cromwell = "cromwell",
                     wdl = "wdl",
                     json = "hapJson/json"))
s3 <- Step(id = "mvOut", run = mvOut,
           In = list(logFile = "HC/log"))

hapCall <- hapCall + s1 + s2 + s3
