##source(system.file("tools", "runWDL.R", package = "RcwlPipelines"))
#' @include tl_mvOut.R tl_runWDL.R
## prepare Json for fq2ubam
rscripts <- "args <- commandArgs(TRUE)
splitList <- function(s)as.list(unlist(strsplit(s, split = ',')))
template <- args[1]
input1 <- jsonlite::fromJSON(template, simplifyVector=FALSE)
input1$ConvertPairedFastQsToUnmappedBamWf.readgroup_name <- splitList(args[2])
input1$ConvertPairedFastQsToUnmappedBamWf.sample_name <- splitList(args[3])
input1$ConvertPairedFastQsToUnmappedBamWf.fastq_1 <- splitList(args[4])
input1$ConvertPairedFastQsToUnmappedBamWf.fastq_2 <- splitList(args[5])
input1$ConvertPairedFastQsToUnmappedBamWf.ubam_list_name <- splitList(args[3])[[1]][1]
input1$ConvertPairedFastQsToUnmappedBamWf.library_name <- splitList(args[6])
input1$ConvertPairedFastQsToUnmappedBamWf.platform_unit <- splitList(args[7])
input1$ConvertPairedFastQsToUnmappedBamWf.run_date <- list(rep(as.character(Sys.Date()), lengths(splitList(args[2]))))
input1$ConvertPairedFastQsToUnmappedBamWf.platform_name <- splitList(args[8])
input1$ConvertPairedFastQsToUnmappedBamWf.sequencing_center <- splitList(args[9])
cat(jsonlite::toJSON(input1, pretty = TRUE, auto_unbox = T))"
rscripts <- gsub("\n", "; ", rscripts)

p1 <- InputParam(id = "tmpl", type = "File", position = 1)
p2 <- InputParam(id = "fastq1", type = "string", position = 4)
p3 <- InputParam(id = "fastq2", type = "string", position = 5)
p4 <- InputParam(id = "readGroup", type = "string", position = 2)
p5 <- InputParam(id = "sampleName", type = "string", position = 3)
p6 <- InputParam(id = "library", type = "string", position = 6)
p7 <- InputParam(id = "platunit", type = "string", position = 7)
p8 <- InputParam(id = "platform", type = "string", position = 8)
p9 <- InputParam(id = "center", type = "string", position = 9)
o1 <- OutputParam(id = "jsonOut", type = "File", glob = "tmpl1.json")
fq2ubamJson <- cwlParam(baseCommand = c("Rscript", "-e", rscripts),
                   inputs = InputParamList(p1, p2, p3, p4, p5, p6, p7, p8, p9),
                   outputs = OutputParamList(o1),
                   stdout = "tmpl1.json")

## prepare json for ubam to bam
rscripts <- "args <- commandArgs(TRUE)
fqlog <- args[1]
tmpl <- args[2]
log1 <- readLines(fqlog)
startn <- grep('Final Outputs:', log1)+1
endn <- grep('}$', log1)[1]
ubamOut <- jsonlite::fromJSON(log1[startn:endn])
sampleName <- sub('.list', '', basename(ubamOut$ConvertPairedFastQsToUnmappedBamWf.unmapped_bam_list))
json1 <- jsonlite::fromJSON(tmpl)
json1$PreProcessingForVariantDiscovery_GATK4.sample_name <- sampleName
json1$PreProcessingForVariantDiscovery_GATK4.flowcell_unmapped_bams_list <- ubamOut$ConvertPairedFastQsToUnmappedBamWf.unmapped_bam_list
cat(jsonlite::toJSON(json1, pretty = TRUE, auto_unbox = T))"
rscripts <- gsub("\n", "; ", rscripts)

p1 <- InputParam(id = "fqlog", type = "File", position = 1)
p2 <- InputParam(id = "template", type = "File", position = 2)
o1 <- OutputParam(id = "json", type = "File", glob = "temp2.json")
ubam2bamJson <- cwlParam(baseCommand = c("Rscript", "-e", rscripts),
                       inputs = InputParamList(p1, p2),
                       outputs = OutputParamList(o1),
                       stdout = "temp2.json")

## GATK alignment pipeline
p1 <- InputParam(id = "fastq1", type = "string")
p2 <- InputParam(id = "fastq2", type = "string")
p3 <- InputParam(id = "readGroup", type = "string")
p4 <- InputParam(id = "sampleName", type = "string")
p5 <- InputParam(id = "library", type = "string")
p6 <- InputParam(id = "platunit", type = "string")
p7 <- InputParam(id = "platform", type = "string")
p8 <- InputParam(id = "center", type = "string")
p9 <- InputParam(id = "tmpl1", type = "File")
p10 <- InputParam(id = "wdl1", type = "File")
p11 <- InputParam(id = "tmpl2", type = "File")
p12 <- InputParam(id = "wdl2", type = "File")
p13 <- InputParam(id = "cromwell", type = "File")
o1 <- OutputParam(id = "bamlog", type = "File", outputSource = "align/log")
o2 <- OutputParam(id = "outdir", type = "Directory", outputSource = "mvOut/OutDir")

GAlign <- cwlStepParam(inputs = InputParamList(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13),
                         outputs = OutputParamList(o1, o2))
s1 <- Step(id = "fqJson", run = fq2ubamJson,
           In = list(tmpl = "tmpl1",
                     fastq1 = "fastq1",
                     fastq2 = "fastq2",
                     readGroup = "readGroup",
                     sampleName = "sampleName",
                     library = "library",
                     platunit = "platunit",
                     platform = "platform",
                     center = "center"))
s2 <- Step(id = "fq2ubam", run = runWDL,
           In = list(cromwell = "cromwell",
                     wdl = "wdl1",
                     json = "fqJson/jsonOut"))
s3 <- Step(id = "ubam2bamJson", run = ubam2bamJson,
           In = list(fqlog = "fq2ubam/log",
                     template = "tmpl2"))
s4 <- Step(id = "align", run = runWDL,
           In = list(cromwell = "cromwell",
                     wdl = "wdl2",
                     json = "ubam2bamJson/json"))
s5 <- Step(id = "mvOut", run = mvOut,
           In = list(logFile = "align/log"))

GAlign <- GAlign + s1 + s2 + s3 + s4 + s5
