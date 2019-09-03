
p1 <- InputParam(id = "tbam", type = "File",
                 prefix = "--tumor_bam", secondaryFiles = ".bai")
p2 <- InputParam(id = "nbam", type = "File",
                 prefix = "--normal_bam", secondaryFiles = ".bai")
p3 <- InputParam(id = "ref", type = "File",
                 prefix = "--reference", secondaryFiles = ".fai")
p4 <- InputParam(id = "region", type = "File", prefix = "--region_bed")
p5 <- InputParam(id = "ensemble", type = "File", prefix = "--ensemble_tsv")
p6 <- InputParam(id = "mapq", type = "int", prefix = "--min_mapq", default = 10L)
p7 <- InputParam(id = "threads", type = "int", prefix = "--num_threads", default = 2L)
o1 <- OutputParam(id = "candidates", type = "File[]",
                  glob = "dataset/work*candidates*.tsv",
                  secondaryFiles = ".idx")
o2 <- OutputParam(id = "fcandidates", type = "File",
                  glob = "work_tumor/filtered_candidates.vcf")
req1 <- list(class = "DockerRequirement",
             dockerPull = "msahraeian/neusomatic")
req2 <- list(class = "ShellCommandRequirement")
neusomatic_preprocess <- cwlParam(baseCommand = c("python",
                                                  "/opt/neusomatic/neusomatic/python/preprocess.py"),
                                  requirements = list(req1, req2),
                                  arguments = list("--mode", "call", "--work", ".",
                                                   "--scan_alignments_binary",
                                                   "/opt/neusomatic/neusomatic/bin/scan_alignments",
                                                   list(valueFrom = "&& for i in `ls dataset/*/*tsv*`; do can=`echo $i | sed 's/\\/candidates/_candidates/'`; cp $i $can;done",
                                                        position = 10L, shellQuote = FALSE)),
                                  inputs = InputParamList(p1, p2, p3, p4, p5, p6, p7),
                                  outputs = OutputParamList(o1, o2))
