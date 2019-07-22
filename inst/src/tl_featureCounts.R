## featureCounts
## Note: output to current directory
f1 <- InputParam(id = "gtf", type = "File", prefix = "-a")
f2 <- InputParam(id = "count", type = "string", prefix = "-o")
f3 <- InputParam(id = "bam", type = "File")
o1 <- OutputParam(id = "Count", type = "File", glob = "$(inputs.count)", secondaryFiles = ".summary")

req1 <- list(class = "DockerRequirement",
             dockerPull = "hubentu/rcwl-rnaseq")
featureCounts <- cwlParam(baseCommand = "featureCounts",
                          requirements = list(req1),
                          inputs = InputParamList(f1, f2, f3),
                          outputs = OutputParamList(o1))

