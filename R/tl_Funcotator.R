## Funcotator
p1 <- InputParam(id = "vcf", type = "File", prefix = "-V")
p2 <- InputParam(id = "ref", prefix = "-R", type = "File",
                 secondaryFiles = c(".fai", "$(self.nameroot).dict"))
p3 <- InputParam(id = "outf", type = "string", prefix = "--output-file-format", default = "MAF")
p4 <- InputParam(id = "dsource", type = "Directory", prefix = "--data-sources-path")
p5 <- InputParam(id = "version", type = "string", prefix = "--ref-version", default = "hg19")
p6 <- InputParam(id = "maf", type = "string", prefix = "-O")
o1 <- OutputParam(id = "mout", type = "File", glob = "$(inputs.maf)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "broadinstitute/gatk:latest")

Funcotator <- cwlParam(baseCommand = c("gatk",
                                       "Funcotator"),
                       requirements = list(req1),
                       arguments = list("--remove-filtered-variants"),
                       inputs = InputParamList(p1, p2, p3, p4, p5, p6),
                       outputs = OutputParamList(o1))
