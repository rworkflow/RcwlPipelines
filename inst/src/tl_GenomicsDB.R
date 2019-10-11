## Create a GenomicsDB
p1 <- InputParam(id = "vcf", type = InputArrayParam(prefix = "-V", items = "File"), secondaryFiles = ".idx")
p2 <- InputParam(id = "Ref", prefix = "-R", type = "File", secondaryFiles = c(".fai", "$(self.nameroot).dict"))
p3 <- InputParam(id = "db", prefix = "--genomicsdb-workspace-path", type = "string", default = "pon_db")
p4 <- InputParam(id = "intervals", prefix = "-L", type = "File")
o1 <- OutputParam(id = "dbout", type = "Directory", glob = "$(inputs.db)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "broadinstitute/gatk:latest")

GenomicsDB <- cwlParam(baseCommand = c("gatk", "GenomicsDBImport"),
                       requirements = list(req1),
                       arguments = list("--merge-input-intervals"),
                       inputs = InputParamList(p1, p2, p3, p4),
                       outputs = OutputParamList(o1))
