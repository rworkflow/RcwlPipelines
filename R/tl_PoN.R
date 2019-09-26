## CreateSomaticPanelOfNormals
p1 <- InputParam(id = "db", type = "Directory", prefix = "gendb://",
                 separate = F, position = 1)
p2 <- InputParam(id = "Ref", prefix = "-R", type = "File",
                 secondaryFiles = c(".fai", "$(self.nameroot).dict"), position = 2)
p3 <- InputParam(id = "pon", prefix = "-O", type = "string", position = 3)
p4 <- InputParam(id = "gresource", type = "File?", prefix = "--germline-resource",
                 secondaryFiles = ".idx", position = 4)
o1 <- OutputParam(id = "pout", type = "File", glob = "$(inputs.pon)",
                  secondaryFiles = ".idx")
req1 <- list(class = "DockerRequirement",
             dockerPull = "broadinstitute/gatk:latest")
## fix bug for lock files on POSIX filesystems
req2 <- list(class = "EnvVarRequirement",
             envDef = list("TILEDB_DISABLE_FILE_LOCKING" = "1"))

PoN <- cwlParam(baseCommand = c("gatk", "CreateSomaticPanelOfNormals"),
                requirements = list(req1, req2),
                arguments = list("-V"),
                inputs = InputParamList(p1, p2, p3, p4),
                outputs = OutputParamList(o1))
