## https://bioinformatics.mdanderson.org/public-software/muse/
p1 <- InputParam(id = "tbam", type = "File",
                 secondaryFiles = ".bai", position = 1)
p2 <- InputParam(id = "nbam", type = "File",
                 secondaryFiles = ".bai", position = 2)
p3 <- InputParam(id = "ref", type = "File", prefix = "-f",
                 secondaryFiles = ".fai", position = 3)
p4 <- InputParam(id = "region", type = "File?",
                 prefix = "-l", position = 4)
p5 <- InputParam(id = "dbsnp", type = "File", prefix = "-D",
                 secondaryFiles = ".tbi", position = 10)
p6 <- InputParam(id = "vcf", type = "string",
                 prefix = "-O", position = 11)
p7 <- InputParam(id = "exome", type = "boolean",
                 prefix = "-E", position = 12, default = TRUE)
p8 <- InputParam(id = "genome", type = "boolean",
                 prefix = "-G", position = 12, default = FALSE)
o1 <- OutputParam(id = "outVcf", type = "File", glob = "$(inputs.vcf)")

req1 <- list(class = "DockerRequirement",
             dockerPull = "marghoob/muse:1.0rc_c")
req2 <- list(class = "ShellCommandRequirement")
MuSE <- cwlParam(baseCommand = c("MuSEv1.0rc_submission_c039ffa", "call"),
                 requirements = list(req1, req2),
                 arguments = list(
                     "-O", "output",
                     list(valueFrom = " && ", position = 5L),
                     list(valueFrom = "MuSEv1.0rc_submission_c039ffa", position = 6L),
                     list(valueFrom = "sump", position = 7L),
                     list(valueFrom = "-I", position = 8L),
                     list(valueFrom = "output.MuSE.txt", position = 9L)
                     ),
                 inputs = InputParamList(p1, p2, p3, p4, p5, p6, p7, p8),
                 outputs = OutputParamList(o1))
