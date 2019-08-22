## https://github.com/Illumina/manta
p1 <- InputParam(id = "tbam", type = "File", prefix = "--tumorBam",
                 secondaryFiles = ".bai", position = 1)
p2 <- InputParam(id = "nbam", type = "File", prefix = "--normalBam",
                 secondaryFiles = ".bai", position = 2)
p3 <- InputParam(id = "ref", type = "File", prefix = "--referenceFasta",
                 secondaryFiles = ".fai", position = 3)
p4 <- InputParam(id = "callRegions", type = "File?", prefix = "--callRegions",
                 secondaryFiles = ".tbi", position = 4)
o1 <- OutputParam(id = "outDir", type = "Directory", glob = "mantaRunDir")

req1 <- list(class = "DockerRequirement",
             dockerPull = "cmopipeline/strelka2_manta")
req2 <- list(class = "ShellCommandRequirement")
manta <- cwlParam(baseCommand = "configManta.py",
                  requirements = list(req1, req2),
                  arguments = list(
                      "--runDir", "mantaRunDir", "--exome",
                      list(valueFrom = " && ", position = 5L),
                      list(valueFrom = "mantaRunDir/runWorkflow.py",
                           position = 6L),
                      list(valueFrom = "-m", position = 7L),
                      list(valueFrom = "local", position = 8L)),
                  inputs = InputParamList(p1, p2, p3, p4),
                  outputs = OutputParamList(o1))
