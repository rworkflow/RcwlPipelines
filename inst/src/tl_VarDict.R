## https://github.com/AstraZeneca-NGS/VarDictJava

p1 <- InputParam(id = "tbam", type = "File",
                 secondaryFiles = ".bai", position = -1)
p2 <- InputParam(id = "nbam", type = "File",
                 secondaryFiles = ".bai", position = -1)
p3 <- InputParam(id = "ref", type = "File", prefix = "-G",
                 secondaryFiles = ".fai", position = 1)
p4 <- InputParam(id = "region", type = "File")
p5 <- InputParam(id = "af", type = "float", default = 0.05, position = -1)
p6 <- InputParam(id = "vcf", type = "string", position = -1)
o1 <- OutputParam(id = "outVcf", type = "stdout")
req1 <- list(class = "DockerRequirement",
             dockerPull = "lethalfang/vardictjava:1.5.1")
req2 <- list(class = "ShellCommandRequirement")
VarDict <- cwlParam(baseCommand = c("/opt/VarDict-1.5.1/bin/VarDict"),
                    requirements = list(req1, req2),
                    arguments = list(
                        list(valueFrom = "-b", position = 2L),
                        list(valueFrom = "$(inputs.tbam.path)|$(inputs.nbam.path)", position = 3L),
                        list(valueFrom = "-f", position  = 4L),
                        list(valueFrom = "$(inputs.af)", position = 5L),
                        "-c", "1", "-S", "2", "-E", "3", "-g", "4",
                        list(valueFrom = " | ", position = 6L),
                        list(valueFrom = "/opt/VarDict/testsomatic.R", position = 7L),
                        list(valueFrom = " | ", position = 8L),
                        list(valueFrom = "/opt/VarDict/var2vcf_paired.pl", position = 9L),
                        list(valueFrom = "-N", position = 10L),
                        list(valueFrom = "TUMOR|NORMAL", position = 11L),
                        list(valueFrom = "-f", position = 12L),
                        list(valueFrom = "$(inputs.af)", position = 13L)
                    ),
                    inputs = InputParamList(p1, p2, p3, p4, p5, p6),
                    outputs = OutputParamList(o1),
                    stdout = "$(inputs.vcf)")
