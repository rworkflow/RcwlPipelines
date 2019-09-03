## http://varscan.sourceforge.net
p1 <- InputParam(id = "vcf", type = "File", position = 1)
p2 <- InputParam(id = "indel", type = "File", prefix = "--indel-file", position = 2)
p3 <- InputParam(id = "outvcf", type = "string", prefix = "--output-file", position = 3)
o1 <- OutputParam(id = "outVcf", type = "File", glob = "$(inputs.outvcf)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "mgibio/varscan-cwl:v2.4.2-samtools1.3.1")
VarScan2_somaticFilter <- cwlParam(baseCommand = c("java", "-jar",
                                                   "/opt/varscan/VarScan.jar", "somaticFilter"),
                                    requirements = list(req1),
                                    inputs = InputParamList(p1, p2, p3),
                                    outputs = OutputParamList(o1))
