## http://varscan.sourceforge.net
p1 <- InputParam(id = "npileup", type = "File", position = 1)
p2 <- InputParam(id = "tpileup", type = "File", position = 2)
p3 <- InputParam(id = "bname", type = "string", position = 3)
p4 <- InputParam(id = "vcfout", type = "boolean", prefix = "--output-vcf",
                 position = 4, default = TRUE)
o1 <- OutputParam(id = "snp", type = "File", glob = "$(inputs.bname).snp.vcf")
o2 <- OutputParam(id = "indel", type = "File", glob = "$(inputs.bname).indel.vcf")

req1 <- list(class = "DockerRequirement",
             dockerPull = "mgibio/varscan-cwl:v2.4.2-samtools1.3.1")

VarScan2_somatic <- cwlParam(baseCommand = c("java", "-jar",
                                             "/opt/varscan/VarScan.jar", "somatic"),
                             requirements = list(req1),
                             inputs = InputParamList(p1, p2, p3, p4),
                             outputs = OutputParamList(o1, o2))
