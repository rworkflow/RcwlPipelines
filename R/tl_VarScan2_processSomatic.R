## http://varscan.sourceforge.net
p1 <- InputParam(id = "vcf", type = "File")
o1 <- OutputParam(id = "somaticHC", type = "File",
                  glob = "$(inputs.vcf.nameroot).Somatic.hc.vcf")
o2 <- OutputParam(id = "somatic", type = "File",
                  glob = "$(inputs.vcf.nameroot).Somatic.vcf")
o3 <- OutputParam(id = "germline", type = "File",
                  glob = "$(inputs.vcf.nameroot).Germline.vcf")
o4 <- OutputParam(id = "germlineHC", type = "File",
                  glob = "$(inputs.vcf.nameroot).Germline.hc.vcf")
o5 <- OutputParam(id = "LOH", type = "File",
                  glob = "$(inputs.vcf.nameroot).LOH.vcf")
o6 <- OutputParam(id = "LOHHC", type = "File",
                  glob = "$(inputs.vcf.nameroot).LOH.hc.vcf")

req1 <- list(class = "DockerRequirement",
             dockerPull = "mgibio/varscan-cwl:v2.4.2-samtools1.3.1")
req2 <- list(class = "InitialWorkDirRequirement",
             listing = list("$(inputs.vcf)"))

VarScan2_processSomatic <- cwlParam(baseCommand = c("java", "-jar",
                                                    "/opt/varscan/VarScan.jar", "processSomatic"),
                                    requirements = list(req1, req2),
                                    inputs = InputParamList(p1),
                                    outputs = OutputParamList(o1, o2, o3, o4, o5, o6))

