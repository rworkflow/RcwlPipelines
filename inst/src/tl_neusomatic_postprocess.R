
p1 <- InputParam(id = "ref", type = "File",
                 prefix = "--reference", secondaryFiles = ".fai")
p2 <- InputParam(id = "tbam", type = "File",
                 prefix = "--tumor_bam", secondaryFiles = ".bai")
p3 <- InputParam(id = "pred", type = "File", prefix = "--pred_vcf")
p4 <- InputParam(id = "fcandidates", type = "File", prefix = "--candidates_vcf")
p5 <- InputParam(id = "ensemble", type = "File", prefix = "--ensemble_tsv")
p6 <- InputParam(id = "ovcf", type = "string", prefix = "--output_vcf")
o1 <- OutputParam(id = "oVcf", type = "File", glob = "$(inputs.ovcf)")

req1 <- list(class = "DockerRequirement",
             dockerPull = "msahraeian/neusomatic")
neusomatic_postprocess <- cwlParam(baseCommand = c("python",
                                                   "/opt/neusomatic/neusomatic/python/postprocess.py"),
                                   requirements = list(req1),
                                   arguments = list("--work", "."),
                                   inputs = InputParamList(p1, p2, p3, p4, p5, p6),
                                   outputs = OutputParamList(o1))
