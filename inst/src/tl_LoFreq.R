## https://csb5.github.io/lofreq/commands/#somatic
p1 <- InputParam(id = "tbam", type = "File", prefix = "-t", secondaryFiles = ".bai")
p2 <- InputParam(id = "nbam", type = "File", prefix = "-n", secondaryFiles = ".bai")
p3 <- InputParam(id = "ref", type = "File", prefix = "-f", secondaryFiles = ".fai")
p4 <- InputParam(id = "region", type = "File", prefix = "-l")
p5 <- InputParam(id = "dbsnp", type = "File", prefix = "-d", secondaryFiles = ".tbi")
p6 <- InputParam(id = "out", type = "string", prefix = "-o")
p7 <- InputParam(id = "threads", type = "int", prefix = "--threads")
o1 <- OutputParam(id = "snp", type = "File",
                  glob = "$(inputs.out)somatic_final.snvs.vcf.gz")
o2 <- OutputParam(id = "snpdb", type = "File",
                  glob = "$(inputs.out)somatic_final_minus-dbsnp.snvs.vcf.gz")
o3 <- OutputParam(id = "indel", type = "File",
                  glob = "$(inputs.out)somatic_final.indels.vcf.gz")
o4 <- OutputParam(id = "indeldb", type = "File",
                  glob = "$(inputs.out)somatic_final_minus-dbsnp.indels.vcf.gz")

req1 <- list(class = "DockerRequirement",
             dockerPull = "andreaswilm/lofreq:v2.1.2")

LoFreq <- cwlParam(baseCommand = c("lofreq", "somatic"),
                   requirements = list(req1),
                   inputs = InputParamList(p1, p2, p3, p4, p5, p6, p7),
                   outputs= OutputParamList(o1, o2, o3, o4))
