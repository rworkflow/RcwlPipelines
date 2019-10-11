
p1 <- InputParam(id = "ref", type = "File",
                 prefix = "--genome-reference",
                 secondaryFiles = c(".fai", "$(self.nameroot).dict"))
p2 <- InputParam(id = "tbam", type = "File", prefix = "--tumor-bam",
                 secondaryFile = ".bai")
p3 <- InputParam(id = "nbam", type = "File", prefix = "--normal-bam",
                 secondaryFiles = ".bai")
p4 <- InputParam(id = "mutect2", type = "File?", prefix = "--mutect2")
p5 <- InputParam(id = "varscanSnv", type = "File?", prefix = "--varscan-snv")
p6 <- InputParam(id = "varscanIndel", type = "File?", prefix = "--varscan-indel")
p7 <- InputParam(id = "sniper", type = "File?", prefix = "--sniper")
p8 <- InputParam(id = "vardict", type = "File?", prefix = "--vardict")
p9 <- InputParam(id = "muse", type = "File?", prefix = "--muse")
p10 <- InputParam(id = "strelkaSnv", type = "File?", prefix = "--strelka-snv")
p11 <- InputParam(id = "strelkaIndel", type = "File?", prefix = "--strelka-indel")
p12 <- InputParam(id = "lofreqSnv", type = "File?", prefix = "--lofreq-snv")
p13 <- InputParam(id = "lofreqIndel", type = "File?", prefix = "--lofreq-indel")
p14 <- InputParam(id = "region", type = "File?", prefix = "--inclusion-region")
p15 <- InputParam(id = "dbsnp", type = "File",
                  prefix = "--dbsnp", secondaryFiles = ".tbi")
o1 <- OutputParam(id = "conSNV", type = "File", glob = "Consensus.sSNV.vcf")
o2 <- OutputParam(id = "conINDEL", type = "File", glob = "Consensus.sINDEL.vcf")
o3 <- OutputParam(id = "EnsSNV", type = "File", glob = "Ensemble.sSNV.tsv")
o4 <- OutputParam(id = "EnsINDEL", type = "File", glob = "Ensemble.sINDEL.tsv")

req1 <- list(class = "DockerRequirement",
             dockerPull = "lethalfang/somaticseq:2.7.2")
SomaticSeq_Wrapper <- cwlParam(baseCommand = "/opt/somaticseq/SomaticSeq.Wrapper.sh",
                               requirements = list(req1),
                               arguments = list("--output-dir", ".", "--gatk",
                                                "/opt/GATK/GenomeAnalysisTK.jar"),
                               inputs = InputParamList(p1, p2, p3, p4, p5, p6, p7, p8,
                                                       p9, p10, p11, p12, p13, p14, p15),
                               outputs = OutputParamList(o1, o2, o3, o4))
