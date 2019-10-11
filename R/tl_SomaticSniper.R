## https://github.com/genome/somatic-sniper
p1 <- InputParam(id = "ref", type = "File", prefix = "-f", secondaryFiles = ".fai", position = 1)
p2 <- InputParam(id = "tbam", type = "File", secondaryFiles = ".bai", position = 2)
p3 <- InputParam(id = "nbam", type = "File", secondaryFiles = ".bai", position = 3)
p4 <- InputParam(id = "vcf", type = "string", position = 4)
o1 <- OutputParam(id = "outVcf", type = "File", glob = "$(inputs.vcf)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "lethalfang/somaticsniper:1.0.5.0-2")

SomaticSniper <- cwlParam(baseCommand = "/opt/somatic-sniper/build/bin/bam-somaticsniper",
                          requirements = list(req1),
                          arguments = list("-q", "10", "-F", "vcf"),
                          inputs = InputParamList(p1, p2 , p3, p4),
                          outputs = OutputParamList(o1))
