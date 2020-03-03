## https://github.com/timothyjamesbecker/SVE
p1 <- InputParam(id = "bam", type = "File",
                 secondaryFiles = ".bai", prefix = "-b")
p2 <- InputParam(id = "ref", type = "File",
                 secondaryFile = c(".fai", "$(self.nameroot).dict"), prefix = "-r")
p3 <- InputParam(id = "outdir", type = "string", prefix = "-o")
p4 <- InputParam(id = "tools", type = "string", prefix = "-s",
                 default = "breakdancer,cnmops,gatk_haplo,delly,lumpy,cnvnator,breakseq,tigra,genome_strip,hydra")
o1 <- OutputParam(id = "outs", type = "Directory", glob = "$(inputs.outdir)")

req1 <- list(class = "DockerRequirement",
             dockerPull = "timothyjamesbecker/sve")
SVE_VP <- cwlParam(baseCommand = "/software/SVE/scripts/variant_processor.py",
                   requirements = list(req1),
                   inputs = InputParamList(p1, p2, p3, p4),
                   outputs = OutputParamList(o1))
