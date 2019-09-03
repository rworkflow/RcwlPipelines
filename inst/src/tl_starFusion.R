## STAR Fusion
## https://github.com/STAR-Fusion/STAR-Fusion/wiki
p1 <- InputParam(id = "fq1", type = "File", prefix = "--left_fq")
p2 <- InputParam(id = "fq2", type = "File?", prefix = "--right_fq")
p3 <- InputParam(id = "genomedir", type = "Directory", prefix = "--genome_lib_dir")
p4 <- InputParam(id = "odir", type = "string", prefix = "--output_dir")
p5 <- InputParam(id = "cpu", type = "int", prefix = "--CPU")
o1 <- OutputParam(id = "sout", type = "Directory", glob = "$(inputs.odir)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "trinityctat/ctatfusion")
starFusion <- cwlParam(baseCommand = "/usr/local/src/STAR-Fusion/STAR-Fusion",
                       requirements = list(req1),
                       inputs = InputParamList(p1, p2, p3, p4, p5),
                                              outputs = OutputParamList(o1))
