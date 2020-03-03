## https://cnvkit.readthedocs.io/en/stable/

p1 <- InputParam(id = "tbams", type = "File[]?", secondaryFiles = ".bai")
p2 <- InputParam(id = "ref", type = "File?", prefix = "--fasta",
                 secondaryFiles = ".fai")
p3 <- InputParam(id = "outdir", type = "string", prefix = "--output-dir")
p4 <- InputParam(id = "normal", type = "File[]?", prefix = "--normal",
                 secondaryFiles = ".bai")
p5a <- InputParam(id = "outref", type = "string?", prefix = "--output-reference")
p5b <- InputParam(id = "reference", type = "File?", prefix = "-r")
p6 <- InputParam(id = "target", type = "File?", prefix = "--targets")
p7 <- InputParam(id = "anti", type = "File?", prefix = "--antitargets")
p8 <- InputParam(id = "access", type = "File?", prefix = "--access")
p9 <- InputParam(id = "annotate", type = "File?", prefix = "--annotate")
p10 <- InputParam(id = "parallel", type = "int",
                 prefix = "-p", default = 1)
p11 <- InputParam(id = "diagram", type = "boolean", prefix = "--diagram",
                  default = TRUE)
p12 <- InputParam(id = "scatter", type = "boolean", prefix = "--scatter",
                  default = TRUE)
o1 <- OutputParam(id = "Outdir", type = "Directory", glob = "$(inputs.outdir)")
o2 <- OutputParam(id = "outRef", type = "File?", glob = "$(inputs.outref)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "etal/cnvkit")
cnvkit_batch <- cwlParam(baseCommand = c("cnvkit.py", "batch"),
                         requirements = list(req1),
                         inputs = InputParamList(p1, p2, p3, p4, p5a, p5b, p6, p7, p8, p9, p10, p11, p12),
                         outputs = OutputParamList(o1, o2))

