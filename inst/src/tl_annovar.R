
p1 <- InputParam(id = "vcf", type = "File", position = 1)
p2 <- InputParam(id = "db", type = "Directory", position = 2)
p3 <- InputParam(id = "build", type = "string", prefix = "-buildver", default = "hg19")
p4 <- InputParam(id = "aout", type = "string", prefix = "-out")
p5 <- InputParam(id = "protocol", type = "string",
                 prefix = "-protocol", default = "refGene,cosmic70")
p6 <- InputParam(id = "operation", type = "string",
                 prefix = "-operation", default = "g,f")
p7 <- InputParam(id = "nastring", type = "string",
                 prefix = "-nastring", default = ".")
o1 <- OutputParam(id = "Aout", type = "File",
                  glob = "$(inputs.aout).$(inputs.build)_multianno.vcf",
                  secondaryFiles = "$(inputs.aout).$(inputs.build)_multianno.txt")

req1 <- list(class = "InlineJavascriptRequirement")
req2 <- list(class = "DockerRequirement",
             dockerPull = "bioinfochrustrasbourg/annovar")

annovar <- cwlParam(baseCommand = "table_annovar.pl",
                    requirements = list(req1, req2),
                    arguments = list("-vcfinput"),
                    inputs = InputParamList(p1, p2, p3, p4, p5, p6, p7),
                    outputs = OutputParamList(o1))
