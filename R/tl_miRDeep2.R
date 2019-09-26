## https://github.com/rajewsky-lab/mirdeep2
## miRDeep2.pl reads genome mappings miRNAs_ref/none miRNAs_other/none precursors/none [options] 2>report.log

p1 <- InputParam(id = "reads", type = "File", position = 1)
p2 <- InputParam(id = "genome", type = "File", position = 2)
p3 <- InputParam(id = "mappings", type = "File", position = 3)
p4 <- InputParam(id = "miRef", type = list("File", "string"),
                 position = 4, default = "none")
p5 <- InputParam(id = "miOther", type = list("File", "string"),
                 position = 5, default = "none")
p6 <- InputParam(id = "precursors", type = list("File", "string"),
                 position = 6, default = "none")
p7 <- InputParam(id = "species", type = "string",
                 prefix = "-t", position = 7)
o1 <- OutputParam(id = "csvfiles", type = "File[]", glob = "*.csv")
o2 <- OutputParam(id = "htmls", type = "File[]", glob = "*.html")
o3 <- OutputParam(id = "bed", type = "File", glob = "*.bed")
o4 <- OutputParam(id = "expression", type = "Directory", glob = "expression_analyses")
o5 <- OutputParam(id = "mirna_results", type = "Directory", glob = "mirna_results*")
o6 <- OutputParam(id = "pdfs", type = "Directory", glob = "pdf*")
req1 <- list(class = "DockerRequirement",
             dockerPull = "hubentu/mirdeep2")
miRDeep2 <- cwlParam(baseCommand = "miRDeep2.pl",
                     requirements = list(req1),
                     inputs = InputParamList(p1, p2, p3, p4, p5, p6, p7),
                     outputs = OutputParamList(o1, o2, o3, o4, o5, o6))
