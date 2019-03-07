## gene body coverage
p1 <- InputParam(id = "bam", type = "File", prefix = "-i", secondaryFiles = ".bai")
p2 <- InputParam(id = "bed", type = "File", prefix = "-r")
p3 <- InputParam(id = "prefix", type = "string", prefix = "-o")
o1 <- OutputParam(id = "gCovPDF", type = "File", glob = "*.geneBodyCoverage.curves.pdf")
o2 <- OutputParam(id = "gCovTXT", type = "File", glob = "*.geneBodyCoverage.txt")
req1 <- list(class = "DockerRequirement",
             dockerPull = "hubentu/rcwl-rnaseq")
geneBody_coverage <- cwlParam(baseCommand = c("geneBody_coverage.py"),
                              requirements = list(req1),
                              inputs = InputParamList(p1, p2, p3),
                              outputs = OutputParamList(o1, o2))
