
#' @include tl_miRMapper.R tl_miRDeep2.R
p1 <- InputParam(id = "reads", type = "File")
p2 <- InputParam(id = "format", type = "string", default = "-c")
p3 <- InputParam(id = "adapter", type = "string")
p4 <- InputParam(id = "len", type = "int", default = 18L)
p5 <- InputParam(id = "genome", type = "File",
                 valueFrom = "$(self.dirname + '/' + self.nameroot)",
                 secondaryFiles = c("$(self.nameroot + '.1.ebwt')",
                                    "$(self.nameroot + '.2.ebwt')",
                                    "$(self.nameroot + '.3.ebwt')",
                                    "$(self.nameroot + '.4.ebwt')",
                                    "$(self.nameroot + '.rev.1.ebwt')",
                                    "$(self.nameroot + '.rev.2.ebwt')"))
p6 <- InputParam(id = "miRef", type = list("File", "string"))
p7 <- InputParam(id = "miOther", type = list("File", "string"))
p8 <- InputParam(id = "precursors", type = list("File", "string"))
p9 <- InputParam(id = "species", type = "string")

s1 <- Step(id = "Mapper", run = miRMapper,
           In = list(reads = "reads",
                     format = "format",
                     adapter = "adapter",
                     genome = "genome",
                     preads = list(
                         valueFrom = "$(inputs.reads.nameroot)_collapsed.fa"),
                     arf = list(
                         valueFrom = "$(inputs.reads.nameroot)_collapsed.arf")))
s2 <- Step(id = "miRDeep2", run = miRDeep2,
           In = list(reads = "Mapper/pReads",
                     genome = "genome",
                     mappings = "Mapper/Arf",
                     miRef = "miRef",
                     miOther = "miOther",
                     precursors = "precursors",
                     species = "species"))

o1 <- OutputParam(id = "csvfiles", type = "File[]", outputSource = "miRDeep2/csvfiles")
o2 <- OutputParam(id = "htmls", type = "File[]", outputSource = "miRDeep2/htmls")
o3 <- OutputParam(id = "bed", type = "File", outputSource = "miRDeep2/bed")
o4 <- OutputParam(id = "expression", type = "Directory", outputSource = "miRDeep2/expression")
o5 <- OutputParam(id = "mirna_results", type = "Directory", outputSource = "miRDeep2/mirna_results")
o6 <- OutputParam(id = "pdfs", type = "Directory", outputSource = "miRDeep2/pdfs")
o7 <- OutputParam(id = "preads", type = "File", outputSource = "Mapper/pReads")
o8 <- OutputParam(id = "arf", type = "File", outputSource = "Mapper/Arf")

req1 <- list(class = "StepInputExpressionRequirement")
req2 <- list(class = "InlineJavascriptRequirement")
miRDeep2PL <- cwlStepParam(
    requirements = list(req1, req2),
    inputs = InputParamList(p1, p2, p3, p4, p5, p6, p7, p8, p9),
    outputs = OutputParamList(o1, o2, o3, o4, o5, o6, o7, o8))
miRDeep2PL <- miRDeep2PL + s1 + s2
