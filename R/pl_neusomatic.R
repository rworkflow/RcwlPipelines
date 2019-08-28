
#' @include tl_neusomatic_preprocess.R tl_neusomatic_call.R tl_neusomatic_postprocess.R
p1 <- InputParam(id = "tbam", type = "File", secondaryFiles = ".bai")
p2 <- InputParam(id = "nbam", type = "File", secondaryFiles = ".bai")
p3 <- InputParam(id = "ref", type = "File", secondaryFiles = ".fai")
p4 <- InputParam(id = "region", type = "File")
p5 <- InputParam(id = "ensemble", type = "File")
p6 <- InputParam(id = "threads", type = "int", default = 2L)
p7 <- InputParam(id = "ovcf", type = "string")
o1 <- OutputParam(id = "outVcf", type = "File", outputSource = "postprocess/oVcf")

s1 <- Step(id = "preprocess", run = neusomatic_preprocess,
           In = list(tbam = "tbam",
                     nbam = "nbam",
                     ref = "ref",
                     ensemble = "ensemble",
                     region = "region",
                     threads = "threads"))

s2 <- Step(id = "call", run = neusomatic_call,
           In = list(candidates = "preprocess/candidates",
                     ref = "ref"))

s3 <- Step(id = "postprocess", run = neusomatic_postprocess,
           In = list(ref = "ref",
                     tbam = "tbam",
                     pred = "call/pred",
                     fcandidates = "preprocess/fcandidates",
                     ensemble = "ensemble",
                     ovcf = "ovcf"))

neusomatic <- cwlStepParam(inputs = InputParamList(p1, p2, p3, p4, p5, p6, p7),
                           outputs = OutputParamList(o1))

neusomatic <- neusomatic + s1 + s2 + s3
