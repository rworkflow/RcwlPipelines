## merge BAMs
p1 <- InputParam(id = "ibam", type = InputArrayParam(items = "File",
                                                    prefix = "I=",
                                                    separate = FALSE))
p2 <- InputParam(id = "obam", type = "string", prefix = "O=", separate = FALSE)
o1 <- OutputParam(id = "oBam", type = "File", glob = "$(inputs.obam)")
req1 <- list(class = "DockerRequirement",
             dockerPull = "biocontainers/picard:2.3.0")
mergeBam <- cwlParam(baseCommand = c("picard", "MergeSamFiles"),
                     requirements = list(req1),
                     inputs = InputParamList(p1, p2),
                     outputs = OutputParamList(o1))

## picard markdup
p1 <- InputParam(id = "ibam", type = "File", prefix = "I=", separate = FALSE)
p2 <- InputParam(id = "obam", type = "string", prefix = "O=", separate = FALSE)
p3 <- InputParam(id = "matrix", type = "string", prefix = "M=", separate = FALSE)
o1 <- OutputParam(id = "mBam", type = "File", glob = "$(inputs.obam)")
o2 <- OutputParam(id = "Mat", type = "File", glob = "$(inputs.matrix)")
markdup <- cwlParam(baseCommand = c("picard", "MarkDuplicates"),
                    requirements = list(req1),
                    inputs = InputParamList(p1, p2, p3),
                    outputs = OutputParamList(o1, o2))

## bam flagstat
p1 <- InputParam(id = "bam", type = "File")
o1 <- OutputParam(id = "flagstat", type = "File", glob = "$(inputs.bam.nameroot).flagstat.txt")
req2 <- list(class = "DockerRequirement",
             dockerPull = "biocontainers/samtools")
samtools_flagstat <- cwlParam(baseCommand = c("samtools", "flagstat"),
                              requirements = list(req2),
                              inputs = InputParamList(p1),
                              outputs = OutputParamList(o1),
                              stdout = "$(inputs.bam.nameroot).flagstat.txt")

## mergeBam + markdup + index + flagstat
p1 <- InputParam(id = "ibam", type = "File[]")
p2 <- InputParam(id = "obam", type = "string")
o1 <- OutputParam(id = "oBam", type = "File", outputSource = "markdup/mBam")
o2 <- OutputParam(id = "matrix", type = "File", outputSource = "markdup/Mat")
o3 <- OutputParam(id = "Idx", type = "File", outputSource = "samtools_index/idx")
o4 <- OutputParam(id = "stat", type = "File",
                  outputSource = "samtools_flagstat/flagstat")
req1 <- list(class = "StepInputExpressionRequirement")
req2 <- list(class = "InlineJavascriptRequirement")
mergeBamDup <- cwlStepParam(requirements = list(req1, req2),
                            inputs = InputParamList(p1, p2),
                            outputs = OutputParamList(o1, o2, o3, o4))
s1 <- Step(id = "mergeBam", run = mergeBam,
           In = list(ibam = "ibam",
                     obam = "obam"))
s2 <- Step(id = "markdup", run = markdup,
           In = list(ibam = "mergeBam/oBam",
                     obam = "obam",
                     matrix = list(valueFrom="$(inputs.ibam.nameroot).markdup.txt")))
s3 <- Step(id = "samtools_index", run = samtools_index,
           In = list(bam = "markdup/mBam"))
s4 <- Step(id = "samtools_flagstat", run = samtools_flagstat,
           In = list(bam = "markdup/mBam"))
mergeBamDup <- mergeBamDup + s1 + s2 + s3 + s4
