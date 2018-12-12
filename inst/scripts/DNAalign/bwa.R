## bwa mem
p1 <- InputParam(id = "threads", type = "int", prefix = "-t", position = 1)
p2 <- InputParam(id = "RG", type = "string", prefix = "-R", position = 2)
p3 <- InputParam(id = "Ref", type = "string", position = 3)
p4 <- InputParam(id = "FQ1", type = "File", position = 4)
p5 <- InputParam(id = "FQ2", type = "File?", position = 5)
o1 <- OutputParam(id = "sam", type = "File", glob = "*.sam")
bwa <- cwlParam(baseCommand = c("bwa", "mem"),
                inputs = InputParamList(p1, p2, p3, p4, p5),
                output = OutputParamList(o1),
                stdout = "bwaOutput.sam")

## samtools sam to bam
p1 <- InputParam(id = "sam", type = "File")
o1 <- OutputParam(id = "bam", type = "File", glob = "$(inputs.sam.basename).bam")
sam2bam <- cwlParam(baseCommand = c("samtools", "view"),
                    arguments = list("-b"),
                    inputs = InputParamList(p1),
                    outputs = OutputParamList(o1),
                    stdout = "$(inputs.sam.basename).bam")

## samtools sort bam
p1 <- InputParam(id = "bam", type = "File")
o1 <- OutputParam(id = "sbam", type = "File", glob = "$(inputs.bam.nameroot).sorted.bam")
sortBam <- cwlParam(baseCommand = c("samtools", "sort"),
                    inputs = InputParamList(p1),
                    outputs = OutputParamList(o1),
                    stdout = "$(inputs.bam.nameroot).sorted.bam")

## Index bam
p1 <- InputParam(id = "bam", type = "File", position = 1)
o1 <- OutputParam(id = "idx", type = "File", glob = "$(inputs.bam.basename)", secondaryFiles = ".bai")
req1 <- list(class = "InitialWorkDirRequirement",
             listing = list("$(inputs.bam)"))
samtools_index <- cwlParam(baseCommand = c("samtools", "index"),
                           requirements = list(req1),
                           inputs = InputParamList(p1),
                           outputs = OutputParamList(o1))

## params
p1 <- InputParam(id = "threads", type = "int")
p2 <- InputParam(id = "RG", type = "string")
p3 <- InputParam(id = "Ref", type = "string")
p4 <- InputParam(id = "FQ1", type = "File")
p5 <- InputParam(id = "FQ2", type = "File")

o1 <- OutputParam(id = "Bam", type = "File", outputSource = "sortBam/sbam")
o2 <- OutputParam(id = "Idx", type = "File", outputSource = "idxBam/idx")

## stepParam
bwaAlign <- cwlStepParam(inputs = InputParamList(p1, p2, p3, p4, p5),
                         outputs = OutputParamList(o1, o2))
## bwa
s1 <- Step(id = "bwa", run = bwa,
           In = list(threads = "threads",
                     RG = "RG",
                     Ref = "Ref",
                     FQ1 = "FQ1",
                     FQ2 = "FQ2"))

## sam to bam
s2 <- Step(id = "sam2bam", run = sam2bam,
           In = list(sam = "bwa/sam"))

## sort bam
s3 <- Step(id = "sortBam", run = sortBam,
           In = list(bam = "sam2bam/bam"))
## index bam
s4 <- Step(id = "idxBam", run = samtools_index,
           In = list(bam = "sortBam/sbam"))

## pipeline
bwaAlign <- bwaAlign + s1 + s2 + s3 + s4
