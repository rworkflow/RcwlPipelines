### bwa_index

p1 <- InputParam(id = "Ref", type = "File")
o1 <- OutputParam(id = "idx", type = "File",
                  valueFrom = "$(self.dirname + '/' + self.basename)",
                  secondaryFiles = c("$(self.basename + '.amb')",
                                    "$(self.basename + '.ann')", 
                                    "$(self.basename + '.bwt')",
                                    "$(self.basename + '.pac')",
                                    "$(self.basename + '.sa')") 
req1 <- list(class = "DockerRequirement",
             dockerPull = "biocontainers/bwa:v0.7.15_cv3")
bwa <- cwlParam(baseCommand = c("bwa", "index"),
                requirements = list(req1),
                arguments = list("-a", "bwtsw")
                inputs = InputParamList(p1),
                output = OutputParamList(o1))
