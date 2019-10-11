## https://pachterlab.github.io/kallisto/

p1 <- InputParam(id = "index", type = "string", prefix = "-i")
p2 <- InputParam(id = "fasta", type = "File")
o1 <- OutputParam(id = "fidx", type = "File", glob = "$(inputs.index)")

req1 <- list(class = "DockerRequirement",
             dockerPull = "zlskidmore/kallisto")
kallisto_index <- cwlParam(baseCommand = c("kallisto", "index"),
                           requirements = list(req1),
                           inputs = InputParamList(p1, p2),
                           outputs = OutputParamList(o1))
