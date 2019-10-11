## https://pachterlab.github.io/kallisto/

p1 <- InputParam(id = "index", type = "File", prefix = "-i")
p2 <- InputParam(id = "fastq", type = "File[]")
p3 <- InputParam(id = "threads", type = "int", prefix = "-t")
o1 <- OutputParam(id = "h5", type = "File", glob = "abundance.h5")
o2 <- OutputParam(id = "tsv", type = "File", glob = "abundance.tsv")
o3 <- OutputParam(id = "info", type = "File", glob = "run_info.json")

req1 <- list(class = "DockerRequirement",
             dockerPull = "zlskidmore/kallisto")
kallisto_quant <- cwlParam(baseCommand = c("kallisto", "quant"),
                           requirements = list(req1),
                           arguments = list("-o", "./"),
                           inputs = InputParamList(p1, p2, p3),
                           outputs = OutputParamList(o1, o2, o3))
