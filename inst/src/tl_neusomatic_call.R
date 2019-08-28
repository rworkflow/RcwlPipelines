
p1 <- InputParam(id = "candidates", type = "File[]",
                 prefix = "--candidates_tsv", secondaryFiles = ".idx")
p2 <- InputParam(id = "ref", type = "File",
                 prefix = "--reference", secondaryFiles = ".fai")
p3 <- InputParam(id = "checkpoint", type = "string", prefix = "--checkpoint",
                 default = "/opt/neusomatic/neusomatic/models/NeuSomatic_v0.1.4_ensemble_SEQC-WGS-Spike.pth")
p4 <- InputParam(id = "threads", type = "int", prefix = "--num_threads", default = 1L)
o1 <- OutputParam(id = "pred", type = "File", glob = "pred.vcf")

req1 <- list(class = "DockerRequirement",
             dockerPull = "msahraeian/neusomatic")
neusomatic_call <- cwlParam(baseCommand = c("python",
                                            "/opt/neusomatic/neusomatic/python/call.py"),
                            requirements = list(req1),
                            arguments = list("--out", ".",
                                             "--ensemble",
                                             "--batch_size", "100"),
                            inputs = InputParamList(p1, p2, p3),
                            outputs = OutputParamList(o1))
