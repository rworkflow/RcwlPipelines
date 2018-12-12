class: CommandLineTool
label: MuSE
cwlVersion: v1.0
baseCommand: [/opt/bin/muse.py, -O, muse.vcf, -w, ./, --muse, MuSEv1.0rc]
requirements:
  - class: "DockerRequirement"
    dockerPull: "opengenomics/muse:latest"
inputs:
  tumor:
    type: File
    inputBinding:
      prefix: --tumor-bam
    secondaryFiles:
      - .bai
  normal:
    type: File
    inputBinding:
      prefix: --normal-bam
    secondaryFiles:
      - .bai
  reference:
    type: File
    inputBinding:
      prefix: -f
  known:
    type: File
    inputBinding:
      prefix: -D
  mode:
    type: {"type": "enum", "name": "Mode", "symbols": ["wgs", "wxs"]}
    inputBinding:
      prefix: --mode
  ncpus:
    type: int
    default: 8
    inputBinding:
      prefix: "--cpus"
outputs:
  mutations:
    type: File
    outputBinding:
      glob: muse.vcf
