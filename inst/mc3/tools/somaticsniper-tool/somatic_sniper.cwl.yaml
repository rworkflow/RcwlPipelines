class: CommandLineTool
label: SomaticSniper
cwlVersion: v1.0
baseCommand: [ python, /opt/SomaticSniper.py, -F, vcf, --workdir, ./ ]
requirements:
  - class: "DockerRequirement"
    dockerPull: "opengenomics/somatic-sniper:latest"
inputs:
  reference:
    type: File
    inputBinding:
      prefix: -f
  tumor_name:
    type: string
    default: "TUMOR"
    inputBinding:
      prefix: -t
  normal_name:
    type: string
    default: "NORMAL"
    inputBinding:
      prefix: -n
  mapq:
    type: int
    inputBinding:
      prefix: -q
    default: 0
  somaticq:
    type: int
    inputBinding:
      prefix: -Q
    default: 40
  loh:
    type: boolean
    default: false
    inputBinding:
      prefix: -L
  gor:
    type: boolean
    default: false
    inputBinding:
      prefix: -G
  dis_priors:
    type: boolean
    default: false
    inputBinding:
      prefix: -p
  use_priorp:
    type: boolean
    default: false
    inputBinding:
      prefix: -J
  prior_p:
    type: float
    default: 0.01
    inputBinding:
      prefix: -s
  tumor:
    type: File
    inputBinding:
      position: 1
  normal:
    type: File
    inputBinding:
      position: 2
  output_name:
    type: string
    default: "somatic_sniper.vcf"
    inputBinding:
      position: 3

outputs:
  mutations:
    type: File
    outputBinding:
      glob: somatic_sniper.vcf
