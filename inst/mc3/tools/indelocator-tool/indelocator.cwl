cwlVersion: v1.0
class: CommandLineTool
label: Indelocator
baseCommand: ["java", "-Xmx7g", "-jar", "/opt/IndelGenotyper.jar", "-T", "IndelGenotyperV2"]
requirements:
  - class: DockerRequirement
    dockerPull: opengenomics/indelocator-tool:latest

inputs:
  somatic_flag:
    type: boolean
    default: true
    inputBinding:
      prefix: --somatic
  quiet_flag:
    type: boolean
    default: true
    inputBinding:
      prefix: -quiet
  window_size:
    type: int
    default: 300
    inputBinding:
      prefix: -ws
  normal:
    type: File
    inputBinding:
      position: 1
      prefix: --input_file:normal
    secondaryFiles:
      - .bai
  tumor:
    type: File
    inputBinding:
      position: 2
      prefix: --input_file:tumor
    secondaryFiles:
      - .bai
  reference:
    type: File
    inputBinding:
      prefix: -R
    secondaryFiles:
      - .fai
      - ^.dict
  bed_file:
    type: File?
    inputBinding:
      prefix: -L
  min_coverage:
    type: int
    default: 3
    inputBinding:
      prefix: --minCoverage
  vcf:
    type: string
    default: indelocator.vcf
    inputBinding:
      prefix: -o

outputs:
  mutations:
    type: File
    outputBinding:
      glob: $(inputs.vcf)
