cwlVersion: v1.0
class: CommandLineTool
label: Reformat Indelocator VCF
baseCommand: ["bash", "/opt/reformat_vcf.sh"]
requirements:
  - class: DockerRequirement
    dockerPull: opengenomics/indelocator-tool:latest

inputs:
  input_vcf:
    type: File
    inputBinding:
      prefix: -i
  output_vcf:
    type: string?
    inputBinding:
      prefix: -o

outputs:
  mutations:
    type: File
    outputBinding:
      glob: $(inputs.output_vcf)
