cwlVersion: v1.0
class: CommandLineTool

label: Merge VCFs

requirements:
  - class: DockerRequirement
    dockerPull: opengenomics/vcftools-tools:latest

baseCommand: [python, /opt/merge_vcfs.py]

#stdout: stdout
#stderr: stderr

inputs:

  output_name:
    type: string
    default: merged.vcf
    inputBinding:
      prefix: --output

  keys:
    type: string[]
    inputBinding:
      prefix: --keys

  vcfs:
    type: File[]
    inputBinding:
      prefix: --input

outputs:
  output_vcf:
    type: File
    outputBinding:
      glob: $(inputs.output_name)
