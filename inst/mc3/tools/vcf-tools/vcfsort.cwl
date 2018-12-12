cwlVersion: v1.0
class: CommandLineTool

label: Sort a VCF file

requirements:
  - class: DockerRequirement
    dockerPull: opengenomics/vcftools-tools:latest

baseCommand: [vcf-sort, "-c"]

inputs:
  vcf:
    type: File
    inputBinding:
      position: 1

  output_name:
    default: sorted.vcf
    type: string

stdout: $(inputs.output_name)

outputs:
  output_vcf:
    type: File
    outputBinding:
      glob: $(inputs.output_name)
