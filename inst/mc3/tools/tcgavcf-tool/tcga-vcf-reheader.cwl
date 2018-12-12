cwlVersion: v1.0
class: CommandLineTool
label: tcga-vcf-reheader
baseCommand: ["bash", "/opt/reheader_wrapper.sh"]
requirements:
  - class: DockerRequirement
    dockerPull: opengenomics/tcgavcf-tool:latest
  - class: InitialWorkDirRequirement
    listing: [ $(inputs.varscani_vcf), $(inputs.varscans_vcf), $(inputs.muse_vcf), $(inputs.mutect_vcf), $(inputs.somsniper_vcf), $(inputs.radia_vcf), $(inputs.pindel_vcf), $(inputs.indelocator_vcf) ]
  - class: InlineJavascriptRequirement
arguments:
  - shellQuote: false
inputs:
  varscani_vcf:
    type: File
    inputBinding:
      prefix: -i
  varscans_vcf:
    type: File
    inputBinding:
      prefix: -i
  muse_vcf:
    type: File
    inputBinding:
      prefix: -i
  mutect_vcf:
    type: File
    inputBinding:
      prefix: -i
  somsniper_vcf:
    type: File
    inputBinding:
      prefix: -i
  radia_vcf:
    type: File
    inputBinding:
      prefix: -i
  pindel_vcf:
    type: File
    inputBinding:
      prefix: -i
  indelocator_vcf:
    type: File
    inputBinding:
      prefix: -i
  tumor_analysis_uuid:
    type: string?
    inputBinding:
      prefix: -T
  tumor_bam_name:
    type: string
    inputBinding:
      prefix: -B
  tumor_aliquot_uuid:
    type: string?
    inputBinding:
      prefix: -X
  tumor_aliquot_name:
    type: string
    inputBinding:
      prefix: -A
  normal_analysis_uuid:
    type: string?
    inputBinding:
      prefix: -n
  normal_bam_name:
    type: string
    inputBinding:
      prefix: -b
  normal_aliquot_uuid:
    type: string?
    inputBinding:
      prefix: -x
  normal_aliquot_name:
    type: string
    inputBinding:
      prefix: -a
  platform:
    type: string?
    inputBinding:
      prefix: -p
  center:
    type: string?
    inputBinding:
      prefix: -c

outputs:
  reheaded_varscani:
    type: File
    outputBinding:
      glob: varscan_indel.reheadered.vcf
  reheaded_varscans:
    type: File
    outputBinding:
      glob: varscan_fpfilter.reheadered.vcf
  reheaded_muse:
    type: File
    outputBinding:
      glob: muse_filtered.reheadered.vcf
  reheaded_mutect:
    type: File
    outputBinding:
      glob: mutect.reheadered.vcf
  reheaded_radia:
    type: File
    outputBinding:
      glob: radia_filtered.reheadered.vcf
  reheaded_somsniper:
    type: File
    outputBinding:
      glob: somatic_sniper_fpfilter.reheadered.vcf
  reheaded_pindel:
    type: File
    outputBinding:
      glob: pindel_filtered.reheadered.vcf
  reheaded_indelocator:
    type: File
    outputBinding:
      glob: indelocator_filtered.reheadered.vcf

