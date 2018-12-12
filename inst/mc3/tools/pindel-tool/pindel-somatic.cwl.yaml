
class: CommandLineTool
label: Pindel-Somatic
cwlVersion: v1.0
baseCommand: ["python", "/opt/pindel.py", "-t", "NORMAL", "-t", "TUMOR"]
requirements:
  - class: DockerRequirement
    dockerPull: opengenomics/pindel:latest

inputs:
  - id: normal
    type: File
    inputBinding:
      position: 1
      prefix: -b

  - id: tumor
    type: File
    inputBinding:
      position: 2
      prefix: -b

  - id: normal_insert_size
    type: int?
    inputBinding:
      position: 1
      prefix: -s

  - id: tumor_insert_size
    type: int?
    inputBinding:
      position: 2
      prefix: -s

  - id: reference
    type: File
    inputBinding:
      prefix: -r

  - id: centromere
    type: File
    inputBinding:
      prefix: -J

  - id: referenceName
    type: string
    default: HG19
    inputBinding:
      prefix: -R

  - id: window_size
    type: float
    default: 5.0
    inputBinding:
      prefix: --window_size

  - id: balance_cutoff
    type: int
    default: 100
    inputBinding:
      prefix: --balance_cutoff

  - id: procs
    type: int
    default: 2
    inputBinding:
      prefix: --number_of_procs

  - id: report_inversions
    type: boolean
    default: false
    inputBinding:
      prefix: --report_inversions

  - id: report_duplications
    type: boolean
    default: false
    inputBinding:
      prefix: --report_duplications

  - id: report_long_insertions
    type: boolean
    default: false
    inputBinding:
      prefix: --report_long_insertions

  - id: report_breakpoints
    type: boolean
    default: false
    inputBinding:
      prefix: --report_breakpoints

  - id: report_only_close_mapped_reads
    type: boolean
    default: false
    inputBinding:
      prefix: -S

  - id: outputRawFile
    type: string
    default: pindel.raw
    inputBinding:
      prefix: -o1
  
  - id: outputVcfFile
    type: string
    default: pindel.vcf
    inputBinding:
      prefix: -o2
  
  - id: outputSomaticVcfFile
    type: string
    default: pindel_somatic.vcf
    inputBinding:
      prefix: -o3
  
  - id: somatic_vaf
    type: float
    default: 0.08
    inputBinding:
      prefix: --somatic_vaf

  - id: somatic_cov
    type: int
    default: 20
    inputBinding:
      prefix: --somatic_cov

  - id: somatic_hom
    type: int
    default: 6
    inputBinding:
      prefix: --somatic_hom

  - id: min_inversion_size
    type: int
    default: 50
    inputBinding:
      prefix: --min_inversion_size

  - id: min_num_matched_bases
    type: int
    default: 30
    inputBinding:
      prefix: -d

  - id: max_range_index
    type: int
    default: 4
    inputBinding:
      prefix: -x

  - id: additional_mismatch
    type: int
    default: 1
    inputBinding:
      prefix: -a
  
  - id: min_perfect_match_around_BP
    type: int
    default: 3
    inputBinding:
      prefix: -m
  
  - id: sequencing_error_rate
    type: float
    default: 0.01
    inputBinding:
      prefix: --sequencing_error_rate

  - id: maximum_allowed_mismatch_rate
    type: float
    default: 0.02
    inputBinding:
      prefix: --maximum_allowed_mismatch_rate

  - id: minimum_support_for_event
    type: int
    default: 1
    inputBinding:
      prefix: -M
  
  - id: sensitivity
    type: float
    default: 0.95
    inputBinding:
      prefix: --sensitivity


outputs:
  - id: vcf
    type: File
    outputBinding:
      glob: pindel.vcf

  - id: somatic_vcf
    type: File
    outputBinding:
      glob: pindel_somatic.vcf
      
  - id: rawFile
    type: File
    outputBinding:
      glob: pindel.raw

        
