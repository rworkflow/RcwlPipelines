cwlVersion: v1.0
class: CommandLineTool
doc: |
  The VCF files that variant callers generate are rarely compliant with VCF 
  specifications. This script fixes the most serious grievances, and creates a 
  VCF with only important fields in INFO and FORMAT, like GT:DP:AD. You may
  optionally specify --add-filters, to use these allele depths and fractions
  to add more tags under FILTER.

requirements:
  - class: DockerRequirement
    dockerPull: "opengenomics/vcf2maf"
        
baseCommand: 
  - "perl"
  - "/opt/vcf2vcf.pl"

inputs:
  inputVCF:
    type: File
    doc: "Path to input file in VCF format"
    inputBinding:
      prefix: "--input-vcf"

  outputVCF:
    type: string
    default: vep.vcf
    doc: "Path to output VCF file [Default: STDOUT]"
    inputBinding:
      prefix: "--output-vcf"

  refFasta:
    type: File
    doc: "Path to reference Fasta file"
    inputBinding:
      prefix: "--ref-fasta"
    secondaryFiles:
      - .fai
      - .gzi
  
  newTumorID:
    type: string?
    doc: "Tumor sample ID to use in the new VCF [Default: --vcf-tumor-id]"
    inputBinding:
      prefix: "--new-tumor-id"

  newNormalID:
    type: string?
    doc: "Matched normal ID to use in the new VCF [Default: --vcf-normal-id]"
    inputBinding:
      prefix: "--new-normal-id"

  vcfTumorID:
    type: string?
    doc: "Tumor sample ID used in VCF's genotype column [Default: TUMOR]"
    inputBinding:
      prefix: "--vcf-tumor-id"

  vcfNormalID:
    type: string?
    doc: "Matched normal ID used in VCF's genotype column [Default: NORMAL]"
    inputBinding:
      prefix: "--vcf-normal-id"

  addFilters:
    type: boolean
    doc: "Use this to add some extra tags under FILTER [Default: 0]"
    default: false
    inputBinding:
      prefix: "--add-filters"

outputs:
  vcf:
    type: File
    outputBinding:
      glob: $(inputs.outputVCF)
