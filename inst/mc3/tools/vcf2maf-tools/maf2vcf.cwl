cwlVersion: v1.0
class: CommandLineTool
doc: |
  This script breaks down variants in a MAF into a multi-sample VCF, in
  preparation for annotation with VEP. Can also create VCFs per-TN pair.

requirements:
  DockerRequirement:
    dockerPull: "opengenomics/vcf2maf"

baseCommand:
  - "perl"
  - "/opt/maf2vcf.pl"

arguments:
  - "--output-dir"
  - "./"

inputs:
  inputMaf:
    type: File
    doc: "Path to input file in MAF format"
    inputBinding:
      prefix: "--input-maf"

  outputVCF:
    type: string?
    doc: "Path to output multi-sample VCF containing all TN-pairs [<output-dir>/<input-maf-name>.vcf]"
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

  perTnVcfs:
    type: boolean?
    doc: "Specify this to generate VCFs per-TN pair, in addition to the multi-sample VCF"
    inputBinding:
        prefix: "--per-tn-vcfs"

  tumDepthCol:
    type: string?
    doc: "Name of MAF column for read depth in tumor BAM [t_depth]"
    inputBinding:
        prefix: "--tum-depth-col"

  tumRadCol:
    type: string?
    doc: "Name of MAF column for reference allele depth in tumor BAM [t_ref_count]"
    inputBinding:
        prefix: "--tum-rad-col"

  tumVadCol:
    type: string?
    doc: "Name of MAF column for variant allele depth in tumor BAM [t_alt_count]"
    inputBinding:
        prefix: "--tum-vad-col"

  nrmDepthCol:
    type: string?
    doc: "Name of MAF column for read depth in normal BAM [n_depth]"
    inputBinding:
        prefix: "--nrm-depth-col"

  nrmRadCol:
    type: string?
    doc: "Name of MAF column for reference allele depth in normal BAM [n_ref_count]"
    inputBinding:
        prefix: "--nrm-rad-col"

  nrmVadCol:
    type: string?
    doc: "Name of MAF column for variant allele depth in normal BAM [n_alt_count]"
    inputBinding:
        prefix: "--nrm-vad-col"

outputs:
  vcf:
    type: File[]
    outputBinding:
      glob: "*.vcf"
  skipped:
    type: File
    outputBinding:
      glob: "*.skipped.tsv"
