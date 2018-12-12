cwlVersion: v1.0
class: CommandLineTool
doc: |
  To convert a VCF into a MAF, each variant must be mapped to only one of all 
  possible gene transcripts/isoforms that it might affect. This selection of a 
  single effect per variant, is often subjective. So this project is an attempt 
  to make the selection criteria smarter, reproducible, and more configurable.
  
  This script needs VEP, a variant annotator that maps effects of a variant on 
  all possible genes and transcripts. For more info, see the README.

requirements:
  - class: DockerRequirement
    dockerPull: "opengenomics/vcf2maf"
  - class: InlineJavascriptRequirement

baseCommand: 
  - "perl"
  - "/opt/vcf2maf.pl"

arguments:
  - "--tmp-dir"
  - "."
  - "--vep-path"
  - "/home/vep"
  - "--filter-vcf"
  - $(inputs.vepData.path)/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz

inputs:
  inputVCF:
    type: File
    doc: "Path to input file in VCF format"
    inputBinding:
      prefix: "--input-vcf"

  outputMAF:
    type: string
    default: vep.maf
    doc: "Path to output MAF file [Default: STDOUT]"
    inputBinding:
      prefix: "--output-maf"

  tumorID:
    type: string?
    doc: "Tumor_Sample_Barcode to report in the MAF"
    inputBinding:
      prefix: "--tumor-id"

  normalID:
    type: string?
    doc: "Matched_Norm_Sample_Barcode to report in the MAF"
    inputBinding:
      prefix: "--normal-id"

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

  vepData:
    type: Directory
    doc: "VEP's base cache/plugin directory"
    inputBinding:
      prefix: "--vep-data"

  vepForks:
    type: int?
    doc: "Number of forked processes to use when running VEP [4]"
    inputBinding:
      prefix: "--vep-forks"

  refFasta:
    type: File
    doc: "Path to reference Fasta file"
    inputBinding:
        prefix: "--ref-fasta"
    secondaryFiles:
      - .fai
      - .gzi
  
  species:
    type: string?
    doc: "Ensembl-friendly name of species (e.g. mus_musculus for mouse) [homo_sapiens]"
    inputBinding:
      prefix: "--species"

  ncbiBuild:
    type: string?
    doc: "NCBI reference assembly of variants MAF (e.g. GRCm38 for mouse) [GRCh37]"
    inputBinding:
      prefix: "--ncbi-build"

  cacheVersion:
    type: int?
    doc: "Version of offline cache to use with VEP (e.g. 75, 82, 86) [Default: Installed version]"
    inputBinding:
      prefix: "--cache-version"

  mafCenter:
    type: string?
    doc: "Variant calling center to report in MAF [.]"
    inputBinding:
      prefix: "--maf-center"

  retainInfo:
    type: string[]?
    doc: "Comma-delimited names of INFO fields to retain as extra columns in MAF []"
    inputBinding:
      prefix: "--retain-info"
      itemSeparator: ","

  customEnst:
    type: File?
    doc: "List of custom ENST IDs that override canonical selection"
    inputBinding:
      prefix: "--custom-enst"

  minHomVaf:
    type: float?
    doc: "If GT undefined in VCF, minimum allele fraction to call a variant homozygous [0.7]"
    inputBinding:
      prefix: "--min-home-vaf"

  remapChain:
    type: File?
    doc: "Chain file to remap variants to a different assembly before running VEP"
    inputBinding:
      prefix: "--remap-chain"

  bufferSize:
    type: int?
    default: 4000
    inputBinding: 
      prefix: "--buffer-size"

outputs:
  maf:
    type: File
    outputBinding:
      glob: $(inputs.outputMAF)
  vepvcf:
    type: File
    outputBinding:
      glob: $(inputs.inputVCF.nameroot).vep.vcf
