cwlVersion: v1.0
class: CommandLineTool
doc: |
  This script runs a given MAF through maf2vcf to generate per-TN-pair VCFs in 
  a temporary folder, and then runs vcf2maf on each VCF to reannotate variant 
  effects and create a new combined MAF

requirements:
  DockerRequirement:
    dockerPull: "opengenomics/vcf2maf"
        
baseCommand: 
  - "perl"
  - "/opt/maf2maf.pl"

arguments:
  - "--tmp-dir"
  - "."
  - "--vep-path"
  - "/home/vep"

inputs:
  inputMAF:
    type: File
    doc: "Path to input file in MAF format"
    inputBinding:
      prefix: "--input-maf"

  outputMAF:
    type: string?
    doc: "Path to output MAF file [Default: STDOUT]"
    inputBinding:
      prefix: "--output-maf"

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

  retainCols:
    type: string[]?
    doc: | 
      Comma-delimited list of columns to retain from the input MAF
      [Center,Verification_Status,Validation_Status,Mutation_Status,
       Sequencing_Phase,Sequence_Source,Validation_Method,Score,BAM_file,
       Sequencer,Tumor_Sample_UUID,Matched_Norm_Sample_UUID]
    inputBinding:
      prefix: "--retain-cols"
      itemSeparator: ","
      separate: false

  customEnst:
    type: File?
    doc: "List of custom ENST IDs that override canonical selection"
    inputBinding:
      prefix: "--custom-enst"

  vepData:
    type: Directory
    doc: "VEP's base cache/plugin directory"
    inputBinding:
      prefix: "--vep-data"
      valueFrom: $(inputs.vepData.dirname)

  vepForks:
    type: int?
    doc: "Number of forked processes to use when running VEP [4]"
    inputBinding:
      prefix: "--vep-forks"

  refFasta:
    type: File
    secondaryFiles:
      - .fai
      - .gzi
    doc: "Path to reference Fasta file"
    inputBinding:
        prefix: "--ref-fasta"
  
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
  tmpvcf:
    type: File
    outputBinding:
      glob: $(inputs.inputMAF.nameroot + ".vep.vcf")
