#!/usr/bin/env cwl-runner
#
# Author: Jeltje van Baren jeltje.van.baren@gmail.com

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [/opt/radia.py]
label: RADIA
doc: "Runs radia on individual chromosomes, then merges output. An input DNA (exome) sample pair is required, additional RNA-Seq data for the same sample is advised but optional. Radia performes significantly better when RNA data is added."

hints:
  DockerRequirement:
    dockerPull: opengenomics/radia:latest

requirements:
  - class: InlineJavascriptRequirement

inputs:
  reference:
    type: File
    doc: |
      genome fasta file that was used to create the input bam files
    inputBinding:
      position: 3
      prefix: --fastaFilename
  # Note: It is possible to specify fasta files for each individual input, see below

  normal:
    type: File
    doc: |
      the name of the normal DNA .bam file
    inputBinding:
      position: 3
      prefix: --dnaNormalFilename
    secondaryFiles:
      - .bai

  tumor:
    type: File
    doc: |
      the name of the tumor DNA .bam file
    inputBinding:
      position: 3
      prefix: --dnaTumorFilename
    secondaryFiles:
      - .bai

  out_vcf:
    type: string
    default: radia.vcf
    doc: |
      the name of the output vcf (radia.vcf)
    inputBinding:
      position: 3
      prefix: --outputFilename

  refId:
    type: string?
    doc: |
      the reference Id - used in the reference VCF meta tag
    inputBinding:
      position: 3
      prefix: --refId

  outputDir:
    type: string?
    doc: |
      the directory where temporary and final filtered output should be stored ('./')
    inputBinding:
      position: 3
      prefix: --outputDir

  patientId:
    type: string?
    default: myPatient
    doc: |
      a unique patient Id that will be added to the SAMPLE and INDIVIDUAL tags in the vcf header (myPatient)
    inputBinding:
      position: 3
      prefix: --patientId

  useChrPrefix:
    type: boolean?
    doc: |
      this argument if the 'chr' prefix should be used in the samtools command for all bams (False)
    inputBinding:
      position: 3
      prefix: --useChrPrefix

  log:
    type: string?
    doc: |
      the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL) (WARNING)
    inputBinding:
      position: 3
      prefix: --log

  refUrl:
    type: string?
    doc: |
      the URL for the reference - used in the reference VCF meta tag
    inputBinding:
      position: 3
      prefix: --refUrl

  refFilename:
    type: string?
    doc: |
      the location of the reference - used in the reference VCF meta tag
    inputBinding:
      position: 3
      prefix: --refFilename

  statsDir:
    type: string?
    doc: |
      a stats directory where some basic stats can be output
    inputBinding:
      position: 3
      prefix: --statsDir

  logFilename:
    type: string?
    doc: |
      the name of the log file, STDOUT by default
    inputBinding:
      position: 3
      prefix: --logFilename

  batchSize:
    type: int?
    doc: |
      the size of the samtools selections that are loaded into memory at one time (250000000)
    inputBinding:
      position: 3
      prefix: --batchSize

  dataSource:
    type: string?
    doc: |
      the source of the data - used in the sample VCF meta tag
    inputBinding:
      position: 3
      prefix: --dataSource

  sequencingPlatform:
    type: string?
    doc: |
      the sequencing platform - used in the sample VCF meta tag
    inputBinding:
      position: 3
      prefix: --sequencingPlatform

  disease:
    type: string?
    doc: |
      a disease abbreviation (i.e. BRCA) for the header
    inputBinding:
      position: 3
      prefix: --disease

  genotypeMinDepth:
    type: int?
    doc: |
      the minimum number of bases required for the genotype (2)
    inputBinding:
      position: 3
      prefix: --genotypeMinDepth

  genotypeMinPct:
    type: float?
    doc: |
      the minimum percentage of reads required for the genotype (.10)
    inputBinding:
      position: 3
      prefix: --genotypeMinPct

  gzip:
    type: boolean?
    doc: |
      this argument if the final VCF should be compressed with gzip (False)
    inputBinding:
      position: 3
      prefix: --gzip

  rnaNormalFilename:
    type: File?
    doc: |
      the name of the normal RNA-Seq .bam file
    inputBinding:
      position: 3
      prefix: --rnaNormalFilename

  rnaNormalBaiFilename:
    type: File?
    doc: |
      the name of the normal RNA-Seq .bai file
    inputBinding:
      position: 3
      prefix: --rnaNormalBaiFilename

  rnaTumorFilename:
    type: File?
    doc: |
      the name of the tumor RNA-Seq .bam file
    inputBinding:
      position: 3
      prefix: --rnaTumorFilename

  rnaTumorBaiFilename:
    type: File?
    doc: |
      the name of the tumor RNA-Seq .bai file
    inputBinding:
      position: 3
      prefix: --rnaTumorBaiFilename

  dnaNormalMinTotalBases:
    type: int?
    doc: |
      the minimum number of overall normal DNA reads covering a position (4)
    inputBinding:
      position: 3
      prefix: --dnaNormalMinTotalBases

  dnaNormalMinAltBases:
    type: int?
    doc: |
      the minimum number of alternative normal DNA reads supporting a variant at a position (2)
    inputBinding:
      position: 3
      prefix: --dnaNormalMinAltBases

  dnaNormalBaseQual:
    type: int?
    doc: |
      the minimum normal DNA base quality (10)
    inputBinding:
      position: 3
      prefix: --dnaNormalBaseQual

  dnaNormalMapQual:
    type: int?
    doc: |
      the minimum normal DNA mapping quality (10)
    inputBinding:
      position: 3
      prefix: --dnaNormalMapQual

  dnaNormalFasta:
    type: File?
    doc: |
      the name of the genome fasta file for the normal DNA .bam file
    inputBinding:
      position: 3
      prefix: --dnaNormalFasta

  dnaNormalDescription:
    type: string?
    doc: |
      the description for the sample in the VCF header (NormalDNASample)
    inputBinding:
      position: 3
      prefix: --dnaNormalDescription

  rnaNormalMinTotalBases:
    type: int?
    doc: |
      the minimum number of overall normal RNA-Seq reads covering a position (4)
    inputBinding:
      position: 3
      prefix: --rnaNormalMinTotalBases

  rnaNormalMinAltBases:
    type: int?
    doc: |
      the minimum number of alternative normal RNA-Seq reads supporting a variant at a position (2)
    inputBinding:
      position: 3
      prefix: --rnaNormalMinAltBases

  rnaNormalBaseQual:
    type: int?
    doc: |
      the minimum normal RNA-Seq base quality (10)
    inputBinding:
      position: 3
      prefix: --rnaNormalBaseQual

  rnaNormalMapQual:
    type: int?
    doc: |
      the minimum normal RNA-Seq mapping quality (10)
    inputBinding:
      position: 3
      prefix: --rnaNormalMapQual

  rnaNormalFasta:
    type: File?
    doc: |
      the name of the fasta file for the normal RNA .bam file
    inputBinding:
      position: 3
      prefix: --rnaNormalFasta

  rnaNormalDescription:
    type: string?
    doc: |
      the description for the sample in the VCF header (NormalRNASample)
    inputBinding:
      position: 3
      prefix: --rnaNormalDescription

  dnaTumorMinTotalBases:
    type: int?
    doc: |
      the minimum number of overall tumor DNA reads covering a position (4)
    inputBinding:
      position: 3
      prefix: --dnaTumorMinTotalBases

  dnaTumorMinAltBases:
    type: int?
    doc: |
      the minimum number of alternative tumor DNA reads supporting a variant at a position (2)
    inputBinding:
      position: 3
      prefix: --dnaTumorMinAltBases

  dnaTumorBaseQual:
    type: int?
    doc: |
      the minimum tumor DNA base quality (10)
    inputBinding:
      position: 3
      prefix: --dnaTumorBaseQual

  dnaTumorMapQual:
    type: int?
    doc: |
      the minimum tumor DNA mapping quality (10)
    inputBinding:
      position: 3
      prefix: --dnaTumorMapQual

  dnaTumorFasta:
    type: File?
    doc: |
      the name of the fasta file for the tumor DNA .bam file
    inputBinding:
      position: 3
      prefix: --dnaTumorFasta

  dnaTumorDescription:
    type: string?
    doc: |
      the description for the sample in the VCF header (TumorDNASample)
    inputBinding:
      position: 3
      prefix: --dnaTumorDescription

  rnaTumorMinTotalBases:
    type: int?
    doc: |
      the minimum number of overall tumor RNA-Seq reads covering a position (4)
    inputBinding:
      position: 3
      prefix: --rnaTumorMinTotalBases

  rnaTumorMinAltBases:
    type: int?
    doc: |
      the minimum number of alternative tumor RNA-Seq reads supporting a variant at a position (2)
    inputBinding:
      position: 3
      prefix: --rnaTumorMinAltBases

  rnaTumorBaseQual:
    type: int?
    doc: |
      the minimum tumor RNA-Seq base quality (10)
    inputBinding:
      position: 3
      prefix: --rnaTumorBaseQual

  rnaTumorMapQual:
    type: int?
    doc: |
      the minimum tumor RNA-Seq mapping quality (10)
    inputBinding:
      position: 3
      prefix: --rnaTumorMapQual

  rnaTumorFasta:
    type: File?
    doc: |
      the name of the fasta file for the tumor RNA .bam file
    inputBinding:
      position: 3
      prefix: --rnaTumorFasta

  rnaTumorDescription:
    type: string?
    doc: |
      the description for the sample in the VCF header (TumorRNASample)
    inputBinding:
      position: 3
      prefix: --rnaTumorDescription

  number_of_procs:
    type: int?
    doc: |
      number of cpus to use (1)
    inputBinding:
      position: 3
      prefix: --number_of_procs

  workdir:
    type: string?
    doc: |
      working directory to use (./)
    inputBinding:
      position: 3
      prefix: --workdir

  no_clean:
    type: boolean?
    doc: |
      do not remove intermediate files (False)
    inputBinding:
      position: 3
      prefix: --no_clean


outputs:

  mutations:
    type: File
    outputBinding:
       glob: $(inputs.out_vcf)
