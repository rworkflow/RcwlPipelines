#!/usr/bin/env cwl-runner
#
# Author: Jeltje van Baren jeltje.van.baren@gmail.com

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [/opt/radia_filter.py]
label: RADIA-Filter

doc: "Filter radia output. If the input is from exomes plus RNA-Seq, set dnaOnly to False"

hints:
  DockerRequirement:
    dockerPull: opengenomics/radia:latest

requirements:
  - class: InlineJavascriptRequirement

inputs:

  inputVCF:
    type: File
    doc: |
      The input Radia vcf file
    inputBinding:
      position: 2
      prefix: --inputVCF

  patientId:
    type: string
    doc: |
      a unique patient Id that will be used to name the output file
    inputBinding:
      position: 2
      prefix: --patientId

  reference:
    type: File
    doc: |
      genome fasta file that was used to create the input bam files
    inputBinding:
      position: 2
      prefix: --fastaFilename

  normal:
    type: File
    doc: |
      the name of the normal DNA .bam file
    inputBinding:
      position: 2
      prefix: --dnaNormalFilename
    secondaryFiles:
      - .bai

  tumor:
    type: File
    doc: |
      the name of the tumor DNA .bam file
    inputBinding:
      position: 2
      prefix: --dnaTumorFilename
    secondaryFiles:
      - .bai

  out_vcf:
    type: string
    default: radia_filtered.vcf
    doc: |
      the name of the output file
    inputBinding:
      position: 2
      prefix: --outputFilename

  outputDir:
    type: string?
    doc: |
      the directory where temporary and final filtered output should be stored (default ./)
    inputBinding:
      position: 2
      prefix: --outputDir

  makeTCGAcompliant:
    type: boolean
    default: true
    doc: |
      Change VCF to make TCGA v1.1 compliant. This changes a lot in the VCF header and corrects the body (True)
    inputBinding:
      position: 2
      prefix: --makeTCGAcompliant

  filter-rejects:
    type: boolean
    default: true
    doc: |
      Filter out rejected calls (True)
    inputBinding:
      position: 2
      prefix: --filter-rejects

  filter-germline:
    type: boolean
    default: true
    doc: |
      Filter out germline calls (True)
    inputBinding:
      position: 2
      prefix: --filter-germline


  dnaNormalFastaFilename:
    type: File?
    doc: |
      the name of the fasta file that was used to create the BAM alignments
    inputBinding:
      position: 2
      prefix: --dnaNormalFastaFilename

  dnaTumorFastaFilename:
    type: File?
    doc: |
      the name of the fasta file that was used to create the BAM alignments
    inputBinding:
      position: 2
      prefix: --dnaTumorFastaFilename

  rnaNormalFilename:
    type: File?
    doc: |
      the name of the normal RNA .bam file
    inputBinding:
      position: 2
      prefix: --rnaNormalFilename

  rnaNormalBaiFilename:
    type: File?
    doc: |
      the name of the normal RNA .bai file
    inputBinding:
      position: 2
      prefix: --rnaNormalBaiFilename

  rnaNormalFastaFilename:
    type: File?
    doc: |
      the name of the fasta file that was used to create the BAM alignments
    inputBinding:
      position: 2
      prefix: --rnaNormalFastaFilename

  rnaTumorFilename:
    type: File?
    doc: |
      the name of the tumor RNA .bam file
    inputBinding:
      position: 2
      prefix: --rnaTumorFilename

  rnaTumorBaiFilename:
    type: File?
    doc: |
      the name of the tumor RNA .bai file
    inputBinding:
      position: 2
      prefix: --rnaTumorBaiFilename

  rnaTumorFastaFilename:
    type: File?
    doc: |
      the name of the fasta file that was used to create the BAM alignments
    inputBinding:
      position: 2
      prefix: --rnaTumorFastaFilename

  blacklistFilename:
    type: File?
    doc: |
      the name of the blacklist bed file
    inputBinding:
      position: 2
      prefix: --blacklistFilename

  targetFilename:
    type: File?
    doc: |
      the name of the exon capture targets file
    inputBinding:
      position: 2
      prefix: --targetFilename

  snpFilename:
    type: File?
    doc: |
      dbSNP vcf file
    inputBinding:
      position: 2
      prefix: --snpFilename

  retroGenesFilename:
    type: File?
    doc: |
      the name of the retrogenes bed file
    inputBinding:
      position: 2
      prefix: --retroGenesFilename

  pseudoGenesFilename:
    type: File?
    doc: |
      the name of the pseudogenes bed file
    inputBinding:
      position: 2
      prefix: --pseudoGenesFilename

  cosmicFilename:
    type: File?
    doc: |
      the name of the Catalogue Of Somatic Mutations In Cancer (COSMIC) annotations file
    inputBinding:
      position: 2
      prefix: --cosmicFilename

  canonical:
    type: boolean?
    doc: |
      include this argument if only the canonical transcripts from snpEff should be used (False)
    inputBinding:
      position: 2
      prefix: --canonical

  rnaGeneBlckFile:
    type: File?
    doc: |
      the RNA gene blacklist file
    inputBinding:
      position: 2
      prefix: --rnaGeneBlckFile

  rnaGeneFamilyBlckFile:
    type: File?
    doc: |
      the RNA gene family blacklist file
    inputBinding:
      position: 2
      prefix: --rnaGeneFamilyBlckFile

  blatFastaFilename:
    type: File?
    doc: |
      the fasta file that can be used during the BLAT filtering
    inputBinding:
      position: 2
      prefix: --blatFastaFilename

  noPositionalBias:
    type: boolean?
    doc: |
      include this argument if the positional bias filter should not be applied (False)
    inputBinding:
      position: 2
      prefix: --noPositionalBias

  dnaOnly:
    type: boolean
    default: True
    doc: |
      this argument if you only have DNA or filtering should only be done on the DNA (True)
    inputBinding:
      position: 2
      prefix: --dnaOnly

  number_of_procs:
    type: int?
    doc: |
      number of cpus to use (1)
    inputBinding:
      position: 2
      prefix: --number_of_procs

  workdir:
    type: string?
    doc: |
      work dir (./)
    inputBinding:
      position: 2
      prefix: --workdir

  no_clean:
    type: boolean?
    doc: |
      do not remove intermediate files (False)
    inputBinding:
      position: 2
      prefix: --no_clean


outputs:

  mutations:
    type: File
    outputBinding:
       glob: $(inputs.out_vcf)
