#' DNASeq alignment, merge and markduplicates
#'
#' The DNASeq pipeline to run bwa alignment, merge and mark
#' duplicates.
#' 
#' @format A `cwlStepParam` object.
#' \describe{
#'  \item{bwaAlign}{to align fastqs with bwa and sort with samtools}
#'  \item{mergeBamDup}{to merge BAMs from different flowcells
#' and then mark duplicates with picard}
#' }
#' @source \url{https://hubentu.github.io/others/Rcwl_DNASeq_Align.html}
#' @export
"alignMerge"

#' DNASeq alignment, merge, markduplicates and recalibration
#'
#' The DNASeq pipeline to run bwa alignment, merge, mark
#' duplicates and recalibration.
#' 
#' @format A `cwlStepParam` object.
#' \describe{
#'  \item{bwaAlign}{to align fastqs with bwa and sort with samtools}
#'  \item{mergeBamDup}{to merge BAMs from different flowcells
#' and then mark duplicates with picard}
#'  \item{BaseRecal}{Base quality recalibration}
#' }
#' @source \url{https://hubentu.github.io/others/Rcwl_DNASeq_Align.html}
#' @export
"bwaMMRecal"

#' DNASeq alignment, markduplicates and recalibration
#'
#' The DNASeq pipeline to run bwa alignment, mark duplicates and
#' recalibration.
#' 
#' @format A `cwlStepParam` object.
#' \describe{
#'  \item{bwaAlign}{to align fastqs with bwa and sort with samtools}
#'  \item{markdup}{to mark duplicates with picard}
#'  \item{BaseRecal}{Base quality recalibration}
#' }
#' @source \url{https://hubentu.github.io/others/Rcwl_DNASeq_Align.html}
#' @export
"bwaMRecal"


#' GATK alignment pipeline
#'
#' Workflows for processing high-throughput sequencing data for
#' variant discovery with GATK4 and related tools. Two workflows from
#' github, seq-format-conversion (last update: 7/13/2018) and
#' gatk4-data-processing (last update: 8/1/2018) were cloned to the
#' package.
#' 
#' @format A `cwlStepParam` object.
#' \describe{
#'  \item{fq2ubam}{To covert fastq to ubam with read group information}
#'  \item{align}{To run BWA alignment and BAM BaseRecalibration.}
#' }
#' @source \url{https://github.com/gatk-workflows/seq-format-conversion}
#' @source \url{https://github.com/gatk-workflows/gatk4-data-processing}
#' @source \url{https://hubentu.github.io/others/Rcwl_GATK4.html}
#' @export
"GAlign"

#' GATK haplotypecaller pipeline
#'
#' The workflow runs HaplotypeCaller from GATK4 in GVCF mode on a
#' single sample according to the GATK Best Practices (June 2016),
#' scattered across intervals. The workflow from github,
#' gatk4-germline-snps-indels (last update: 7/23/2018) was cloned to
#' this package.
#' 
#' @format A `cwlStepParam` object.
#' \describe{
#'  \item{HC}{HaplotypeCaller from GATK4}
#' }
#' @source \url{https://github.com/gatk-workflows/gatk4-germline-snps-indels}
#' @source \url{https://hubentu.github.io/others/Rcwl_GATK4.html}
#' @export
"hapCall"

#' GATK joint discovery pipeline
#'
#' The joint discovery and VQSR filtering portion of the GATK Best
#' Practices (June 2016) for germline SNP and Indel discovery in human
#' whole-genome sequencing (WGS) and exome sequencing data.
#'
#' @format A `cwlStepParam` object.
#' \describe{
#'  \item{JD}{variant joint genotyping}
#' }
#' @source \url{https://github.com/gatk-workflows/gatk4-germline-snps-indels}
#' @source \url{https://hubentu.github.io/others/Rcwl_GATK4.html}
#' @export
"jdCall"

#' GATK4: create a panel of normals
#'
#' The Panel of Normals Workflow
#'
#' @format A `cwlStepParam` object.
#' \describe{
#'  \item{GPoN}{The best practice pipeline to create a panel of normals.}
#' }
#' @source \url{https://software.broadinstitute.org/gatk/documentation/article?id=24057}
#' @export
"GPoN"

#' GATK4: Mutect2
#'
#' Somatic short variant discovery (SNVs + Indels)
#'
#' @format A `cwlStepParam` object.
#' \describe{
#'  \item{Mutect2PL}{The best practice pipeline to Identify somatic short
#' variants (SNVs and Indels).}
#' }
#' @source \url{https://software.broadinstitute.org/gatk/best-practices/workflow?id=11146}
#' @export
"Mutect2PL"

#' RNASeq pipeline with STAR and featureCounts
#'
#' RNASeq pipeline by STAR and featureCounts.
#' @format A `cwlStepParam` object.
#' \describe{
#'  An RNASeq alignment and quantification pipeline built by
#' `Rcwl`, which contains steps:
#'  \item{fastqc}{The reads QC step by fastQC}
#'  \item{STAR}{The alignment step by STAR}
#'  \item{samtools_index}{Index bam file by samtools}
#'  \item{samtools_flagstat}{Flag stat by samtools}
#'  \item{featureCounts}{Gene level quantification by featureCounts}
#'  \item{RSeQC}{QC for RNASeq alignments by RSeQC}
#' }
#' @source \url{https://hubentu.github.io/others/Rcwl_RNASeq.html}
#' @export
"rnaseq_Sf"

## #' TCGA MC3 pipeline
## #' #'
## #' #' The MC3 somatic variant calling pipeline. The mc3 workflow (last
## #' #' update: 11/5/2018) from github was cloned to the package.
## #' #'
## #' #' @format A `cwlStepParam` object.
## #' #' \describe{
## #' #'  \item{call_variants}{To call somatic variants by multiple pipelines}
## #' #'  \item{convert}{To convert and merge VCFs to MAF}
## #' #' }
## #' #' @source \url{https://github.com/OpenGenomics/mc3}
## #' #' @source \url{https://hubentu.github.io/others/Rcwl_MC3.html}
## #' @export
## "mc3"

#' RNASeq quality control by RSeQC
#'
#' RNASeq pipeline by STAR and featureCounts.
#' @format A `cwlStepParam` object.
#' \describe{
#'  An RNASeq QC pipeline by RSeQC
#'  which contains steps:
#'  \item{gtfToGenePred}{GTF to GenePred format}
#'  \item{genePredToBed}{GenePred format to Bed format}
#'  \item{read_distribution}{Reads distribution over genome feature}
#'  \item{geneBody_coverage}{Reads coverage over gene body}
#' }
#' @source \url{http://rseqc.sourceforge.net/}
#' @export
"RSeQC"

#' VarScan2 somatic caller
#'
#' VarScan2 Somatic caller pipeline.
#' @format A `cwlStepParam` object.
#' \describe{
#'  VarScan2 Somatic caller pipeline,
#'  which contains steps:
#'  \item{mpileup}{mpileup by samtools}
#'  \item{somatic}{somatic calling by VarScan2 somatic}
#'  \item{processSomatic}{processSomatic by VarScan2}
#'  \item{somaticFilter}{Filter by VarScan2}
#' }
#' @source \url{http://varscan.sourceforge.net}
#' @export
"VarScan2Somatic"
