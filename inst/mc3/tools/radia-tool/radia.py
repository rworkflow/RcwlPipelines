#!/usr/bin/env python

import sys
import argparse
import os
import re
import shutil
import subprocess
import tempfile
import time
from multiprocessing import Pool


def execute(cList):
    import shlex
    """ function to execute a cmd and report if an error occurs. Takes in a list with one or two arguments, the second one optionally being the output file"""
    cmd = cList[0]
    print(cmd)
    try:
        output = cList[1]
    except:
        output = None
    try:
        process = subprocess.Popen(args=shlex.split(
            cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
    except Exception as e:  # error from my command : stderr
        sys.stderr.write("problem doing : %s\n%s\n" % (cmd, e))
        return 1
    if output:
        output = open(output, 'w')
        output.write(stdout)
        output.close()
    if stderr != '':  # internal program error : stdout
        sys.stdout.write(
            "warning or error while doing : %s\n-----\n%s-----\n\n" %
            (cmd, stderr))
        return 1
    return 0


def indexBam(workdir, prefix, inputBamFile, inputBamFileIndex=None):
    inputBamLink = os.path.join(os.path.abspath(workdir), prefix + ".bam")
    os.symlink(inputBamFile, inputBamLink)
    if os.path.exists(inputBamFile + ".bai"):
        inputBamFileIndex = inputBamFile + ".bai"
    if inputBamFileIndex is None or inputBamFileIndex == "None":
        cmd = "samtools index %s" % (inputBamLink)
        execute([cmd])
    else:
        os.symlink(inputBamFileIndex, inputBamLink + ".bai")
    return inputBamLink


def indexFasta(
        workdir,
        inputFastaFile,
        inputFastaFileIndex=None,
        prefix="dna"):
    """Checks if fasta index exists. If so, creates link. If not, creates index"""
    inputFastaLink = os.path.join(
        os.path.abspath(workdir), prefix + "_reference.fa")
    os.symlink(inputFastaFile, inputFastaLink)
    inputFastaFileIndex = inputFastaFile + ".fai"
    if os.path.exists(inputFastaFileIndex):
        os.symlink(inputFastaFileIndex, inputFastaLink + ".fai")
    else:
        cmd = "samtools faidx %s" % (inputFastaLink)
        execute([cmd])
    return inputFastaLink


def idxStats(bamfile):
    """runs samtools idxstats"""
    samtools = which("samtools")
    cmd = [samtools, "idxstats", bamfile]
    process = subprocess.Popen(args=cmd, stdout=subprocess.PIPE)
    stdout, stderr = process.communicate()
    return stdout


def mitName(idx):
    """Returns the mitochondrion chromosome ID used in the bam file if it starts with M (usually it is called M or MT)"""
    for line in idx.split("\n"):
        tmp = line.split("\t")
        if len(tmp) == 4 and tmp[0].startswith("chrM"):
            return tmp[0][3:]  # remove chr
        if len(tmp) == 4 and tmp[0].startswith("M"):
            return tmp[0]
    return 'M'  # not found, so does not matter


def bamChrScan(idx):
    """Checks if the bam chromosome IDs start with chr"""
    for line in idx.split("\n"):
        tmp = line.split("\t")
        if len(tmp) == 4 and tmp[0].startswith("chr"):
            return True
    return False


def addNumsAndQuals(args, cmd, sample):
    """Append mapping quality parameters to radia command"""
    if sample == 'dnaNormal':
        cmd += " --dnaNormalDescription %s --dnaNormalMinTotalBases %d --dnaNormalMinAltBases %d --dnaNormalBaseQual %d --dnaNormalMapQual %d" % (
            args.dnaNormalDesc, args.dnaNormalMinTotalBases, args.dnaNormalMinAltBases, args.dnaNormalMinBaseQual, args.dnaNormalMinMappingQual)
        return cmd

    if sample == 'dnaTumor':
        cmd += " --dnaTumorDescription %s --dnaTumorMinTotalBases %d --dnaTumorMinAltBases %d --dnaTumorBaseQual %d --dnaTumorMapQual %d" % (
            args.dnaTumorDesc, args.dnaTumorMinTotalBases, args.dnaTumorMinAltBases, args.dnaTumorMinBaseQual, args.dnaTumorMinMappingQual)
        return cmd

    if sample == 'rnaNormal':
        cmd += " --rnaNormalDescription %s --rnaNormalMinTotalBases %d --rnaNormalMinAltBases %d --rnaNormalBaseQual %d --rnaNormalMapQual %d" % (
            args.rnaNormalDesc, args.rnaNormalMinTotalBases, args.rnaNormalMinAltBases, args.rnaNormalMinBaseQual, args.rnaNormalMinMappingQual)
        return cmd

    if sample == 'rnaTumor':
        cmd += " --rnaTumorDescription %s --rnaTumorMinTotalBases %d --rnaTumorMinAltBases %d --rnaTumorBaseQual %d --rnaTumorMapQual %d" % (
            args.rnaTumorDesc, args.rnaTumorMinTotalBases, args.rnaTumorMinAltBases, args.rnaTumorMinBaseQual, args.rnaTumorMinMappingQual)
        return cmd


def radia(
        chrom,
        args,
        outputDir,
        dnaNormalFilename=None,
        rnaNormalFilename=None,
        dnaTumorFilename=None,
        rnaTumorFilename=None,
        dnaNormalFastaFilename=None,
        rnaNormalFastaFilename=None,
        dnaTumorFastaFilename=None,
        rnaTumorFastaFilename=None):

    # python radia.py id chrom [Options]

    # -i GRCh37
    # -m Homo_sapiens_assembly19.fasta
    # -d CGHub
    # -q Illumina
    # --disease GBM

    # quadruplets
    if (rnaNormalFilename is not None and rnaTumorFilename is not None):
        cmd = "python %s/radia.py %s %s -n %s -x %s  -t %s -r %s  --dnaNormalFasta %s --rnaNormalFasta %s --dnaTumorFasta %s --rnaTumorFasta %s " % (
            args.scriptsDir,
            args.patientId, chrom,
            dnaNormalFilename,
            rnaNormalFilename,
            dnaTumorFilename,
            rnaTumorFilename,
            dnaNormalFastaFilename, rnaNormalFastaFilename, dnaTumorFastaFilename, rnaTumorFastaFilename)
        cmd = addNumsAndQuals(args, cmd, "dnaNormal")
        cmd = addNumsAndQuals(args, cmd, "dnaTumor")
        cmd = addNumsAndQuals(args, cmd, "rnaTumor")
        cmd = addNumsAndQuals(args, cmd, "rnaNormal")
    # triplets
    elif (rnaTumorFilename is not None):
        cmd = "python %s/radia.py %s %s -n %s -t %s -r %s --dnaNormalFasta %s --dnaTumorFasta %s --rnaTumorFasta %s " % (
            args.scriptsDir,
            args.patientId, chrom,
            dnaNormalFilename,
            dnaTumorFilename,
            rnaTumorFilename,
            dnaNormalFastaFilename, dnaTumorFastaFilename, rnaTumorFastaFilename)
        cmd = addNumsAndQuals(args, cmd, "dnaNormal")
        cmd = addNumsAndQuals(args, cmd, "dnaTumor")
        cmd = addNumsAndQuals(args, cmd, "rnaTumor")
    # pairs
    else:
        cmd = "python %s/radia.py %s %s -n %s -t %s --dnaNormalFasta %s --dnaTumorFasta %s " % (
            args.scriptsDir,
            args.patientId, chrom,
            dnaNormalFilename,
            dnaTumorFilename,
            dnaNormalFastaFilename, dnaTumorFastaFilename)
        cmd = addNumsAndQuals(args, cmd, "dnaNormal")
        cmd = addNumsAndQuals(args, cmd, "dnaTumor")

    # determine naming for chromosomes (with or without 'chr') and
    # mitochondrion (M or something starting with M)
    if dnaNormalFilename is not None:
        idx = idxStats(dnaNormalFilename)
        cmd += ' --dnaNormalMitochon=' + mitName(idx)
    if rnaNormalFilename is not None:
        idx = idxStats(rnaNormalFilename)
        cmd += ' --rnaNormalMitochon=' + mitName(idx)
    if dnaTumorFilename is not None:
        idx = idxStats(dnaTumorFilename)
        cmd += ' --dnaTumorMitochon=' + mitName(idx)
    if rnaTumorFilename is not None:
        idx = idxStats(rnaTumorFilename)
        cmd += ' --rnaTumorMitochon=' + mitName(idx)

    # add genotype parameters
    cmd += ' --genotypeMinDepth %d --genotypeMinPct %.3f' % (
        args.genotypeMinDepth,
        args.genotypeMinPct)

    # vcf header arguments
    if args.refId:
        cmd += ' --refId %s' % args.refId
    if args.refUrl:
        cmd += ' --refUrl %s' % args.refUrl
    if args.refFilename:
        cmd += ' --refFilename %s' % args.refFilename
    if args.dataSource:
        cmd += ' --dataSource %s' % args.dataSource
    if args.sequencingPlatform:
        cmd += ' --sequencingPlatform %s' % args.sequencingPlatform
    if args.disease:
        cmd += ' --disease %s' % args.disease

    outfile = os.path.join(outputDir, args.patientId + "_chr" + chrom + ".vcf")
    if args.gzip:
        cmd += ' --gzip '
        outfile += '.gz'
    return cmd, outfile


def radiaMerge(args, inputDir):
    """Merges vcf files if they follow the pattern patientID_chr<N>.vcf(.gz)"""

    # python mergeChroms.py patientId /radia/filteredChroms/ /radia/filteredPatients/ --gzip
    #  -h, --help            show this help message and exit
    #  -o OUTPUT_FILE, --outputFilename=OUTPUT_FILE
    #                   the name of the output file, <id>.vcf(.gz) by default

    #  -l LOG, --log=LOG     the logging level (DEBUG, INFO, WARNING, ERROR,
    #                        CRITICAL), WARNING by default
    #  -g LOG_FILE, --logFilename=LOG_FILE
    #                        the name of the log file, STDOUT by default
    #  --gzip                include this argument if the final VCF should be
    #                        compressed with gzip
    # radia works in the workdir
    cmd = "python %s/mergeChroms.py %s %s %s -o %s" % (
        args.scriptsDir,
        args.patientId, inputDir, args.workdir,
        args.outputFilename)
    return cmd


def which(cmd):
    cmd = ["which", cmd]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    res = p.stdout.readline().rstrip()
    if len(res) == 0:
        return None
    return res


def identicalName(inputList):
    """returns duplicate name if two inputs have the same name and are not None"""
    dup = set(x for x in inputList if inputList.count(x) >= 2)
    dup.discard(None)  # this doesn't complain if None is not in the set
    if dup:
        print "ERROR: found duplicate input %s" % dup.pop()
        return True
    return False


def removeSpaces(mystring):
    if mystring:
        return ("_").join(mystring.split(" "))
    return False


def get_bam_seq(inputBamFile, exclude):
    samtools = which("samtools")
    cmd = [samtools, "idxstats", inputBamFile]
    print "calling", cmd
    process = subprocess.Popen(args=cmd, stdout=subprocess.PIPE)
    stdout, stderr = process.communicate()
    seqs = []
    for line in stdout.split("\n"):
        tmp = line.split("\t")
        if len(tmp) == 4 and tmp[2] != "0":
            seqs.append(tmp[0])
    pattern = "|".join(exclude)
    r = re.compile(pattern)
    return filter(lambda i: not r.match(i), seqs)


def __main__():
    # small hack, sometimes it seems like docker file systems are avalible
    # instantly
    time.sleep(1)
    parser = argparse.ArgumentParser(
        description="RNA and DNA Integrated Analysis (RADIA)")

    #############################
    #    RADIA params    #
    #############################
    parser.add_argument(
        "-o",
        "--outputFilename",
        dest="outputFilename",
        metavar="OUTPUT_FILE",
        default='out.vcf',
        help="the name of the output file")
    parser.add_argument(
        "--outputDir",
        dest="outputDir",
        default='./',
        metavar="FILTER_OUT_DIR",
        help="the directory where temporary and final filtered output should be stored")
    parser.add_argument(
        "--scriptsDir",
        dest="scriptsDir",
        default='/opt/radia-1.1.5/scripts',
        metavar="SCRIPTS_DIR",
        help="the directory that contains the RADIA filter scripts")
    parser.add_argument(
        "--patientId",
        dest="patientId",
        required=True,
        metavar="PATIENT_ID",
        help="a unique patient Id that will be used to name the output file")

    parser.add_argument(
        "-f",
        "--fastaFilename",
        dest="fastaFilename",
        metavar="FASTA_FILE",
        help="the name of the fasta file that can be used on all .bams, see below for specifying individual fasta files for each .bam file")
    parser.add_argument(
        "-i",
        "--refId",
        dest="refId",
        metavar="REF_ID",
        help="the reference Id - used in the reference VCF meta tag")
    parser.add_argument(
        "-u",
        "--refUrl",
        dest="refUrl",
        metavar="REF_URL",
        help="the URL for the reference - used in the reference VCF meta tag")
    parser.add_argument(
        "-m",
        "--refFilename",
        dest="refFilename",
        metavar="REF_FILE",
        help="the location of the reference - used in the reference VCF meta tag")
##### not implemented in this wrapper  ####
    parser.add_argument(
        "-a",
        "--startCoordinate",
        type=int,
        default=int(1),
        dest="startCoordinate",
        metavar="START_COORDINATE",
        help="the start coordinate for testing small regions, default 1")
    parser.add_argument(
        "-z",
        "--stopCoordinate",
        type=int,
        default=int(0),
        dest="stopCoordinate",
        metavar="STOP_COORDINATE",
        help="the stop coordinate for testing small regions, default 0")
    parser.add_argument(
        "-s",
        "--statsDir",
        dest="statsDir",
        metavar="STATS_DIR",
        help="a stats directory where some basic stats can be output")
    parser.add_argument(
        "-g",
        "--logFilename",
        dest="logFilename",
        metavar="LOG_FILE",
        help="the name of the log file, STDOUT by default")
    parser.add_argument(
        "-b",
        "--batchSize",
        type=int,
        dest="batchSize",
        default=int(250000000),
        metavar="BATCH_SIZE",
        help="the size of the samtool selections that are loaded into memory at one time, 250000000 by default")
############################################
    parser.add_argument(
        "-d",
        "--dataSource",
        dest="dataSource",
        metavar="DATA_SOURCE",
        help="the source of the data - used in the sample VCF meta tag")
    parser.add_argument(
        "-q",
        "--sequencingPlatform",
        dest="sequencingPlatform",
        metavar="SEQ_PLATFORM",
        help="the sequencing platform - used in the sample VCF meta tag")
    parser.add_argument(
        "--disease",
        dest="disease",
        metavar="DISEASE",
        help="a disease abbreviation (i.e. BRCA) for the header")
    parser.add_argument(
        "--genotypeMinDepth",
        type=int,
        default=int(2),
        dest="genotypeMinDepth",
        metavar="GT_MIN_DP",
        help="the minimum number of bases required for the genotype, default 2")
    parser.add_argument(
        "--genotypeMinPct",
        type=float,
        default=float(.10),
        dest="genotypeMinPct",
        metavar="GT_MIN_PCT",
        help="the minimum percentage of reads required for the genotype, default .10")
    parser.add_argument(
        "--gzip",
        action="store_true",
        default=False,
        dest="gzip",
        help="include this argument if the final VCF should be compressed with gzip")

    # params for normal DNA
    parser.add_argument(
        "-n",
        "--dnaNormalFilename",
        dest="dnaNormalFilename",
        metavar="DNA_NORMAL_FILE",
        help="the name of the normal DNA .bam file")
    parser.add_argument(
        "--dnaNormalBaiFilename",
        dest="dnaNormalBaiFilename",
        metavar="DNA_NORMAL_BAI_FILE",
        help="the name of the normal DNA .bai file")
    parser.add_argument(
        "--dnaNormalMinTotalBases",
        type=int,
        default=int(4),
        dest="dnaNormalMinTotalBases",
        metavar="DNA_NOR_MIN_TOTAL_BASES",
        help="the minimum number of overall normal DNA reads covering a position, default 4")
    parser.add_argument(
        "--dnaNormalMinAltBases",
        type=int,
        default=int(2),
        dest="dnaNormalMinAltBases",
        metavar="DNA_NOR_MIN_ALT_BASES",
        help="the minimum number of alternative normal DNA reads supporting a variant at a position, default 2")
    parser.add_argument(
        "--dnaNormalBaseQual",
        type=int,
        default=int(10),
        dest="dnaNormalMinBaseQual",
        metavar="DNA_NOR_BASE_QUAL",
        help="the minimum normal DNA base quality, default 10")
    parser.add_argument(
        "--dnaNormalMapQual",
        type=int,
        default=int(10),
        dest="dnaNormalMinMappingQual",
        metavar="DNA_NOR_MAP_QUAL",
        help="the minimum normal DNA mapping quality, default 10")
    parser.add_argument(
        "--dnaNormalFasta",
        dest="dnaNormalFastaFilename",
        metavar="DNA_NOR_FASTA_FILE",
        help="the name of the fasta file for the normal DNA .bam file")
    parser.add_argument(
        "--dnaNormalDescription",
        default="NormalDNASample",
        dest="dnaNormalDesc",
        metavar="DNA_NOR_DESC",
        help="the description for the sample in the VCF header, default NormalDNASample")

    # params for normal RNA
    parser.add_argument(
        "-x",
        "--rnaNormalFilename",
        dest="rnaNormalFilename",
        metavar="RNA_NORMAL_FILE",
        help="the name of the normal RNA-Seq .bam file")
    parser.add_argument(
        "--rnaNormalBaiFilename",
        dest="rnaNormalBaiFilename",
        metavar="RNA_NORMAL_BAI_FILE",
        help="the name of the normal RNA .bai file")
    parser.add_argument(
        "--rnaNormalMinTotalBases",
        type=int,
        default=int(4),
        dest="rnaNormalMinTotalBases",
        metavar="RNA_NOR_MIN_TOTAL_BASES",
        help="the minimum number of overall normal RNA-Seq reads covering a position, default 4")
    parser.add_argument(
        "--rnaNormalMinAltBases",
        type=int,
        default=int(2),
        dest="rnaNormalMinAltBases",
        metavar="RNA_NOR_MIN_ALT_BASES",
        help="the minimum number of alternative normal RNA-Seq reads supporting a variant at a position, default 2")
    parser.add_argument(
        "--rnaNormalBaseQual",
        type=int,
        default=int(10),
        dest="rnaNormalMinBaseQual",
        metavar="RNA_NOR_BASE_QUAL",
        help="the minimum normal RNA-Seq base quality, default 10")
    parser.add_argument(
        "--rnaNormalMapQual",
        type=int,
        default=int(10),
        dest="rnaNormalMinMappingQual",
        metavar="RNA_NOR_MAP_QUAL",
        help="the minimum normal RNA-Seq mapping quality, default 10")
    parser.add_argument(
        "--rnaNormalFasta",
        dest="rnaNormalFastaFilename",
        metavar="RNA_NOR_FASTA_FILE",
        help="the name of the fasta file for the normal RNA .bam file")
    parser.add_argument(
        "--rnaNormalDescription",
        default="NormalRNASample",
        dest="rnaNormalDesc",
        metavar="RNA_NOR_DESC",
        help="the description for the sample in the VCF header,i default NormalRNASample")

    # params for tumor DNA
    parser.add_argument(
        "-t",
        "--dnaTumorFilename",
        dest="dnaTumorFilename",
        metavar="DNA_TUMOR_FILE",
        help="the name of the tumor DNA .bam file")
    parser.add_argument(
        "--dnaTumorBaiFilename",
        dest="dnaTumorBaiFilename",
        metavar="DNA_TUMOR_BAI_FILE",
        help="the name of the tumor DNA .bai file")
    parser.add_argument(
        "--dnaTumorMinTotalBases",
        type=int,
        default=int(4),
        dest="dnaTumorMinTotalBases",
        metavar="DNA_TUM_MIN_TOTAL_BASES",
        help="the minimum number of overall tumor DNA reads covering a position, default 4")
    parser.add_argument(
        "--dnaTumorMinAltBases",
        type=int,
        default=int(2),
        dest="dnaTumorMinAltBases",
        metavar="DNA_TUM_MIN_ALT_BASES",
        help="the minimum number of alternative tumor DNA reads supporting a variant at a position, default 2")
    parser.add_argument(
        "--dnaTumorBaseQual",
        type=int,
        default=int(10),
        dest="dnaTumorMinBaseQual",
        metavar="DNA_TUM_BASE_QUAL",
        help="the minimum tumor DNA base quality, default 10")
    parser.add_argument(
        "--dnaTumorMapQual",
        type=int,
        default=int(10),
        dest="dnaTumorMinMappingQual",
        metavar="DNA_TUM_MAP_QUAL",
        help="the minimum tumor DNA mapping quality, default 10")
    parser.add_argument(
        "--dnaTumorFasta",
        dest="dnaTumorFastaFilename",
        metavar="DNA_TUM_FASTA_FILE",
        help="the name of the fasta file for the tumor DNA .bam file")
    parser.add_argument(
        "--dnaTumorDescription",
        default="TumorDNASample",
        dest="dnaTumorDesc",
        metavar="DNA_TUM_DESC",
        help="the description for the sample in the VCF header, default TumorDNASample")

    # params for tumor RNA
    parser.add_argument(
        "-r",
        "--rnaTumorFilename",
        dest="rnaTumorFilename",
        metavar="RNA_TUMOR_FILE",
        help="the name of the tumor RNA-Seq .bam file")
    parser.add_argument(
        "--rnaTumorBaiFilename",
        dest="rnaTumorBaiFilename",
        metavar="RNA_TUMOR_BAI_FILE",
        help="the name of the tumor RNA .bai file")
    parser.add_argument(
        "--rnaTumorMinTotalBases",
        type=int,
        default=int(4),
        dest="rnaTumorMinTotalBases",
        metavar="RNA_TUM_MIN_TOTAL_BASES",
        help="the minimum number of overall tumor RNA-Seq reads covering a position, default 4")
    parser.add_argument(
        "--rnaTumorMinAltBases",
        type=int,
        default=int(2),
        dest="rnaTumorMinAltBases",
        metavar="RNA_TUM_MIN_ALT_BASES",
        help="the minimum number of alternative tumor RNA-Seq reads supporting a variant at a position, default 2")
    parser.add_argument(
        "--rnaTumorBaseQual",
        type=int,
        default=int(10),
        dest="rnaTumorMinBaseQual",
        metavar="RNA_TUM_BASE_QUAL",
        help="the minimum tumor RNA-Seq base quality, default 10")
    parser.add_argument(
        "--rnaTumorMapQual",
        type=int,
        default=int(10),
        dest="rnaTumorMinMappingQual",
        metavar="RNA_TUM_MAP_QUAL",
        help="the minimum tumor RNA-Seq mapping quality, default 10")
    parser.add_argument(
        "--rnaTumorFasta",
        dest="rnaTumorFastaFilename",
        metavar="RNA_TUM_FASTA_FILE",
        help="the name of the fasta file for the tumor RNA .bam file")
    parser.add_argument(
        "--rnaTumorDescription",
        default="TumorRNASample",
        dest="rnaTumorDesc",
        metavar="RNA_TUM_DESC",
        help="the description for the sample in the VCF header, default TumorRNASample")

    # some extra stuff
    parser.add_argument('--number_of_procs', dest='procs', type=int, default=1)
    parser.add_argument('--workdir', default="./")
    parser.add_argument('--no_clean', action="store_true", default=False)

    parser.add_argument('--exclude',
                        type=str,
                        nargs="+",
                        default=["hs37d5", "GL.*"],
                        help="chromosomes/contigs matching these patterns will be excluded from analysis")

    args = parser.parse_args()
    tempDir = tempfile.mkdtemp(dir="./", prefix="radia_work_")

    try:
        # if a universal fasta file is specified, then use it
        if (args.fastaFilename is not None):
            universalFastaFile = indexFasta(
                args.workdir, args.fastaFilename, prefix="universal")

        # if individual fasta files are specified, they over-ride the universal
        # one
        if (args.dnaNormalFastaFilename is not None):
            i_dnaNormalFastaFilename = indexFasta(
                args.workdir, args.dnaNormalFastaFilename, prefix="dnaN")
        else:
            i_dnaNormalFastaFilename = universalFastaFile
        if (args.rnaNormalFastaFilename is not None):
            i_rnaNormalFastaFilename = indexFasta(
                args.workdir, args.rnaNormalFastaFilename, prefix="rnaN")
        else:
            i_rnaNormalFastaFilename = universalFastaFile
        if (args.dnaTumorFastaFilename is not None):
            i_dnaTumorFastaFilename = indexFasta(
                args.workdir, args.dnaTumorFastaFilename, prefix="dnaT")
        else:
            i_dnaTumorFastaFilename = universalFastaFile
        if (args.rnaTumorFastaFilename is not None):
            i_rnaTumorFastaFilename = indexFasta(
                args.workdir, args.rnaTumorFastaFilename, prefix="rnaT")
        else:
            i_rnaTumorFastaFilename = universalFastaFile

        # sanity check: input bam files should all be different
        if identicalName([args.dnaNormalFilename,
                          args.dnaTumorFilename,
                          args.rnaNormalFilename,
                          args.rnaTumorFilename]):
            raise Exception("ERROR: Found duplicate input bam file")

        if (args.dnaNormalFilename is not None):
            i_dnaNormalFilename = indexBam(
                workdir=args.workdir,
                inputBamFile=args.dnaNormalFilename,
                inputBamFileIndex=args.dnaNormalBaiFilename,
                prefix="dnaNormal")
        else:
            i_dnaNormalFilename = None

        if (args.dnaTumorFilename is not None):
            i_dnaTumorFilename = indexBam(
                workdir=args.workdir,
                inputBamFile=args.dnaTumorFilename,
                inputBamFileIndex=args.dnaTumorBaiFilename,
                prefix="dnaTumor")
        else:
            i_dnaTumorFilename = None

        if (args.rnaNormalFilename is not None):
            i_rnaNormalFilename = indexBam(
                workdir=args.workdir,
                inputBamFile=args.rnaNormalFilename,
                inputBamFileIndex=args.rnaNormalBaiFilename,
                prefix="rnaNormal")
        else:
            i_rnaNormalFilename = None

        if (args.rnaTumorFilename is not None):
            i_rnaTumorFilename = indexBam(
                workdir=args.workdir,
                inputBamFile=args.rnaTumorFilename,
                inputBamFileIndex=args.rnaTumorBaiFilename,
                prefix="rnaTumor")
        else:
            i_rnaTumorFilename = None

        # clean input descriptions (this matters if we want to create TCGA
        # compliant headers in radia_filter)
        args.dnaNormalDesc = removeSpaces(args.dnaNormalDesc)
        args.dnaTumorDesc = removeSpaces(args.dnaTumorDesc)
        args.rnaNormalDesc = removeSpaces(args.rnaNormalDesc)
        args.rnaTumorDesc = removeSpaces(args.rnaTumorDesc)

        radiaOuts = []
        chroms = get_bam_seq(i_dnaNormalFilename, args.exclude)
        if args.procs == 1:
            for chrom in chroms:
                cmd, radiaOutput = radia(chrom, args, tempDir,
                                         dnaNormalFilename=i_dnaNormalFilename, rnaNormalFilename=i_rnaNormalFilename,
                                         dnaTumorFilename=i_dnaTumorFilename, rnaTumorFilename=i_rnaTumorFilename,
                                         dnaNormalFastaFilename=i_dnaNormalFastaFilename, rnaNormalFastaFilename=i_rnaNormalFastaFilename,
                                         dnaTumorFastaFilename=i_dnaTumorFastaFilename, rnaTumorFastaFilename=i_rnaTumorFastaFilename)
                if execute([cmd, radiaOutput]):
                    raise Exception("Radia Call failed")
                radiaOuts.append(radiaOutput)
        else:
            cmds = []
            for chrom in chroms:
                # create the RADIA commands
                cmd, radiaOutput = radia(chrom, args, tempDir,
                                         i_dnaNormalFilename, i_rnaNormalFilename, i_dnaTumorFilename, i_rnaTumorFilename,
                                         i_dnaNormalFastaFilename, i_rnaNormalFastaFilename, i_dnaTumorFastaFilename, i_rnaTumorFastaFilename)
                cmds.append(cmd)
                radiaOuts.append(radiaOutput)
            p = Pool(args.procs)
            # pool.map only accepts one input, so make that a list
            combiCmds = zip(cmds, radiaOuts)
            values = p.map(execute, combiCmds, 1)

        # even though we have a list of radia output files, we don't really need it:
        # the radiaMerge command only uses the output directory and patient
        # name
        cmd = radiaMerge(args, tempDir)
        if execute([cmd]):
            raise Exception("RadiaMerge Call failed")
    finally:
        if not args.no_clean and os.path.exists(tempDir):
            shutil.rmtree(tempDir)


if __name__ == "__main__":
    __main__()
