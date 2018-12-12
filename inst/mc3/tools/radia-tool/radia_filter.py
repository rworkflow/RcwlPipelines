#!/usr/bin/env python

import argparse
import os
import re
import shutil
import subprocess
import tempfile
import time
import gzip
import zipfile
import sys
from multiprocessing import Pool
from collections import OrderedDict


def execute(cmd, output=None):
    import shlex
    # function to execute a cmd and report if an error occurs
    print(cmd)
    try:
        process = subprocess.Popen(
            args=shlex.split(cmd),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
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


def move(avant, apres):
    if os.path.exists(avant):
        execute("mv %s %s" % (avant, apres))


def correctLineCount(number, vcfFile):
    """Checks number of non header lines in input VCF."""
    f = get_read_fileHandler(vcfFile)
    count = 0
    for line in f:
        if line.startswith('#'):
            continue
        count += 1
    f.close()
    if count != number:
        sys.stdout.write(
            "ERROR expected %s non header lines in %s, got %s\n" %
            (vcfFile, count, number))
        return False
    return True


def indexBam(workdir, prefix, inputBamFile, inputBamFileIndex=None):
    inputBamLink = os.path.join(os.path.abspath(workdir), prefix + ".bam")
    os.symlink(inputBamFile, inputBamLink)
    if os.path.exists(inputBamFile + ".bai"):
        inputBamFileIndex = inputBamFile + ".bai"
    if inputBamFileIndex is None:
        cmd = "samtools index %s" % (inputBamLink)
        execute(cmd)
    else:
        os.symlink(inputBamFileIndex, inputBamLink + ".bai")
    return inputBamLink


def indexFasta(workdir, inputFastaFile,
               inputFastaFileIndex=None, prefix="dna"):
    """Checks if fasta index exists. If so, creates link. If not, creates index"""
    inputFastaLink = os.path.join(
        os.path.abspath(workdir),
        prefix + "_reference.fa")
    os.symlink(inputFastaFile, inputFastaLink)
    inputFastaFileIndex = inputFastaFile + ".fai"
    if os.path.exists(inputFastaFileIndex):
        os.symlink(inputFastaFileIndex, inputFastaLink + ".fai")
    else:
        cmd = "samtools faidx %s" % (inputFastaLink)
        execute(cmd)
    return inputFastaLink


def get_read_fileHandler(aFilename):
    """ Open aFilename for reading and return the file handler.  The file can be gzipped or not."""
    if aFilename.endswith('.gz'):
        return gzip.open(aFilename, 'rb')
    else:
        return open(aFilename, 'r')


def get_write_fileHandler(aFilename):
    """ Open aFilename for writing and return the file handler.  The file can be gzipped or not."""
    if aFilename.endswith('.gz'):
        return gzip.open(aFilename, 'wb')
    else:
        return open(aFilename, 'w')


def rewriteVcfGenerator(vcfline, files):
    """Replace filenames in vcfGenerator field to match local files."""
    fields = vcfline.split(',')
    for i in xrange(len(fields)):
        try:
            key, value = fields[i].split('=')
        except:
            continue
        if key == 'dnaNormalFilename':
            if files.dnaNormalFilename is None:
                sys.stderr.write(
                    "VCF header contains DNA normal bam, please input corresponding bam file")
                sys.exit(1)
            fields[i] = ('').join([key, "=<", files.dnaNormalFilename, ">"])
        elif key == 'dnaTumorFilename':
            if files.dnaTumorFilename is None:
                sys.stderr.write(
                    "VCF header contains DNA tumor bam, please input corresponding bam file")
                sys.exit(1)
            fields[i] = ('').join([key, "=<", files.dnaTumorFilename, ">"])
        elif key == 'rnaNormalFilename':
            if files.rnaNormalFilename is None:
                sys.stderr.write(
                    "VCF header contains RNA normal bam, please input corresponding bam file")
                sys.exit(1)
            fields[i] = ('').join([key, "=<", files.rnaNormalFilename, ">"])
        elif key == 'rnaTumorFilename':
            if files.rnaTumorFilename is None:
                sys.stderr.write(
                    "VCF header contains RNA tumor bam, please input corresponding bam file")
                sys.exit(1)
            fields[i] = ('').join([key, "=<", files.rnaTumorFilename, ">"])
        elif key == 'dnaNormalFastaFilename':
            if files.dnaNormalFastaFilename is None:
                sys.stderr.write(
                    "VCF header contains DNA normal fasta, please input corresponding fasta file")
                sys.exit(1)
            fields[i] = ('').join(
                [key, "=<", files.dnaNormalFastaFilename, ">"])
        elif key == 'dnaTumorFastaFilename':
            if files.dnaTumorFastaFilename is None:
                sys.stderr.write(
                    "VCF header contains DNA tumor fasta, please input corresponding fasta file")
                sys.exit(1)
            fields[i] = ('').join(
                [key, "=<", files.dnaTumorFastaFilename, ">"])
        elif key == 'rnaNormalFastaFilename':
            if files.rnaNormalFastaFilename is None:
                sys.stderr.write(
                    "VCF header contains RNA normal fasta, please input corresponding fasta file")
                sys.exit(1)
            fields[i] = ('').join(
                [key, "=<", files.rnaNormalFastaFilename, ">"])
        elif key == 'rnaTumorFastaFilename':
            if files.rnaTumorFastaFilename is None:
                sys.stderr.write(
                    "VCF header contains RNA tumor fasta, please input corresponding fasta file")
                sys.exit(1)
            fields[i] = ('').join(
                [key, "=<", files.rnaTumorFastaFilename, ">"])
    newline = (',').join(fields)
    return newline


def splitVcf(infile, outdir, files=None, expected=None):
    """Splits up VCF file in chromosome files, a dict of chromosome names with non header linecounts and a dict of chromosome names (without chr) and corresponding vcf files (with chr)"""
    chrNames = dict()
    chrLines = dict()
    header = ""
    f = get_read_fileHandler(infile)
    headFlag = True
    for line in f:
        if headFlag:
            if files is not None and line.startswith('##vcfGenerator'):
                line = rewriteVcfGenerator(line, files)
            header += line
            if line.startswith('#CHROM'):
                headFlag = False
        else:
            fields = line.split("\t")
            chrom = fields[0].replace('chr', '')  # remove chr if present
            if chrom not in chrNames:
                try:
                    o.close()
                except:
                    pass
                outfile = os.path.join(outdir, 'chr' + chrom + ".vcf.gz")
                o = get_write_fileHandler(outfile)      # append not necessary
                o.write(header)
                chrLines[chrom] = 0
                chrNames[chrom] = outfile
            o.write(line)
            chrLines[chrom] += 1
    o.close
    f.close()
    if expected is not None:
        wanted = expected.keys()
        wanted = [s.replace('chr', '') for s in wanted]  # remove chr
        created = set(chrNames.keys())
        for chrom in set(wanted).difference(created):
            sys.stderr.write(
                "WARNING, missing chromosome %s in filter %s, ignoring...\n" %
                (chrom, outdir))
            # creating empty file
            outfile = os.path.join(outdir, "chr" + chrom + ".vcf")
            o = get_write_fileHandler(outfile)
            o.close
    return chrNames, chrLines


def splitBed(bedfile, outdir, chromDict):
    """Splits bed file in chromosome files, issued warnings on missing files and creates empties"""
    wanted = chromDict.keys()
    wanted = [s.replace('chr', '') for s in wanted]  # remove chr
    wanted = ['chr' + s for s in wanted]  # add chr
    chrNames = set()
    header = ""
    f = get_read_fileHandler(bedfile)
    for line in f:
        fields = line.split("\t")
        chrom = fields[0]
        if not chrom.startswith('chr'):
            chrom = 'chr' + chrom
        if chrom not in chrNames:
            try:
                o.close()
            except:
                pass
            outfile = os.path.join(outdir, chrom + ".bed.gz")
            o = get_write_fileHandler(outfile)      # append not necessary
            chrNames.add(chrom)
        o.write(line)
    o.close
    f.close()
    for chrom in set(wanted).difference(chrNames):
        sys.stderr.write(
            "WARNING, missing chromosome %s in filter %s, ignoring...\n" %
            (chrom, outdir))
        # creating empty file
        outfile = os.path.join(outdir, chrom + ".bed")
        o = get_write_fileHandler(outfile)
        o.close


def identicalName(inputList):
    """returns duplicate name if two inputs have the same name and are not None"""
    dup = set(x for x in inputList if inputList.count(x) >= 2)
    dup.discard(None)   # this doesn't complain if None is not in the set
    if dup:
        print "ERROR: found duplicate input %s" % dup.pop()
        return True
    return False


def makeSnpEffConfig(workdir, genome, datadir):
    """Creates a short config file for snpEff. Assumes a human genome."""
    configFile = os.path.join(workdir, "snpEff.config")
    f = open(configFile, 'w')
    f.write("data.dir = %s\n" % datadir)
    f.write("lof.ignoreProteinCodingAfter : 0.95\n")
    f.write("lof.ignoreProteinCodingBefore : 0.05\n")
    f.write("lof.deleteProteinCodingBases : 0.50\n")
    f.write("codon.Standard : TTT/F, TTC/F, TTA/L, TTG/L+, TCT/S, TCC/S, TCA/S, TCG/S, TAT/Y, TAC/Y, TAA/*, TAG/*, TGT/C, TGC/C, TGA/*, TGG/W, CTT/L, CTC/L, CTA/L, CTG/L+, CCT/P, CCC/P, CCA/P, CCG/P, CAT/H, CAC/H, CAA/Q, CAG/Q, CGT/R, CGC/R, CGA/R, CGG/R, ATT/I, ATC/I, ATA/I, ATG/M+, ACT/T, ACC/T, ACA/T, ACG/T, AAT/N, AAC/N, AAA/K, AAG/K, AGT/S, AGC/S, AGA/R, AGG/R, GTT/V, GTC/V, GTA/V, GTG/V, GCT/A, GCC/A, GCA/A, GCG/A, GAT/D, GAC/D, GAA/E, GAG/E, GGT/G, GGC/G, GGA/G, GGG/G\n")
    f.write("codon.Vertebrate_Mitochondrial : TTT/F, TTC/F, TTA/L, TTG/L, TCT/S, TCC/S, TCA/S, TCG/S, TAT/Y, TAC/Y, TAA/*, TAG/*, TGT/C, TGC/C, TGA/W, TGG/W, CTT/L, CTC/L, CTA/L, CTG/L, CCT/P, CCC/P, CCA/P, CCG/P, CAT/H, CAC/H, CAA/Q, CAG/Q, CGT/R, CGC/R, CGA/R, CGG/R, ATT/I+, ATC/I+, ATA/M+, ATG/M+, ACT/T, ACC/T, ACA/T, ACG/T, AAT/N, AAC/N, AAA/K, AAG/K, AGT/S, AGC/S, AGA/*, AGG/*, GTT/V, GTC/V, GTA/V, GTG/V+, GCT/A, GCC/A, GCA/A, GCG/A, GAT/D, GAC/D, GAA/E, GAG/E, GGT/G, GGC/G, GGA/G, GGG/G\n")
    f.write("%s.genome : Homo_sapiens\n" % genome)
    f.close()
    return configFile


def radiaFilter(filterDirs, snpEffGenome, snpEffConfig,
                args, chrom, inputVcf, outputDir, logFile):

    # python filterRadia.pyc id chrom inputFile outputDir scriptsDir [Options]

    # python filterRadia.pyc TCGA-02-0047-10A-01D-1490-08_TCGA-02-0047-01A-01D-1490-08 1
    # TCGA-02-0047-10A-01D-1490-08_TCGA-02-0047-01A-01D-1490-08_chr1.vcf.gz
    # /radia/finalChromVCFs/
    # /rnaEditing/scripts/
    # --blacklistDir /rnaEditing/data/hg19/blacklists/1000Genomes/phase1/
    # --dbSnpDir /rnaEditing/data/hg19/snp135/
    # --retroGenesDir /rnaEditing/data/hg19/retroGenes/
    # --pseudoGenesDir /rnaEditing/data/hg19/pseudoGenes/
    # --cosmicDir /rnaEditing/data/hg19/cosmic/
    # --targetDir /rnaEditing/data/hg19/broadTargets/
    # --snpEffDir /snpEff/
    # --rnaGeneBlckFile /rnaEditing/data/rnaGeneBlacklist.tab
    # --rnaGeneFamilyBlckFile /rnaEditing/data/rnaGeneFamilyBlacklist.tab
    # --blatFastaFilename hg19.fasta
    # --canonical
    # --log=INFO
    # --gzip
    cmd = "python %s/filterRadia.py %s %s %s %s %s --gzip --log=WARNING -g %s" % (
        args.scriptsDir,
        args.patientId, chrom, inputVcf,
        outputDir, args.scriptsDir,
        logFile)

    if args.blatFastaFilename is not None:
        cmd += " --blatFastaFilename %s" % args.blatFastaFilename
    else:
        cmd += ' --noBlat'

    if "blacklist" in filterDirs:
        cmd += " --blacklistDir %s" % filterDirs["blacklist"]
    else:
        cmd += ' --noBlacklist'

    if "target" in filterDirs:
        cmd += " --targetDir %s" % filterDirs["target"]
    else:
        cmd += ' --noTargets'

    if "snp" in filterDirs:
        cmd += " --dbSnpDir %s" % filterDirs["snp"]
    else:
        cmd += ' --noDbSnp'

    if "pseudoGenes" in filterDirs:
        cmd += " --pseudoGenesDir %s" % filterDirs["pseudoGenes"]
    else:
        cmd += ' --noPseudoGenes'

    if "retroGenes" in filterDirs:
        cmd += " --retroGenesDir %s" % filterDirs["retroGenes"]
    else:
        cmd += ' --noRetroGenes'

    if "cosmic" in filterDirs:
        cmd += " --cosmicDir %s" % filterDirs["cosmic"]
    else:
        cmd += ' --noCosmic'

    if snpEffGenome is not None:
        cmd += " --snpEffDir %s --snpEffGenome %s --snpEffConfig %s" % (
            args.snpEffDir, snpEffGenome, snpEffConfig)
        if args.canonical:
            cmd += ' --canonical '
        if args.rnaGeneBlckFile:
            cmd += ' --rnaGeneBlckFile %s --rnaGeneFamilyBlckFile %s' % (
                args.rnaGeneBlckFile, args.rnaGeneFamilyBlckFile)
        else:
            cmd += ' --noRnaBlacklist'
    else:
        cmd += ' --noSnpEff --noRnaBlacklist'

    if args.noPositionalBias:
        cmd += ' --noPositionalBias'
    if args.dnaOnly:
        cmd += ' --dnaOnly'
#    if args.rnaOnly:
#        cmd += ' --rnaOnly '
#    if args.gzip:
#        cmd += ' --gzip '
    outfile = os.path.join(
        outputDir,
        args.patientId +
        "_chr" +
        chrom +
        ".vcf.gz")
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
    mergeOutput = os.path.join(args.workdir, "mergeOut.vcf")
    cmd = "python %s/mergeChroms.py %s %s %s -o %s" % (
        args.scriptsDir,
        args.patientId, inputDir, args.workdir,
        mergeOutput)
#    if args.gzip:
#        cmd += ' --gzip'
    return cmd, mergeOutput


class localFiles(object):
    """Formats input fasta and bam files and creates local filenames"""

    def __init__(self):
        self.dnaNormalFilename, self.dnaTumorFilename, self.rnaNormalFilename, self.rnaTumorFilename, self.dnaNormalFastaFilename, self.dnaTumorFastaFilename, self.rnaNormalFastaFilename, self.rnaTumorFastaFilename = (
            None for i in xrange(8))

    def universalFasta(self, args):
        if (args.fastaFilename is not None):
            universalFastaFile = indexFasta(
                args.workdir, args.fastaFilename, prefix="universal")
            self.dnaNormalFastaFilename, self.dnaTumorFastaFilename, self.rnaNormalFastaFilename, self.rnaTumorFastaFilename = (
                universalFastaFile for i in xrange(4))

    def doFasta(self, args):
        # if individual fasta files are specified, they over-ride the universal
        # one
        if (args.dnaNormalFastaFilename is not None):
            self.dnaNormalFastaFilename = indexFasta(
                args.workdir, args.dnaNormalFastaFilename, prefix="dnaN")

        if (args.rnaNormalFastaFilename is not None):
            self.rnaNormalFastaFilename = indexFasta(
                args.workdir, args.rnaNormalFastaFilename, prefix="rnaN")
        if (args.dnaTumorFastaFilename is not None):
            self.dnaTumorFastaFilename = indexFasta(
                args.workdir, args.dnaTumorFastaFilename, prefix="dnaT")
        if (args.rnaTumorFastaFilename is not None):
            self.rnaTumorFastaFilename = indexFasta(
                args.workdir, args.rnaTumorFastaFilename, prefix="rnaT")

    def doBam(self, args):
        # index bam files and return True if there's only DNA files
        if (args.dnaNormalFilename is not None):
            self.dnaNormalFilename = indexBam(
                workdir=args.workdir,
                inputBamFile=args.dnaNormalFilename,
                inputBamFileIndex=args.dnaNormalBaiFilename,
                prefix="dnaNormal")

        if (args.rnaNormalFilename is not None):
            self.rnaNormalFilename = indexBam(
                workdir=args.workdir,
                inputBamFile=args.rnaNormalFilename,
                inputBamFileIndex=args.rnaNormalBaiFilename,
                prefix="rnaNormal")

        if (args.dnaTumorFilename is not None):
            self.dnaTumorFilename = indexBam(
                workdir=args.workdir,
                inputBamFile=args.dnaTumorFilename,
                inputBamFileIndex=args.dnaTumorBaiFilename,
                prefix="dnaTumor")

        if (args.rnaTumorFilename is not None):
            self.rnaTumorFilename = indexBam(
                workdir=args.workdir,
                inputBamFile=args.rnaTumorFilename,
                inputBamFileIndex=args.rnaTumorBaiFilename,
                prefix="rnaTumor")
        if (args.rnaTumorFilename is None and args.rnaNormalFilename is None):
            return True
        return False

#############################
# TCGA compliance modules   #
#############################


class VCFFormatError(Exception):

    def __init__(self, text):
        Exception.__init__(self, text)


class Info:
    """Object that holds information on VCF INFO header tag"""

    def __init__(self, sampleid, number, sampletype, description):
        self.id = sampleid
        self.number = number
        self.type = sampletype
        self.description = description

    def __str__(self):
        return '##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % (
            self.id, self.number, self.type, self.description)


def parse_info(info):
    """Extracts keys and values from VCF INFO header tag, returns dict"""
    info_dict = OrderedDict()
    for i in info.split(";"):
        if "=" in i:
            key, value = i.split("=")
            info_dict[key] = value
        else:
            info_dict[i] = True
    return info_dict


def format_info(info_dict):
    """Formats INFO tag values for printing to VCF header"""
    line_list = []
    for key, value in info_dict.items():
        if value is True:
            line_list.append(key)
        else:
            line_list.append("%s=%s" % (key, value))
    return ";".join(line_list)


class Data:
    """Holds VCF body data"""

    def __init__(self, chrom, pos, sampleid, ref, alt, qual,
                 filters, info, genotype_data=None, extra_headers=None):
        self.chrom = chrom
        self.pos = int(pos)
        self.id = sampleid.split(";")
        self.ref = ref
        self.alt = alt.split(",")
        self.qual = int(float(qual))
        self.filter = filters.split(";")
        self.info = parse_info(info)

        self.genotype_data_order = genotype_data[0].split(":")
        self.samples = extra_headers[1:]

        self.genotype_data = OrderedDict()
        for sample, genotype_values in zip(self.samples, genotype_data[1:]):
            self.genotype_data[sample] = OrderedDict()
            for i, g in zip(self.genotype_data_order,
                            genotype_values.split(":")):
                self.genotype_data[sample][i] = g
        #self.genotype = genotype_data

    def add_genotype_data(self, label, value_list):
        #assert(len(value_list) == (len(self.samples)))
        self.genotype_data_order.append(label)
        for sample, value in zip(self.samples, value_list):
            self.genotype_data[sample][label] = value

    def __str__(self):
        output = map(
            str, [
                self.chrom, self.pos, ";".join(
                    self.id), self.ref, ",".join(
                    self.alt), self.qual, ";".join(
                    self.filter), format_info(
                        self.info)])
        output.append(":".join(self.genotype_data_order))
        for sample in self.samples:
            output.append(":".join(self.genotype_data[sample].values()))
        return "\t".join(output)


class VCF:
    """Holds VCF info (header and body)"""

    def __init__(self):
        self.meta = []
        self.infos = []
        self.filters = []
        self.formats = []
        self.headers = []
        self.data = []

    def make_info(self, value):

        # TCGA doesn't allow a space before the final quote, and that is the default format from SnpEff
        # so remove the space here
        value = value.rstrip("\r\n")
        if (value.endswith(' ">')):
            value = value.replace(' ">', '">')

        # make re match spec
        match = re.match(
            r"""<ID=(.*),Number=(.*),Type=(.*),Description=['"](.*)['"]""",
            value)
        if match:
            return Info(*match.groups())
        else:
            raise VCFFormatError(
                "Improperly formatted INFO metadata: " + value)

    def set_headers(self, header_list):
        self.headers = header_list

    def make_data(self, data_list):
        baseData = data_list[:8]
        genotypeData = data_list[8:]
        data = baseData + [genotypeData] + [self.headers[8:]]
        return Data(*data)


def format_vcf(filename, outFile, filterRejects, filterGermline):
    """Reformats VCF output to comply with TCGA specs"""
    if (filename.endswith(".gz")):
        vcf_file = gzip.open(filename, 'rb')
    else:
        vcf_file = open(filename, "r")
    vcf_out = open(outFile, 'w')
    currVCF = VCF()
    line = ""

    rsid_dict = {}

    # first deal with the header info
    for line in vcf_file:
        # if it doesn't start with ##, then we're done with the header, so
        # break out
        if not line.startswith("##"):
            break

        key, value = line[2:].split("=", 1)
        if key == "tcgaversion":
            vcf_out.write("##tcgaversion=1.1\n")
        elif key == "INDIVIDUAL":
            continue
        elif key == "INFO":
            info = currVCF.make_info(value)
            if info.id == "Gene":
                continue
            elif info.id == "VC":
                continue
            elif info.id == "SS":
                info.type = "Integer"
                vcf_out.write(str(info))
            elif info.id == "VT":
                info.description = "Variant type, can be SNP, INS or DEL"
                vcf_out.write(str(info))
            elif info.id == "DP":
                info.description = "Total Depth across samples"
                vcf_out.write(str(info))
            elif info.id == "SOMATIC":
                info.description = "Indicates if record is a somatic mutation"
                vcf_out.write(str(info))
            elif info.id == "DB":
                info.description = "dbSNP membership"
                vcf_out.write(str(info))
            elif info.id == "NS":
                info.description = "Number of Samples With Data"
                vcf_out.write(str(info))
            elif info.id == "AN":
                info.description = "Total number of alleles in called genotypes"
                vcf_out.write(str(info))
            elif info.id == "AF":
                info.description = "Allele Frequency in primary data, for each ALT allele, in the same order as listed"
                vcf_out.write(str(info))
            elif info.id == "BQ":
                info.description = "RMS base quality"
                info.type = "Integer"
                vcf_out.write(str(info))
            elif info.id == "SB":
                info.description = "Strand bias"
                vcf_out.write(str(info))
            else:
                vcf_out.write(str(info))
        elif key == "FILTER":
            if "ID=perfectsbias" in value:
                vcf_out.write(
                    '##FILTER=<ID=perfectsbias,Description="A strand bias exists on the perfect reads.">' +
                    "\n")
            else:
                vcf_out.write(line)
        elif key == "FORMAT":
            if "ID=AD" in value:
                vcf_out.write(
                    '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Depth of reads supporting alleles 0/1/2/3...">' +
                    "\n")
            elif "ID=MQ" in value:
                vcf_out.write(
                    '''##FORMAT=<ID=MQ,Number=1,Type=Integer,Description="Phred style probability score that the variant is novel with respect to the genome's ancestor">''' +
                    "\n")
                vcf_out.write(
                    '##FORMAT=<ID=MQA,Number=.,Type=Float,Description="Average mapping quality for reads supporting alleles">' +
                    "\n")
            elif "ID=BQ" in value:
                vcf_out.write(
                    '##FORMAT=<ID=BQ,Number=.,Type=Integer,Description="Average base quality for reads supporting alleles">' +
                    "\n")
            elif "ID=DP" in value:
                vcf_out.write(
                    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth at this position in the sample">' +
                    "\n")
            elif "ID=SSC" in value:
                continue
            elif "ID=SS" in value:
                continue
            else:
                vcf_out.write(line)
        else:
            vcf_out.write(line)
    vcf_out.write('##INFO=<ID=SS,Number=1,Type=Integer,Description="Variant status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post-transcriptional modification,5=unknown">' + "\n")
    vcf_out.write('##FORMAT=<ID=SS,Number=1,Type=Integer,Description="Variant status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post-transcriptional modification,5=unknown">' + "\n")
    vcf_out.write(
        '##FORMAT=<ID=SSC,Number=1,Type=Integer,Description="Somatic score between 0 and 255">' +
        "\n")

    if line.startswith("#"):
        headers = line[1:].strip().split("\t")
        currVCF.set_headers(headers)
        if headers[:8] != ["CHROM", "POS", "ID",
                           "REF", "ALT", "QUAL", "FILTER", "INFO"]:
            print >> sys.stderr, "Traceback: ERROR: The headers don't seem to be right, expected " + \
                '["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"]' + " but got " + str(headers)
            sys.exit(1)
    else:
        print >> sys.stderr, "Traceback: ERROR: expected # followed by headers"
        sys.exit(1)

    vcf_out.write(line)

    # Now look at the vcf body
    for line in vcf_file:
        dataAsList = line.strip().split("\t")

        # note from Amie:
        # this is awful.  there is sth wrong with the bambam and radia merging script so that
        # some lines don't have the proper header (excluding the RNA_TUMOR header) and then
        # the data has a . in the RNA_TUMOR column.  as a hack, add 1 to length
        # here.
        if (len(dataAsList) == len(headers) or len(
                dataAsList) == (len(headers) + 1)):
            curr_data = currVCF.make_data(dataAsList)

            # handle SS, how should I do wildtype?
            if curr_data.info["VT"] == "LOH":
                curr_data.info["SS"] = "3"
            elif curr_data.info["SS"] == "Germline":
                curr_data.info["SS"] = "1"
            elif curr_data.info["SS"] == "Somatic":
                curr_data.info["SS"] = "2"

            # handle VT
            if curr_data.info["VT"] in ["SDel", "GDel", "Del"]:
                curr_data.info["VT"] = "DEL"
            if curr_data.info["VT"] in ["SIns", "GIns", "Ins"]:
                curr_data.info["VT"] = "INS"
            elif curr_data.info["VT"] == "LOH":
                if len(curr_data.ref) > len(curr_data.alt[0]):
                    curr_data.info["VT"] = "DEL"
                if len(curr_data.ref) == len(curr_data.alt[0]):
                    curr_data.info["VT"] = "SNP"
                if len(curr_data.ref) < len(curr_data.alt[0]):
                    curr_data.info["VT"] = "INS"

            # add SS to genotype information
            if "RNA_NORMAL" in curr_data.genotype_data:
                curr_data.add_genotype_data(
                    "SS", ["0", "0", curr_data.info["SS"], curr_data.info["SS"]])
            elif "RNA_TUMOR" in curr_data.genotype_data:
                curr_data.add_genotype_data(
                    "SS", ["0", curr_data.info["SS"], curr_data.info["SS"]])
            else:
                curr_data.add_genotype_data("SS", ["0", curr_data.info["SS"]])

            # handle ID
            for rsid in curr_data.id:
                if rsid != ".":
                    if rsid in rsid_dict:
                        if rsid_dict[rsid] == "many":
                            curr_data.id[curr_data.id.index(
                                rsid)] = rsid + "_%s_%d" % (curr_data.chrom, curr_data.pos)
                        else:
                            # only one other has been found, so we need to
                            # change the original as well
                            prevdata = rsid_dict[rsid]
                            rsid_dict[rsid] = "many"
                            prevdata.id[prevdata.id.index(
                                rsid)] = rsid + "_%s_%d" % (prevdata.chrom, prevdata.pos)
                            curr_data.id[curr_data.id.index(
                                rsid)] = rsid + "_%s_%d" % (curr_data.chrom, curr_data.pos)
                    else:
                        rsid_dict[rsid] = curr_data
            # delete Gene
            if "Gene" in curr_data.info:
                del curr_data.info["Gene"]

            # delete VC
            if "VC" in curr_data.info:
                del curr_data.info["VC"]

            # fix y chromosome stuff
            if curr_data.chrom == "Y":
                for sample in curr_data.genotype_data:
                    if '/' in curr_data.genotype_data[sample]["GT"]:
                        g = curr_data.genotype_data[sample]["GT"]
                        left, right = re.match(r'(\d)/(\d)', g).groups()
                        curr_data.genotype_data[sample][
                            "GT"] = max(left, right)

            # if there is no info, make sure it is ./.
            if curr_data.chrom != "Y":
                for sample in curr_data.genotype_data:
                    if curr_data.genotype_data[sample]["GT"] == ".":
                        curr_data.genotype_data[sample]["GT"] = "./."

            # add somatic score
            if "RNA_NORMAL" in curr_data.genotype_data:
                curr_data.add_genotype_data("SSC", [str(curr_data.qual), str(
                    curr_data.qual), str(curr_data.qual), str(curr_data.qual)])
            elif "RNA_TUMOR" in curr_data.genotype_data:
                curr_data.add_genotype_data(
                    "SSC", [str(curr_data.qual), str(curr_data.qual), str(curr_data.qual)])
            else:
                curr_data.add_genotype_data(
                    "SSC", [str(curr_data.qual), str(curr_data.qual)])

            # add MQA and fix MQ
            if "MQ" in curr_data.genotype_data["DNA_NORMAL"]:
                if "RNA_NORMAL" in curr_data.genotype_data and "MQ" in curr_data.genotype_data[
                        "RNA_NORMAL"]:
                    curr_data.add_genotype_data("MQA",
                                                [curr_data.genotype_data["DNA_NORMAL"]["MQ"],
                                                 curr_data.genotype_data["RNA_NORMAL"]["MQ"],
                                                    curr_data.genotype_data["DNA_TUMOR"]["MQ"],
                                                    curr_data.genotype_data["RNA_TUMOR"]["MQ"]])
                elif "RNA_TUMOR" in curr_data.genotype_data and "MQ" in curr_data.genotype_data["RNA_TUMOR"]:
                    curr_data.add_genotype_data("MQA",
                                                [curr_data.genotype_data["DNA_NORMAL"]["MQ"],
                                                 curr_data.genotype_data["DNA_TUMOR"]["MQ"],
                                                    curr_data.genotype_data["RNA_TUMOR"]["MQ"]])
                else:
                    curr_data.add_genotype_data("MQA", [curr_data.genotype_data["DNA_NORMAL"][
                                                "MQ"], curr_data.genotype_data["DNA_TUMOR"]["MQ"]])

            for sample in curr_data.genotype_data:
                if "MQ" not in curr_data.genotype_data[
                        sample] or curr_data.genotype_data[sample]["MQ"] == ".":
                    continue
                depths = map(int, curr_data.genotype_data[
                             sample]["DP"].split(","))
                quals = map(float, curr_data.genotype_data[
                            sample]["MQ"].split(","))
                tot = 0
                for depth, qual in zip(depths, quals):
                    tot += depth * qual
                #curr_data.genotype_data[sample]["MQ"] = str(int(tot / sum(depths)))
                curr_data.genotype_data[sample]["MQ"] = "0"

            # make BQ int
            for sample in curr_data.genotype_data:
                if "BQ" not in curr_data.genotype_data[
                        sample] or curr_data.genotype_data[sample]["BQ"] == ".":
                    continue
                quals = map(
                    int, map(
                        float, curr_data.genotype_data[sample]["BQ"].split(",")))
                curr_data.genotype_data[sample][
                    "BQ"] = ",".join(map(str, quals))
            if "BQ" in curr_data.info:
                curr_data.info["BQ"] = str(int(float(curr_data.info["BQ"])))

            # fix rtbias
            if "rtbias" in curr_data.filter:
                curr_data.filter.remove("rtbias")
                curr_data.filter.append("rtsbias")
            # empty out a set of genotype data if we don't know the genotype
            for sample in curr_data.genotype_data:
                if curr_data.genotype_data[sample][
                        "GT"] == "." or curr_data.genotype_data[sample]["GT"] == "./.":
                    gt = OrderedDict()
                    for formatItem in curr_data.genotype_data_order:
                        gt[formatItem] = "."
                    curr_data.genotype_data[sample] = gt

            if not filterRejects or "PASS" in curr_data.filter:
                if not filterGermline or curr_data.info["SS"] != "1":
                    vcf_out.write(str(curr_data) + "\n")


def __main__():
    # small hack, sometimes it seems like docker file systems are avalible
    # instantly
    time.sleep(1)
    parser = argparse.ArgumentParser(description="RADIA filter")

    #############################
    #    RADIA filter params    #
    #############################
    parser.add_argument(
        "--inputVCF",
        dest="inputVCF",
        required=True,
        metavar="INPUT_VCF",
        help="The input Radia vcf file")
    parser.add_argument(
        "--patientId",
        dest="patientId",
        required=True,
        metavar="PATIENT_ID",
        help="a unique patient Id that will be used to name the output file")
    parser.add_argument(
        "-o",
        "--outputFilename",
        dest="outputFilename",
        default="filtered_out.vcf",
        metavar="OUTPUT_FILE",
        help="the name of the output file (filtered_out.vcf)")
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
        "-f",
        "--fastaFilename",
        dest="fastaFilename",
        metavar="FASTA_FILE",
        help="the name of the fasta file that can be used on all .bams, see below for specifying individual fasta files for each .bam file")
    parser.add_argument(
        "--makeTCGAcompliant",
        action="store_true",
        default=False,
        dest="makeTCGAcompliant",
        help="Change VCF to make TCGA v1.1 compliant")
    parser.add_argument(
        "--filter-rejects",
        action="store_true",
        default=False,
        dest="filterRejects",
        help="Filter out rejected calls")
    parser.add_argument(
        "--filter-germline",
        action="store_true",
        default=False,
        dest="filterGermline",
        help="Filter out germline calls")
    # normal DNA
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
        "--dnaNormalFastaFilename",
        dest="dnaNormalFastaFilename",
        metavar="DNA_NORMAL_FASTA_FILE",
        help="the name of the fasta file that was used to create the BAM alignments")

    # tumor DNA
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
        "--dnaTumorFastaFilename",
        dest="dnaTumorFastaFilename",
        metavar="DNA_TUMOR_FASTA_FILE",
        help="the name of the fasta file that was used to create the BAM alignments")

    # normal RNA
    parser.add_argument(
        "-x",
        "--rnaNormalFilename",
        dest="rnaNormalFilename",
        metavar="RNA_NORMAL_FILE",
        help="the name of the normal RNA .bam file")
    parser.add_argument(
        "--rnaNormalBaiFilename",
        dest="rnaNormalBaiFilename",
        metavar="RNA_NORMAL_BAI_FILE",
        help="the name of the normal RNA .bai file")
    parser.add_argument(
        "--rnaNormalFastaFilename",
        dest="rnaNormalFastaFilename",
        metavar="RNA_NORMAL_FASTA_FILE",
        help="the name of the fasta file that was used to create the BAM alignments")

    # tumor RNA
    parser.add_argument(
        "-r",
        "--rnaTumorFilename",
        dest="rnaTumorFilename",
        metavar="RNA_TUMOR_FILE",
        help="the name of the tumor RNA .bam file")
    parser.add_argument(
        "--rnaTumorBaiFilename",
        dest="rnaTumorBaiFilename",
        metavar="RNA_TUMOR_BAI_FILE",
        help="the name of the tumor RNA .bai file")
    parser.add_argument(
        "--rnaTumorFastaFilename",
        dest="rnaTumorFastaFilename",
        metavar="RNA_TUMOR_FASTA_FILE",
        help="the name of the fasta file that was used to create the BAM alignments")

    parser.add_argument(
        "--blacklistFilename",
        dest="blacklistFilename",
        metavar="BLACKLIST_FILE",
        help="the name of the blacklist bed file")

    parser.add_argument(
        "--targetFilename",
        dest="targetFilename",
        metavar="TARGET_FILE",
        help="the name of the exon capture targets file")
    parser.add_argument(
        "--snpFilename",
        dest="snpFilename",
        metavar="SNP_FILE",
        help="dbSNP vcf file")
    parser.add_argument(
        "--retroGenesFilename",
        dest="retroGenesFilename",
        metavar="RETRO_FILE",
        help="the name of the retrogenes bed file")
    parser.add_argument(
        "--pseudoGenesFilename",
        dest="pseudoGenesFilename",
        metavar="PSEUDO_FILE",
        help="the name of the pseudogenes bed file")
    parser.add_argument(
        "--cosmicFilename",
        dest="cosmicFilename",
        metavar="COSMIC_FILE",
        help="the name of the Catalogue Of Somatic Mutations In Cancer (COSMIC) annotations file")
    parser.add_argument(
        "--snpEffDir",
        dest="snpEffDir",
        default='/opt/snpEff',
        metavar="SNP_EFF_DIR",
        help="the path to the snpEff directory")
    parser.add_argument(
        "--snpEffFilename",
        dest="snpEffFilename",
        metavar="SNP_EFF_FILE",
        help="the snpEff input database zip file")
    parser.add_argument(
        "--canonical",
        action="store_true",
        default=False,
        dest="canonical",
        help="include this argument if only the canonical transcripts from snpEff should be used")
    parser.add_argument(
        "--rnaGeneBlckFile",
        dest="rnaGeneBlckFile",
        metavar="RNA_GENE_FILE",
        help="the RNA gene blacklist file")
    parser.add_argument(
        "--rnaGeneFamilyBlckFile",
        dest="rnaGeneFamilyBlckFile",
        metavar="RNA_GENE_FAMILY_FILE",
        help="the RNA gene family blacklist file")
    parser.add_argument(
        "--blatFastaFilename",
        dest="blatFastaFilename",
        metavar="FASTA_FILE",
        help="the fasta file that can be used during the BLAT filtering")
    parser.add_argument(
        "--noPositionalBias",
        action="store_false",
        default=True,
        dest="noPositionalBias",
        help="include this argument if the positional bias filter should not be applied")

    parser.add_argument(
        "--dnaOnly",
        action="store_true",
        default=False,
        dest="dnaOnly",
        help="include this argument if you only have DNA or filtering should only be done on the DNA")
#    parser.add_argument("--rnaOnly", action="store_true", default=False, dest="rnaOnly", help="include this argument if the filtering should only be done on the RNA")
#    parser.add_argument("--gzip", action="store_true", default=False, dest="gzip", help="include this argument if the final VCF should be compressed with gzip")

    # some extra stuff
    parser.add_argument('--number_of_procs', dest='procs', type=int, default=1)
    parser.add_argument('--workdir', default="./")
    parser.add_argument('--no_clean', action="store_true", default=False)

    args = parser.parse_args()
    tempDir = tempfile.mkdtemp(dir="./", prefix="radia_work_")

    # sanity checks
    if identicalName([args.dnaNormalFilename, args.dnaTumorFilename,
                      args.rnaNormalFilename, args.rnaTumorFilename]):
        raise Exception("ERROR: Found duplicate input bam file")
    if identicalName([args.blacklistFilename,
                      args.targetFilename,
                      args.retroGenesFilename,
                      args.pseudoGenesFilename,
                      args.cosmicFilename]):
        raise Exception("ERROR: Found duplicate input bed file")
    if (args.rnaGeneFamilyBlckFile and not args.rnaGeneBlckFile) or (
            args.rnaGeneBlckFile and not args.rnaGeneFamilyBlckFile):
        raise Exception("ERROR: Must input two RNA blacklist files")
    if identicalName([args.rnaGeneFamilyBlckFile, args.rnaGeneBlckFile]):
        raise Exception("ERROR: Found duplicate input RNA blacklist file")

    files = localFiles()  # prepares and holds fasta and bam files
    filterDirs = dict()		# holds filters and file locations
    try:
        files.universalFasta(args)
        files.doFasta(args)
        args.dnaOnly = files.doBam(args)
        # split vcf in chromosomes
        chromDict, chromLines = splitVcf(
            args.inputVCF, args.workdir, files=files)

        # All files come in as complete genome files, so first split them
        # split blacklist
        if (args.blacklistFilename is not None):
            blacklistDir = os.path.join(args.workdir, "blacklistDir")
            os.mkdir(blacklistDir)
            splitBed(args.blacklistFilename, blacklistDir, chromDict)
            filterDirs["blacklist"] = blacklistDir
        # split target
        if (args.targetFilename is not None):
            targetDir = os.path.join(args.workdir, "targetDir")
            os.mkdir(targetDir)
            splitBed(args.targetFilename, targetDir, chromDict)
            filterDirs["target"] = targetDir
        # split snp
        if (args.snpFilename is not None):
            snpDir = os.path.join(args.workdir, "snpDir")
            os.mkdir(snpDir)
            splitVcf(args.snpFilename, snpDir, expected=chromDict)
            filterDirs["snp"] = snpDir
        # split retrogenes
        if (args.retroGenesFilename is not None):
            retroGenesDir = os.path.join(args.workdir, "retroGenesDir")
            os.mkdir(retroGenesDir)
            splitBed(args.retroGenesFilename, retroGenesDir, chromDict)
            filterDirs["retroGenes"] = retroGenesDir
        # split pseudogenes
        if (args.pseudoGenesFilename is not None):
            pseudoGenesDir = os.path.join(args.workdir, "pseudoGenesDir")
            os.mkdir(pseudoGenesDir)
            splitBed(args.pseudoGenesFilename, pseudoGenesDir, chromDict)
            filterDirs["pseudoGenes"] = pseudoGenesDir
        # split cosmic
        if (args.cosmicFilename is not None):
            cosmicDir = os.path.join(args.workdir, "cosmicDir")
            os.mkdir(cosmicDir)
            splitBed(args.cosmicFilename, cosmicDir, chromDict)
            filterDirs["cosmic"] = cosmicDir
        # setup snpEff database
        if (args.snpEffFilename):
            with zipfile.ZipFile(args.snpEffFilename, "r") as z:
                z.extractall(args.workdir)
            # this creates a directory named data, which holds the genome directory
            # the name of that directory is used by snpEff
            datadir = os.path.join(args.workdir, "data")
            snpEffGenome = os.listdir(datadir)[0]
            snpEffConfig = makeSnpEffConfig(
                args.workdir, snpEffGenome, datadir)
        else:
            snpEffGenome = None
            snpEffConfig = None

        rfOuts = []
        if args.procs == 1:
            for chrom in chromDict:
                logFile = os.path.join(args.workdir, "log." + chrom)
                cmd, outfile = radiaFilter(
                    filterDirs, snpEffGenome, snpEffConfig, args, chrom, chromDict[chrom], tempDir, logFile)
                # the output is generated by the filter, not on stdout
                if execute(cmd):
                    raise Exception("RadiaFilter Call failed")
                if not correctLineCount(chromLines[chrom], outfile):
                    errmsg = "RadiaFilter sanity check failed on chrom %s\n" % (
                        chrom)
                    raise Exception(errmsg)
                with open(logFile, 'r') as f:
                    print >>sys.stderr, f.read()
        else:
            cmds = []
            rawOuts = dict()
            for chrom in chromDict:
                logFile = os.path.join(args.workdir, "log." + chrom)
                cmd, outfile = radiaFilter(
                    filterDirs, snpEffGenome, snpEffConfig, args, chrom, chromDict[chrom], tempDir, logFile)
                cmds.append(cmd)
                rawOuts[chrom] = outfile
            p = Pool(args.procs)
            values = p.map(execute, cmds, 1)
            # check if all output files are the same size as inputs
            # and print logging info to stderr
            for chrom in chromDict:
                logFile = os.path.join(args.workdir, "log." + chrom)
                with open(logFile, 'r') as f:
                    print >>sys.stderr, f.read()
                if not correctLineCount(chromLines[chrom], rawOuts[chrom]):
                    errmsg = "RadiaFilter sanity check failed on chrom %s\n" % (
                        chrom)
                    raise Exception(errmsg)

        # the radiaMerge command only uses the output directory and patient
        # name
        cmd, mergeOut = radiaMerge(args, tempDir)
        if execute(cmd):
            raise Exception("RadiaMerge Call failed")

        # tcga compliance (note that this does NOT create the vcfLog tag)
        if args.makeTCGAcompliant:
            format_vcf(
                mergeOut,
                args.outputFilename,
                args.filterRejects,
                args.filterGermline)
        else:
            shutil.move(mergeOut, args.outputFilename)

    finally:
        args.no_clean = True
        if not args.no_clean and os.path.exists(tempDir):
            shutil.rmtree(tempDir)

if __name__ == "__main__":
    __main__()
