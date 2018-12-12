#!/usr/bin/env python

import sys
import re
import logging
import pysam
from argparse import ArgumentParser

class MetaData:
    def __init__(self, line):
        res = re.search(r'^\#\#(\w+)=(.*)$', line)
        if res is None:
            raise IOException("Bad Metadata line")
        self.name = res.group(1)
        self.record = res.group(2)

    def __str__(self):
        return "##%s=%s" % (self.name, self.record)

class Record:
    def __init__(self, line):
        self.line = line.rstrip("\r\n")
        tmp = self.line.split("\t")
        self.seq = tmp[0]
        self.pos = long(tmp[1])
        self.id = tmp[2]
        self.ref = tmp[3]
        self.alt = tmp[4]
        self.qual = tmp[5]
        self.filter = tmp[6]
        self.info = tmp[7]
        self.format = ""
        self.samples = []

    def __str__(self):
        return "\t".join( [
            self.seq, "%d" % (self.pos), self.id, self.ref, self.alt,
            self.qual, self.filter, self.info, self.format] + self.samples )

def atoi(text):
    return int(text) if text.isdigit() else text.lower()

class VCF:

    def __init__(self):
        self.metadata = []
        self.records = []
        self.header = None
        self.samples = []

    def parse(self, handle):
        for line in handle:
            if line.startswith("##"):
                self.metadata.append( MetaData(line) )
            elif line.startswith("#CHROM"):
                self.header = line.rstrip("\n\r")
                self.samples = line.rstrip("\n\r").split("\t")[9:]
            else:
                self.records.append( Record(line) )

    def adjust_format(self, bam_files):
        bams = []
        for b in bam_files:
            bams.append(pysam.AlignmentFile(b, "rb"))

        self.records.sort(key=lambda x: [atoi(x.seq), x.pos])

        #loop through each record
        for rec in self.records:
            all_reads = []
            good_reads = []
            all_quals = []
            good_quals = []
            for sample in bams:
                #for each sample grab the pileup
                pileup = sample.pileup(rec.seq, rec.pos-1, rec.pos)
                for pile in pileup:
                    #we're only concered about one specific position
                    if pile.reference_pos == rec.pos-1:
                        all_sread = []
                        all_squal = []
                        good_sread = []
                        good_squal = []
                        #grab all the reads and quality scores that map to this postion
                        for row in pile.pileups:
                            if not row.indel and not row.is_del:
                                good_read = True
                                c_read = row.alignment.query_sequence[row.query_position]
                                c_qual = row.alignment.query_qualities[row.query_position]
                                if row.alignment.mapping_quality <= 0:
                                    good_read = False
                                if c_qual < 5:
                                    good_read = False
                                #remove reads that have been mapped by BWA 'mate rescue'
                                if row.alignment.has_tag('XT') and row.alignment.get_tag('XT') == 'M':
                                    good_read = False
                                soft_clip_size = sum( a[1] for a in row.alignment.cigar if a[0] == 4 )
                                if float(soft_clip_size) / float(row.alignment.query_length) >= 0.3:
                                    good_read = False
                                if good_read:
                                    good_sread.append( c_read )
                                    good_squal.append( c_qual )
                                all_sread.append( c_read )
                                all_squal.append( c_qual )
                        all_reads.append(all_sread)
                        all_quals.append(all_squal)
                        good_reads.append(good_sread)
                        good_quals.append(good_squal)
            sample_fields = []
            #for each of the samples, count the reads and average the quailiy scores for the alt reads
            for sample_all_reads, sample_all_quals, sample_good_reads, sample_good_quals in zip(all_reads, all_quals, good_reads, good_quals):
                out = {}
                #create True/False arrays for the reads
                refs = map(lambda x: x == rec.ref, sample_all_reads)
                alts = map(lambda x: x == rec.alt, sample_all_reads)
                good_alts = map(lambda x: x == rec.alt, sample_good_reads)
                good_alt_quals = list( q * a for q,a in zip(sample_good_quals, good_alts) )
                if sum(good_alts):
                    out['BQ'] = sum(good_alt_quals) / float(sum(good_alts))
                else:
                    out['BQ'] = 0
                out['AD'] = sum(good_alts)
                out['DP'] = len(sample_all_reads)
                out['DP4'] = len(sample_good_reads)
                sample_fields.append(out)

            #FIXME: This may also need the SS field to be TCGA compliant
            rec.format = "DP:DP4:AD:BQ" #"GT:DP:AD:BQ"
            rec.samples = []
            for i, c in enumerate(sample_fields):
                info = "%d:%d:%d:%d" % ( c['DP'], c['DP4'], c['AD'], c['BQ'] )
                rec.samples.append(info)

        for b in bams:
            b.close()

        #FIXME: This assumes that there are not existing FORMAT fields
        self.metadata.append( MetaData('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">') )
        self.metadata.append( MetaData('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth at this position in the sample">') )
        self.metadata.append( MetaData('##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Depth of reads supporting alleles 0/1/2/3...">') )
        self.metadata.append( MetaData('##FORMAT=<ID=BQ,Number=.,Type=Integer,Description="Average base quality for reads supporting alleles">') )


    def write(self, handle):
        for meta in self.metadata:
            handle.write("%s\n" % (meta) )
        handle.write("%s\n" % ( "\t".join(self.header.split("\t")[:8] + ['FORMAT'] + self.samples ) ) )
        for record in self.records:
            handle.write("%s\n" % (record) )

def run_adjust(vcf, samples, out):

    input = VCF()
    with open(vcf) as handle:
        input.parse(handle)

    if len(samples) != len(input.samples):
        logging.info("Sample MisMatch")
        if len(input.samples ) < len(samples):
            input.samples += [None] * (len(samples) - len(input.samples))

    s = []
    for b,c in zip(list(a[0] for a in samples), input.samples):
        if c is not None:
            s.append(c)
        else:
            s.append(b)
    input.samples = s
    input.adjust_format( list(a[1] for a in samples) )
    if out is None:
        input.write(sys.stdout)
    else:
        with open(out, "w") as handle:
            input.write(handle)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("vcf")
    parser.add_argument("-b", dest="samples", nargs=2, action="append")
    parser.add_argument("-o", "--out", default=None)

    args = parser.parse_args()
    v = vars(args)
    run_adjust(**v)
