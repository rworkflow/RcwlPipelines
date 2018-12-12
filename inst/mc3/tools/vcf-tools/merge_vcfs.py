
import os, os.path, sys
import glob
import hgsc_vcf
import logging
from collections import *

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

##
# Goal; to merge multiple vcf files for a subject.
# This will merge files first buffering calls together and then resolving the call.
#
# we will need:
#   1. object to read the files and perform chunking
#   2. routines to handle the merger of the events (subroutine of the class or independent)
#   3. writer (hgsc_vcf.Writer class)
#   4. something to merge the headers?

class MetaRecord(object):
    def __init__(self, caller, record):
        self.caller = caller
        self.record = record
    
    def _seq_pos(self, record):
        if record['CHROM'] == 'X':
            return 30 * 1000000000 + record['POS']
        elif record['CHROM'] == 'Y':
            return 40 * 1000000000 + record['POS']
        elif record['CHROM'] in ('M', 'MT'):
            return 50 * 1000000000 + record['POS']
        else:
            return int(record['CHROM']) * 1000000000 + record['POS']

    def __cmp__(self, other):
        si = self._seq_pos(self.record)
        oi = self._seq_pos(other.record)
        return si - oi
        '''
        block commented out because we don't actually care that the alts aren't the same???
        if si != oi:
            return si - oi
        else:
            return sum([hash(a) for a in self.record['ALT']]) - sum([hash(a) for a in other.record['ALT']])
        '''
    def __repr__(self):
        return "<record chr=%s, pos=%s, alt=%s>" % (self.record['CHROM'], self.record['POS'], self.record['ALT'])

class MetaReader(object):
    def __init__(self, fobj):
        self.reader = hgsc_vcf.Reader(fobj)
        self.caller = fobj.name
        # get the normal and primary sample ids
        sampleMapping = {l.fields.get('ID'):l.fields.get('SampleTCGABarcode') for l in self.reader.header.get_headers('SAMPLE')}
        if 'PRIMARY' not in sampleMapping and 'METASTATIC' in sampleMapping:
            sampleMapping['PRIMARY'] = sampleMapping['METASTATIC']
        elif 'PRIMARY' not in sampleMapping and 'RECURRANCE' in sampleMapping:
            sampleMapping['PRIMARY'] = sampleMapping['RECURRANCE']
        logger.info("Sample mapping for %s: %s", fobj.name, sampleMapping)
        self.normal = sampleMapping.get('NORMAL', "Unknown normal sample ID")
        self.primary = sampleMapping.get('PRIMARY', "Unknown primary sample ID")
        self._next = None
        self.take() # call to take this time will return None but will also fast forward the reader to the next position
    def __cmp__(self, other):
        return self._next.__cmp__(other._next)

    def peek(self):
        return self._next
    
    def take(self):
        old = self._next
        # fast forward to the next
        new = None
        while True:
            try:
                n = self.reader.next()
                if 'GL' in n['CHROM']:
                    logger.info("GL in chrom %s", n['CHROM'])
                    logger.info("Closing %s", self.caller)
                    #self.reader.fobj.close()
                    new = None
                else:
                    if 'NORMAL' not in n['SAMPLES']:
                        n['SAMPLES']['NORMAL'] = n['SAMPLES'][self.normal]

                    if 'PRIMARY' not in n['SAMPLES']: 
                        if self.primary in n['SAMPLES']:
                            n['SAMPLES']['PRIMARY'] = n['SAMPLES'][self.primary]
                        elif 'METASTATIC' in n['SAMPLES']:
                            n['SAMPLES']['PRIMARY'] = n['SAMPLES']['METASTATIC']
                        elif 'RECURRANCE' in n['SAMPLES']:
                            n['SAMPLES']['PRIMARY'] = n['SAMPLES']['RECURRANCE']
                        else:
                            raise ValueError("Can't find the PRIMARY sample")

                    if n['SAMPLES']['NORMAL']['GT'][0] not in  ('0/0', '0', '.', './.'):
                        continue

                    if 'PASS' not in n['FILTER']:
                        continue

                    new = n
                    break

            except StopIteration: # swallow the error and just set to None
                logger.info("Stopped iteration")
                logger.info("Closing %s", self.caller)
                self.reader.fobj.close()
                new = None
                break
        if new is None:
            self._next = None
        else:
            self._next = MetaRecord(self.caller, new)
        return old

    def __repr__(self):
        if self._next is None:
            return "<%s, closed>" % self.caller
        else:
            return "<%s, chr=%s, pos=%s>" % (self.caller, self._next.record['CHROM'], self._next.record['POS'])

    

class MultiVCFReader(object):
    def __init__(self, infiles, outfile, keys):
        self.buffer = 10
        self.infiles = {f:MetaReader(open(f, 'r')) for f in infiles}
        self.outfile = outfile
        self.outwriter = hgsc_vcf.Writer(open(self.outfile, 'w'), self.generate_header())

        # get the normal and primary sample ids
        sampleMapping = {l.fields.get('ID'):l.fields.get('SampleTCGABarcode') for l in self.infiles.values()[0].reader.header.get_headers('SAMPLE')}

        if 'PRIMARY' not in sampleMapping and 'METASTATIC' in sampleMapping:
            sampleMapping['PRIMARY'] = sampleMapping['METASTATIC']
        elif 'PRIMARY' not in sampleMapping and 'RECURRANCE' in sampleMapping:
            sampleMapping['PRIMARY'] = sampleMapping['RECURRANCE']

        self.normal = sampleMapping.get('NORMAL', "Unknown normal sample ID")
        self.primary = sampleMapping.get('PRIMARY', "Unknown primary sample ID")
        self.keymap = dict(zip(infiles, keys))

    # lets make this a generator so that we can keep up with the sorting
    # very likely that once something is in sorted order we can keep it that way
    def get_next(self):
        while True:
            nextsort = sorted([mr for mr in self.infiles.values() if mr.peek() is not None])
            if len(nextsort) < 1:
                raise StopIteration() # we can break here
            _raises = True
            for r in self._get_sorted_next(nextsort):
                _raises = False
                yield r
            if _raises:
                logger.info("Sorted %s", nextsort)
                if len(nextsort) > 1:
                    logger.info("n0 <= n1: %s", nextsort[0] <= nextsort[1])
                    logger.info("nextsort[0]._next: %s", nextsort[0]._next)
                raise ValueError("Entered sorted but did not yield")
    # helper generator 
    def _get_sorted_next(self, nextsort):
        if len(nextsort) < 2:
            n0 = nextsort[0]
            while n0.peek() is not None:
                yield n0.take()
        else:
            n0 = nextsort[0]
            n1 = nextsort[1]
            while n0.peek() is not None and n0 <= n1:
                yield n0.take()
            if n0.peek() is None:
                logger.info("Reached the end of %s", n0.caller)


    def generate_header(self):
        newHeader = hgsc_vcf.VCFHeader()
        newHeader.samples = ['NORMAL', 'PRIMARY'] # deterministic sample names now
        for infile in self.infiles.values():
            newHeader.headers += infile.reader.header.headers # append all of the headers together, who cares, we can sort out later
        newHeader.add_header('##COMMAND=<ID=vcf-merge>')
        return newHeader

    ##
    # yields batches of MetaRecord's
    def chunk(self):
        batch = []
        cpos = 0
        for r in self.get_next(): # run through the generator
            if len(batch) < 1:
                batch.append(r)
                cpos = r.record['POS'] + max(len(r.record['REF']), max([len(a) for a in r.record['ALT']])) - 1 + self.buffer
            elif r.record['CHROM'] != batch[0].record['CHROM']:
                yield batch
                batch = [r]
            elif r.record['POS'] > (batch[-1].record['POS'] + self.buffer):
                yield batch
                batch = [r]
                cpos = r.record['POS'] + max(len(r.record['REF']), max([len(a) for a in r.record['ALT']])) - 1 + self.buffer
            else:
                batch.append(r)
                cpos = max(cpos, r.record['POS'] + max(len(r.record['REF']), max([len(a) for a in r.record['ALT']])) - 1 + self.buffer)
        yield batch # yield the last batch

    # make this an iterable
    def __iter__(self):
        return self.chunk()

def parseInfo(merge, type):
    sdp = sad = 0.0
    for m in merge:
        s = m.record['SAMPLES'][type]
        if 'BCOUNT' in s:
            aindex = ['A', 'C', 'G', 'T'].index(m.record['ALT'][0])
            sad += int(s['BCOUNT'][aindex])
        else:
            sad += int(s['AD'][-1])
        sdp += int(s['DP'][0])
    result = {
            'GT': ['0/1'] if type == 'PRIMARY' else ['0/0'],
            'DP': [str(int(sdp / len(merge)))],
            'AD': [str(int( (sdp - sad) / len(merge) )), str(int( sad / len(merge) ))]
            }
    return result

def resolve_merge(merge, callermap):
    newRecord = OrderedDict()
    for k in ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'):
        newRecord[k] = merge[0].record[k]
    newRecord['FORMAT'] = ['GT', 'DP', 'AD']
    newRecord['SAMPLES'] = OrderedDict()
    normalInfo = parseInfo(merge, 'NORMAL')
    tumorInfo = parseInfo(merge, 'PRIMARY')
    newRecord['SAMPLES']['NORMAL']  = OrderedDict([(k, normalInfo[k]) for k in newRecord['FORMAT']])
    newRecord['SAMPLES']['PRIMARY'] = OrderedDict([(k, tumorInfo[k]) for k in newRecord['FORMAT']])
    return newRecord, '|'.join([callermap[m.caller] for m in merge])

def contains_pindel(batch):
    for r in batch:
        if 'pindel' in r.caller:
            return True
    return False

def resolve_records(batch, callermap):
    if len(batch) == 1:
        yield resolve_merge(batch, callermap)
    else:
        logger.info("Processing batch size: %s", len(batch))
        if contains_pindel(batch):
            # it's a pindel call, merge all and yield the pindel call
            pc = [r for r in batch if 'pindel' in r.caller][0] # this must be true, but there might be more than one???
            callset = []
            for r in batch:
                if 'pindel' in r.caller:
                    callset.append(callermap[r.caller])
                else:
                    callset.append(callermap[r.caller] + '*')
            logger.info("Merged pindel call with %s", callset)
            yield resolve_merge([pc], callermap)[0], '|'.join(callset)
        else:
            merge = []
            for r in batch:
                if len(merge) < 1:
                    merge.append(r)
                else:
                    p = merge[0]
                    if r != p: # takes "advantage" of the cmp method, convenient for this purpose but otherwise dangerous
                        logger.info("Splitting batch because not the same")
                        logger.info(r)
                        logger.info(p)
                        logger.info("Yielding batch of %s", len(merge))
                        yield resolve_merge(merge, callermap)
                        merge = [r]
                    else:
                        merge.append(r)
            logger.info("Yielding batch of %s", len(merge))
            yield resolve_merge(merge, callermap)


def main(args):
    reader = MultiVCFReader(args.input, args.output, args.keys)
    reader.outwriter.header.add_header('##INFO=<ID=CENTERS,Number=1,Type=String,Description="Center files that made the call">')
    reader.outwriter.write_header()
    
    callermap = dict(zip(args.input, args.keys))

    for chunk in reader:
        for r, c in resolve_records(chunk, callermap):
            r['INFO'] = {'CENTERS':[c]}
            reader.outwriter.write_record(r)
    logger.info("Done")
            

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--keys', type = str, help = 'caller keys', nargs = '+')
    parser.add_argument('--output', type = str, help = 'output file')
    parser.add_argument('--input', nargs='+', type = str, help = 'input files')

    args = parser.parse_args()
    main(args)
