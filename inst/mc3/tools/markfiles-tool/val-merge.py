##
# This file merges filter marks into a MAF file.
#
# A mark key is generated for each record in the MAF file.
# The input MAF is then sorted by that mark key.
# Simultaneously input filter files are sorted by their mark keys.
#
# The mark keys are Tumor_Sample_Barcode|Chromosome|Start_Position|End_Position|Reference_Allele|Tumor_Seq_Allele2.
# 
# Sorted files are then streamed into the merger subroutines.  As they are sorted, 
# we do not need to hold all of the data in memory at one time.
#
# Input keys are read in batches such that all records with the same key are grouped.
# The marks for these keys are extracted, split by comma, merged, sorted and joined by comma.
# If the length of the marks for a variant key is zero, PASS is reported.

import os, os.path, sys
import csv, json

def mafkeyfun(record):
    return '|'.join([record['Tumor_Sample_Barcode'][:12], record['Chromosome'], 
                     record['Start_Position'], record['End_Position'], 
                     record['Reference_Allele'], record['Tumor_Seq_Allele2']])

def batch(iter):
    lastkey = None
    _batch = []
    mafrecord = []
    for r in iter:
        try:
            thiskey = r[0]
            thistype = r[1]
            thisval = {k:v for k,v in [vv.split('|',1) for vv in r[2:]]}
        except Exception as inst:
            log("Error parsing value: %s" % str(inst))
            log(str(r))
            continue
        if lastkey and thiskey != lastkey:
            # print 'yield batch %s' % lastkey
            yield mafrecord, _batch
            _batch = []
            mafrecord = []
        lastkey = thiskey
        if thistype == 'mafr':
            mafrecord.append(thisval)
        elif thistype == 'valr' and thisval['Validation_Status'] == 'Valid':
            _batch.append(thisval)
    yield mafrecord, _batch


def mainreduce(args):
    reader = csv.reader(args.INPUT, delimiter = '\t')
    # use the maf arg to get the fieldnames
    mafreader = csv.DictReader(args.maf, delimiter = '\t')
    fields = mafreader.fieldnames
    if 'MERGESOURCE' not in fields:
        fields.append('MERGESOURCE')
    writer = csv.DictWriter(args.output, fieldnames = fields, delimiter = '\t', extrasaction = 'ignore', restval = '.')
    writer.writeheader()
    
    # we read in a single stream for this application
    for mrs, vset in batch(reader):
        if not mrs:
            # This happens if there are no primary maf records for this validation
            # in this case we write out the aggergated record
            allkeys = set([k for v in vset for k in v.keys()])
            mr = {k: ','.join(list(set([vs[k] for vs in vset]))) for k in allkeys}
            mr['MERGESOURCE'] = 'COPY'
            writer.writerow(mr)
        else:
            for mr in mrs:
                if vset:
                    mr.update({'Verification_Status': ','.join(list(set([vs['Verification_Status'] for vs in vset]))),
                        'Validation_Status': ','.join(list(set([vs['Validation_Status'] for vs in vset]))),
                        'Mutation_Status': ','.join(list(set([vs['Mutation_Status'] for vs in vset]))),
                        'Validation_Method': ','.join(list(set([vs['Validation_Method'] for vs in vset])))})
                    mr['MERGESOURCE'] = 'MERGE'
                else:
                    mr['MERGESOURCE'] = 'PRIMARY'
                writer.writerow(mr)

def log(m):
    sys.stderr.write(m + '\n')

def mainmap(args):
    if not args.maf:
        log("--maf is a required argument")
        sys.exit(2)
    if not args.output:
        log('--output is a required argument')
        sys.exit(2)
    if not args.type in (0,1):
        log('--type must be either 0 or 1')
        sys.exit(2)
    writer = csv.writer(args.output, delimiter = '\t')
    
    log("Starting mark writing")
    for markfile in args.MARKFILES:
        with open(markfile, 'r') as fi:
            reader = csv.DictReader(fi, delimiter = '\t')
            for i, data in enumerate(reader):
                if data['Validation_Status'] == 'Valid':
                    writer.writerow([mafkeyfun(data), 'valr'] + ['%s|%s' % (k,v) for k, v in data.items()])
        log("wrote %s" % markfile)
    mafreader = csv.DictReader(args.maf, delimiter = '\t')
    for i, record in enumerate(mafreader):
        writer.writerow([mafkeyfun(record), 'mafr'] + ['%s|%s' % (k, v) for k, v in record.items()])
        if i % 100000 == 0:
            log("processed %s maf records" % str(i))
    log("Done with mapping")

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parsers = parser.add_subparsers()
    map_parser = parsers.add_parser('map', help = 'run the map function')
    map_parser.add_argument('--maf', type = argparse.FileType('r'), help = 'input maf file')
    map_parser.add_argument('--type', type = int, default = 0, help = 'should the normal be used in the key, 0 = no, 1 = yes')
    map_parser.add_argument('--output', type = argparse.FileType('w'), default = sys.stdout, help = 'output maf file')
    map_parser.add_argument('MARKFILES', nargs = '+', help = 'input mark files')
    map_parser.set_defaults(func=mainmap)
    
    reduce_parser = parsers.add_parser('reduce', help = 'run the reduce function')
    reduce_parser.add_argument('--maf', type = argparse.FileType('r'), help = 'maf file for headers')
    reduce_parser.add_argument('INPUT', type = argparse.FileType('r'), nargs = '?', default = sys.stdin, help = 'mapped input')
    reduce_parser.add_argument('--output', type = argparse.FileType('w'), help = 'output file')
    reduce_parser.add_argument('--append', action = 'store_true', help = 'should the filters be appended?')
    reduce_parser.set_defaults(func=mainreduce)

    args = parser.parse_args()
    
    args.func(args)
    
