import collections
import itertools
from collections import OrderedDict
from hgsc_vcf.io import *
import logging

logger = logging.getLogger('hgsc_vcf')
logger.addHandler(logging.NullHandler())
logger.setLevel(logging.WARN)
## Exceptions


##
# merge records
#
# merges a set of records, specifically, making a new record with 
# the CHROM, POS, REF, and ALT data, FILTER data must be specified.
# Samples are organized based on the samples_helper, which is 
# either a samples index (dict of key: name value: index or 
# Reader, whereby the Reader._samples_index is checked for the same.
# Old filter fields are added as FL:
#  ##FORMAT=<ID=FL,Number=1,Type=String,Description="Filter field prior to merger">
def merge_records(records, samples_helper, filter):
    ##
    # records are all checked to ensure that they have the same 
    # CHROM, POS, REF, and ALT.  If this is not the case then 
    # None is returned.
    if not (checkEqualIvo([r.CHROM for r in records]) and 
            checkEqualIvo([r.POS for r in records]) and 
            checkEqualIvo([r.REF for r in records]) and
            checkEqualIvo([r.ALT for r in records if r.ALT is not None]) and 
            checkEqualIvo([r.FORMAT for r in records])):
        return None
    chrom = records[0].CHROM
    pos = records[0].POS
    id = records[0].ID
    ref = records[0].REF 
    alt = records[0].ALT 
    format = records[0].FORMAT
    quals = [r.QUAL for r in records if r.QUAL is not None] 
    if len(quals) < 1:
        qual = None
    else:
        qual = max(quals)
    if isinstance(samples_helper, Reader) or isinstance(samples_helper, Writer):
        sample_indexes = samples_helper._sample_indexes
    else:
        sample_indexes = samples_helper
    info = {}
    allsamples = []
    for r in records:
        for s in r.samples:
            append_format(s, 'FL', '|'.join(r.FILTER), 'String', 1)
    allsamples = {sample.sample:sample for record in records for sample in record.samples}
    samples = [allsamples[s] for s, i in sorted(sample_indexes.items(), key = lambda x: x[1])]
    
    return _Record(chrom, pos, id, ref, alt, qual, filter, info, format,
            sample_indexes, samples)

##
# check that all values in a list are equal
# http://stackoverflow.com/q/3844948/
def checkEqualIvo(lst):
    return not lst or lst.count(lst[0]) == len(lst)

##
# check that the GT field is identical for all samples in the record
#
# used generally for working with Wheeljack output.
# @param vcf.Record
# @returns logical 
def check_all_gt(samples):
    return checkEqualIvo([sample['GT'] for sample in samples.values()])

##
# split and return the gt info for a record
#
# @param vcf.Record
# @returns list of gt inicies
def split_gt(record):
    samples = record.get('SAMPLES', {})
    if not check_all_gt(samples):
        raise ValueError("Record does not contain all the same GT info")
    gt = [int(i) for i in samples.values()[0]['GT'][0].split('/')]
    return gt

##
# return indexes of ref and alt sorted by the sum of ac in the records
#
# @param vcf.Record
# @returns list of indicies
def sum_ac(record):
    # if not check_all_gt(record):
    #     raise ValueError("Record does not contain all the same GT info")
    samples = record.get('SAMPLES', {})
    ac_sum = [0]*len(samples.values()[0]['AC'])
    for sample in samples.values():
        for i, val in enumerate(sample['AC']):
            ac_sum[i] += int(val)
    return ac_sum

##                                                                                                                                 
# select the "best" alt for use as an alt                                                                                          
#                                                                                                                                  
# reference is always selected (0).  This utility function returns the index of the highest alt                                    
def best_alt_index(record, gt = None):
    if gt is None:
        gt = split_gt(record)                                                                                                 
    if len(gt) < 2:                                                                                                                
        return -1 # special return value indicating that this is a reference only location (homozygous reference or no coverage)   
    if len(gt) == 2:
        return 0 if gt[0] != 0 else 1
    acsum = sum_ac(record)                                                                                                
    acorder = [e[0] for e in sorted(enumerate(acsum), key = lambda x: x[1], reverse = True)]
    for i in acorder:                                           
        if gt[i] != 0:                                                                                                             
            return i                                                                                                               
    raise ValueError("Theoretically unreachable code, length of GT is greater than 1 but all gt seem to point to 0")               

def ref_index(record, gt = None):
    if gt is None:
        gt = split_gt(record)
    if 0 not in gt:
        # this record is not like wheeljack records, the index is therefore 0
        return 0
    for i, g in enumerate(gt):
        if g == 0:
            return i
    raise ValueError("Theoretically unreachable code, no gt == 0")

##
# filter vcf
#
# yields lines from a vcf_handle after filtering
# note that no other function should execute the next function of the reader
def filter_vcf(reader, filter_function):
    for record in reader:
        if filter_function(record):
            yield record
       
       
## 
# process a vcf by applying the process_funciton to all records
#
# an optional filter_function can also be supplied for read record filtering
def process_vcf(reader, process_function, filter_function = None):
    if filter_function is None:
        def filter_function(record):
            return True
    for record in filter_vcf(reader, filter_function):
        yield process_function(record)

##
# simplifies an allele
#
# This strips off leading and lagging identical bases.
# Stripping takes place from the right end first and then the left.
# The position is modified as needed if there is any stripping from the left.
def _simplify_allele(ref, alt, pos):
    assert isinstance(alt, list), "Alt must be in hgsc_vcf alt format, which is a list"
    rc = []
    lc = []
    for l, r in [_get_slice_indicies(ref, a) for a in alt]:
        rc.append(r)
        lc.append(l)
    l = min(lc)
    r = min(rc)
    return ref[l:len(ref)-r], [a[l:len(a)-r] for a in alt], pos + l

def _get_slice_indicies(ref, alt):
    alt = str(alt)
    ref = str(ref)
    # strip the right
    r_len = len(ref)
    a_len = len(alt)
    minlen = min(r_len, a_len)
    if minlen == 1:
        logger.debug("%s and %s are of length 1", ref, alt)
        return 0, 0
    l = r = 0
    for i in range(minlen):
        r = i
        if ref[r_len - 1 - i] != alt[a_len - 1 - i]:
            break
    ref = ref[:(r_len - r)]
    alt = alt[:(a_len - r)]
    logger.debug("clipped at %s to get %s and %s", r, ref, alt)
    # strip the left
    r_len = len(ref)
    a_len = len(alt)
    minlen = min(a_len, r_len)
    if minlen == 1:
        return l, r
    for i in range(minlen):
        l = i
        if ref[i] != alt[i]:
            break
    logger.debug("final clip is %s and %s", l, r)
    return l, r

##
# select allele
#
# selects the "best" alt from a set of multi-alt calls (like from Wheeljack)
# This is actually a convenience function, because the selection is based 
# on the result of a selection_function parameter provided to the function.
# The "best" alt (index) is then simplified and returned.  The result of
# this function is a list allowing multiple alts to be determined to be the
# "best" in case there is a tie (or near tie, or whatever).  The selection_function must take
# the record as it's input.
def select_allele(record, selection_function, simplify = False):
    if len(record['ALT']) < 2:
        return [record]
    gt = split_gt(record)
    ref_i = ref_index(record, gt) 
    newrecords = []
    for alt_index in selection_function(record):
        alt = record['ALT'][gt[alt_index] - 1] # since ALT is a list but the gt index is off by one.
        # make a new record
        if simplify: # it may or may not be a good idea to simplify
            ref, alt, pos = _simplify_allele(record['REF'], [alt], record['POS'])
        else:
            ref = record['REF']
            pos = record['POS']
            alt = [alt]
        new_record = OrderedDict()
        for k, v in record.items():
            if k == 'ALT':
                new_record[k] = alt
            elif k == 'REF':
                new_record[k] = ref
            elif k == 'POS':
                new_record[k] = str(pos)
            elif k == 'SAMPLES':
                # need to split the samples and update their info, sample info will always be in ref/alt order
                new_samples = OrderedDict()
                for s, sinfo in v.items():
                    new_sinfo = OrderedDict()
                    for fk, finfo in sinfo.items():
                        if fk == 'GT':
                            # since ref/alt order the GT should be 0/1 in our strange cancer convention
                            new_sinfo[fk] = ['0/1']
                        elif len(finfo) == len(gt):
                            # ensure that the data are in ref/alt order
                            new_sinfo[fk] = [finfo[ref_i], finfo[alt_index]]
                        else:
                            new_sinfo[fk] = finfo
                    new_samples[s] = new_sinfo
                new_record[k] = new_samples
            else:
                new_record[k] = v
        newrecords.append(new_record)
    ##
    # records are processed in batches and sorted to maintain sorted order of the VCF file
    newrecords = sorted(newrecords, key = lambda x: x['POS'])
    return newrecords

# extract allele fraction and remainder informaiton for sliced alleles from a record
