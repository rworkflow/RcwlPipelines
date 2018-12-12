import argparse
import os, os.path, sys
import hgsc_vcf

SORTED_TIERS = ['Tier1', 'Tier2', 'Tier3', 'Tier4', 'Tier5']

def convert_record(record, convert, filter):
    # converts in place
    for con in convert:
        if con in record['FILTER']:
            # convert
            record['INFO']['OF'] = ['|'.join(record['FILTER'])]
            record['FILTER'] = ['PASS']


parser = argparse.ArgumentParser()

parser.add_argument(
    '--level',
    default = '5',
    type = str,
    help = 'one of 1,2,3,4,5,all.  level represents the last level that will still be filtered, all implies that all levels will be retained (their filter columns will be converted to PASS).')

parser.add_argument('INPUT', type = argparse.FileType('r'), help = 'input file')
parser.add_argument('OUTPUT', type = argparse.FileType('w'), help = 'output file')

args = parser.parse_args()

if args.level == 'all':
    convert = SORTED_TIERS
    filter = []
    
elif args.level in ['1','2','3','4','5']:
    cut = int(args.level) - 1
    convert = SORTED_TIERS[:cut]
    filter = SORTED_TIERS[cut:]

else:
    raise ValueError("%s is not a valid level" % args.level)

reader = hgsc_vcf.Reader(args.INPUT)
header = reader.header

header.add_header('##INFO=<ID=OF,Number=1,Type=String,Description="original tiering call for this variant in this sample">')
header.add_header('##COMMAND=<ID=filter_muse.py,Params="%s">' % ' '.join(sys.argv))

writer = hgsc_vcf.Writer(args.OUTPUT, header)
writer.write_header()
for record in reader:
    convert_record(record, convert, filter)
    writer.write_record(record)
