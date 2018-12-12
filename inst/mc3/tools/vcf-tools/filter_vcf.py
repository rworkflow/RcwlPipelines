#!/usr/bin/env python

import argparse
import vcf

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tumor", default="TUMOR")
    parser.add_argument("--normal", default="NORMAL")
    parser.add_argument("--no-ad", action="store_true")
    parser.add_argument("--cutoff", type=int, default=None)
    parser.add_argument("vcf")
    parser.add_argument("out")

    args = parser.parse_args()

    vcf_reader = vcf.Reader(open(args.vcf, 'r'))
    vcf_writer = vcf.Writer(open(args.out, 'w'), vcf_reader)

    for record in vcf_reader:
        keep = True
        if args.cutoff is not None:
            for call in record.samples:
                if call.sample == args.tumor and args.no_ad == False:
                    for n, d in call.data._asdict().items():
                        if n == "AD" and int(d[1]) < args.cutoff:
                            keep = False
                elif call.sample == args.tumor and args.no_ad == True:
                    for n, d in record.INFO.items():
                        if n == "T_DP" and int(d) < args.cutoff:
                            keep = False
        if keep:
            vcf_writer.write_record(record)
