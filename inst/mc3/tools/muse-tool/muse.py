#!/usr/bin/env python

import sys
import re
import os
import string
import shutil
import logging
import subprocess
import tempfile
from multiprocessing import Pool
from argparse import ArgumentParser

def which(cmd):
    cmd = ["which",cmd]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    res = p.stdout.readline().rstrip()
    if len(res) == 0: return None
    return res

def fai_chunk(path, blocksize):
    seq_map = {}
    with open( path ) as handle:
        for line in handle:
            tmp = line.split("\t")
            seq_map[tmp[0]] = long(tmp[1])

    for seq in seq_map:
        l = seq_map[seq]
        for i in xrange(1, l, blocksize):
            yield (seq, i, min(i+blocksize-1, l))

def cmd_caller(cmd):
    logging.info("RUNNING: %s" % (cmd))
    print "running", cmd
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if len(stderr):
        print stderr
    return p.returncode

def cmds_runner(cmds, cpus):
    p = Pool(cpus)
    values = p.map(cmd_caller, cmds, 1)
    return values

def call_cmd_iter(muse, ref_seq, block_size, tumor_bam, normal_bam, contamination, output_base):
    contamination_str = ""
    if contamination is not None:
        contamination_str = "-p %s" % (contamination)
    template = string.Template("${MUSE} call -f ${REF_SEQ} ${CONTAMINATION} -r '${INTERVAL}' ${TUMOR_BAM} ${NORMAL_BAM} -O ${OUTPUT_BASE}.${BLOCK_NUM}")
    for i, block in enumerate(fai_chunk( ref_seq + ".fai", block_size ) ):
            cmd = template.substitute(
                dict(
                    REF_SEQ=ref_seq,
                    CONTAMINATION=contamination_str,
                    BLOCK_NUM=i,
                    INTERVAL="%s:%s-%s" % (block[0], block[1], block[2]) ),
                    MUSE=muse,
                    TUMOR_BAM=tumor_bam,
                    NORMAL_BAM=normal_bam,
                    OUTPUT_BASE=output_base
            )
            yield cmd, "%s.%s.MuSE.txt" % (output_base, i)

def run_muse(args):

    mode_flag = ""
    if args.muse.endswith("MuSEv1.0rc"):
        args.p = None
        if args.mode == "wgs":
            mode_flag = "-G"
        else:
            mode_flag = "-E"

    if not os.path.exists(args.muse):
        args.muse = which(args.muse)

    workdir = os.path.abspath(tempfile.mkdtemp(dir=args.workdir, prefix="muse_work_"))

    if not os.path.exists(args.f + ".fai"):
        new_ref = os.path.join(workdir, "ref_genome.fasta")
        os.symlink(os.path.abspath(args.f),new_ref)
        subprocess.check_call( ["/usr/bin/samtools", "faidx", new_ref] )
        args.f = new_ref

    if args.normal_bam_index is None:
        if not os.path.exists(args.normal_bam + ".bai"):
            new_bam = os.path.join(os.path.abspath(workdir), "normal.bam")
            os.symlink(os.path.abspath(args.normal_bam),new_bam)
            subprocess.check_call( ["/usr/bin/samtools", "index", new_bam] )
            args.normal_bam = new_bam
    else:
        new_bam = os.path.join(os.path.abspath(workdir), "normal.bam")
        os.symlink(os.path.abspath(args.normal_bam), new_bam)
        os.symlink(os.path.abspath(args.normal_bam_index), new_bam + ".bai")
        args.normal_bam = new_bam

    if args.tumor_bam_index is None:
        if not os.path.exists(args.tumor_bam + ".bai"):
            new_bam = os.path.join(os.path.abspath(workdir), "tumor.bam")
            os.symlink(os.path.abspath(args.tumor_bam),new_bam)
            subprocess.check_call( ["/usr/bin/samtools", "index", new_bam] )
            args.tumor_bam = new_bam
    else:
        new_bam = os.path.join(workdir, "tumor.bam")
        os.symlink(os.path.abspath(args.tumor_bam), new_bam)
        os.symlink(os.path.abspath(args.tumor_bam_index), new_bam + ".bai")
        args.tumor_bam = new_bam

    cmds = list(call_cmd_iter(ref_seq=args.f,
        muse=args.muse,
        block_size=args.b,
        tumor_bam=args.tumor_bam,
        normal_bam=args.normal_bam,
        contamination=args.p,
        output_base=os.path.join(workdir, "output.file"))
    )

    rvals = cmds_runner(list(a[0] for a in cmds), args.cpus)
    if any(rvals):
        raise Exception("MuSE CALL failed")
    #check if rvals is ok
    first = True
    merge = os.path.join(workdir, "merge.output")
    with open(merge, "w") as ohandle:
        for cmd, out in cmds:
            with open(out) as handle:
                for line in handle:
                    if first or not line.startswith("#"):
                        ohandle.write(line)
            first = False
            if not args.no_clean:
                os.unlink(out)

    dbsnp_file = None
    if args.D:
        new_dbsnp = os.path.join(workdir, "db_snp.vcf.gz")
        os.symlink(args.D,new_dbsnp)
        #subprocess.check_call( ["/usr/bin/bgzip", new_dbsnp] )
        subprocess.check_call( ["/usr/bin/tabix", "-p", "vcf", new_dbsnp ])
        dbsnp_file = new_dbsnp
        sump_template = string.Template("${MUSE} sump -I ${MERGE} -O ${OUTPUT} -D ${DBSNP} ${MODE}")
    else:
        sump_template = string.Template("${MUSE} sump -I ${MERGE} -O ${OUTPUT} ${MODE}")

    tmp_out = os.path.join(workdir, "tmp.vcf")
    sump_cmd = sump_template.substitute( dict (
        MUSE=args.muse,
        MERGE=merge,
        OUTPUT=tmp_out,
        DBSNP=dbsnp_file,
        MODE=mode_flag
    ))
    cmd_caller(sump_cmd)

    if args.muse.endswith("MuSEv0.9.9.5"):
        subprocess.check_call( ["/opt/bin/vcf_reformat.py", tmp_out, "-o", args.O,
            "-b", "TUMOR", args.tumor_bam, "-b", "NORMAL", args.normal_bam] )
    else:
        shutil.copy(tmp_out, args.O)

    if not args.no_clean:
        shutil.rmtree(workdir)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-m", "--muse", help="Which Copy of MuSE", choices=["MuSEv0.9.9.5", "MuSEv1.0rc"], default="MuSEv0.9.9.5")
    parser.add_argument("-f", help="faidx indexed reference sequence file", required=True)
    #parser.add_argument("-r", help="single region (chr:pos-pos) where somatic mutations are called")
    #parser.add_argument("-l", help="list of regions (chr:pos-pos or BED), one region per line")
    parser.add_argument("-p", type=float, help="normal data contamination rate [0.050]", default=0.05)
    parser.add_argument("-b", type=long, help="Parallel Block Size", default=50000000)
    parser.add_argument("-O", help="output file name (VCF)", default="out.vcf")
    parser.add_argument("-D", help="""dbSNP vcf file that should be bgzip compressed,
tabix indexed and based on the same reference
genome used in 'MuSE call'""")

    parser.add_argument("-n", "--cpus", type=int, default=8)
    parser.add_argument("-w", "--workdir", default="/tmp")
    parser.add_argument("--no-clean", action="store_true", default=False)
    parser.add_argument("--mode", choices=["wgs", "wxs"], default="wgs")
    parser.add_argument("--tumor-bam", dest="tumor_bam", required=True)
    parser.add_argument("--tumor-bam-index", dest="tumor_bam_index", default=None)
    parser.add_argument("--normal-bam", dest="normal_bam", required=True)
    parser.add_argument("--normal-bam-index", dest="normal_bam_index", default=None)
    args = parser.parse_args()
    run_muse(args)
