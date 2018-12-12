#!/bin/bash

usage()
{
        echo "reformat_vcf.sh"
        echo "Kami E. Chiotti 06.05.18"
        echo
        echo "Reformats IndelGenotyper VCF by converting and moving the INFO fields N_DP, N_SC, T_DP, and T_SC to"
        echo "AD, DP, ADF, and ADR FORMAT fields. If present, leaves SOMATIC flag in INFO."
        echo
        echo "Usage: $0 [-i input_vcf ] [-o output_vcf]"
        echo
        echo "  [-i input_vcf]            - Full path to Indelocator-generated VCF file for reformatting."
        echo "  [-o output_vcf]           - Output filename"
        echo "  [-h help]                 - Display this page"
        exit
}

input_vcf=()
output_vcf=""

while getopts ":i:o:h" Option
do
    case $Option in
        i ) input_vcf="$OPTARG" ;;
        o ) output_vcf="$OPTARG" ;;
        h ) usage ;;
        * ) echo "Unrecognized argument. Use '-h' for usage information."; exit 255 ;;
    esac
done
shift $(($OPTIND -1))

if [[ "$input_vcf" == "" ]]
then
    usage
fi

if [[ "$output_vcf" == "" ]]
then
    output_vcf="${input_vcf%.*}.reformatted.vcf"
fi

if [[ ! -r "$input_vcf" ]] 
then
    echo " Error: Can't open input VCF ($1)." >&2
    exit 1
fi

IFS=$'\n'
for LINE in $(cat $input_vcf)
do
    if [[ $LINE =~ ^# ]]
    then
        echo $LINE >> $output_vcf
    else
        CHROM=`echo $LINE | cut -f 1`
        POS=`echo $LINE | cut -f 2`
        ID=`echo $LINE | cut -f 3`
        REF=`echo $LINE | cut -f 4`
        ALT=`echo $LINE | cut -f 5`
        QUAL=`echo $LINE | cut -f 6`
        FILTER=`echo $LINE | cut -f 7`
        inINFO=`echo $LINE | cut -f 8`
        TUM=`echo $LINE | cut -f 10`
        NORM=`echo $LINE | cut -f 11`

        if [[ $inINFO =~ .*SOMATIC*. ]]
        then
           NDP=`echo $inINFO | cut -d\; -f 2 | cut -d= -f 2`
           nsc=`echo $inINFO | cut -d\; -f 7 | cut -d= -f 2`
           INFO=`echo $inINFO | cut -d\; -f 8`
           TDP=`echo $inINFO | cut -d\; -f 10 | cut -d= -f 2`
           tsc=`echo $inINFO | cut -d\; -f 15 | cut -d= -f 2`
        else
           NDP=`echo $inINFO | cut -d\; -f 2 | cut -d= -f 2` 
           nsc=`echo $inINFO | cut -d\; -f 7 | cut -d= -f 2`
           INFO="."
           TDP=`echo $inINFO | cut -d\; -f 9 | cut -d= -f 2`
           tsc=`echo $inINFO | cut -d\; -f 14 | cut -d= -f 2`
        fi

        ncf=`echo $nsc | cut -d, -f 1`
        ncr=`echo $nsc | cut -d, -f 2`
        nrf=`echo $nsc | cut -d, -f 3`
        nrr=`echo $nsc | cut -d, -f 4`
        tcf=`echo $tsc | cut -d, -f 1`
        tcr=`echo $tsc | cut -d, -f 2`
        trf=`echo $tsc | cut -d, -f 3`
        trr=`echo $tsc | cut -d, -f 4`

        nadr=`expr $nrf + $nrr`
        nadc=`expr $ncf + $ncr`
        tadr=`expr $trf + $trr`
        tadc=`expr $tcf + $tcr`
        NAD=`echo $nadr,$nadc`
        TAD=`echo $tadr,$tadc`

        NADF=`expr $ncf + $nrf`
        NADR=`expr $ncr + $nrr`
        TADF=`expr $tcf + $trf`
        TADR=`expr $tcr + $trr`

        NGT=`echo $NORM | cut -d: -f 1`
        TGT=`echo $TUM | cut -d: -f 1`

        FORMAT="GT:DP:AD:ADF:ADR"
        NORMAL=`echo $NGT:$NDP:$NAD:$NADF:$NADR`
        TUMOR=`echo $TGT:$TDP:$TAD:$TADF:$TADR`

        if [[ $FILTER == "." ]]
        then
            FILTER="PASS"
        else
            FILTER=`echo "PASS;$FILTER"`
        fi


        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$CHROM" "$POS" "$ID" "$REF" "$ALT" "$QUAL" "$FILTER" "$INFO" "$FORMAT" "$TUMOR" "$NORMAL" >> $output_vcf
    fi
done
unset $IFS
