task maf2vcf {
    File inputMAF
    File refFasta
    File refFastaFai
    String? tumDepthCol
    String? tumRadCol
    String? tumVadCol
    String? nrmDepthCol
    String? nrmRadCol
    String? nrmVadCol
    Boolean perTumNrmVcfs = false
    String outputFilePrefix

    command {
        perl /opt/maf2vcf.pl --input-maf ${inputMAF} \                           
                              --ref-fasta ${refFasta} \
                              ${"--tum-depth-col " + tumDepthCol} \
                              ${"--tum-rad-col " + tumRadCol} \
                              ${"--tum-vad-col " + tumVadCol} \
                              ${"--nrm-depth-col " + nrmDepthCol} \
                              ${"--nrm-rad-col " + nrmRadCol} \
                              ${"--nrm-vad-col " + nrmVadCol} \
                              ${true="--per-tn-vcfs" false="" perTumNrmVcfs} \
                              --output-vcf ${outputFilePrefix}.vcf
    }

    output {
        File vcf = "${outputFilePrefix}.vcf"
    }

    runtime {
        docker: "opengenomics/vcf2maf"
    }
}

workflow run {
    call maf2vcf
}
