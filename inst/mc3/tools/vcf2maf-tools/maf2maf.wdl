task maf2maf {
    File inputMAF
    File? vepAnnotatedInputVCF
    File? customEnst
    String? tumDepthCol
    String? tumRadCol
    String? tumVadCol
    String? nrmDepthCol
    String? nrmRadCol
    String? nrmVadCol
    Array[String]? retainCols
    File vepOfflineCacheDir
    File refFasta
    File refFastaFai
    String ncbiBuild = "GRCh37"
    String species = "homo_sapiens"
    String outputFilePrefix
    String tmpDir = "."

    command {
        if [ -n '${vepAnnotatedInputVCF}' ]; then
            ln -s ${vepAnnotatedInputVCF} ${tmpDir}/$(basename ${inputMAF} .maf).vep.vcf
        fi
    
        # Handling the optional arrays
        if [ -n '${sep="," retainCols}' ]; then
            RETAIN_COLS_FLAG=--apply-filters ${sep="," retainCols}
        else
            RETAIN_COLS_FLAG=''
        fi

        perl /opt/maf2maf.pl --input-maf ${inputMAF} \
                              --output-maf ${outputFilePrefix}.maf
                              --tmp-dir ${tmpDir} \
                              ${"--tum-depth-col " + tumDepthCol} \
                              ${"--tum-rad-col " + tumRadCol} \
                              ${"--tum-vad-col " + tumVadCol} \
                              ${"--nrm-depth-col " + nrmDepthCol} \
                              ${"--nrm-rad-col " + nrmRadCol} \
                              ${"--nrm-vad-col " + nrmVadCol} \
                              ${"--custom-enst " + customEnst} \
                              $RETAIN_COLS_FLAG \
                              --vep-data ${vepOfflineCacheDir} \
                              --ref-fasta ${refFasta} \                              
                              --species ${species} \
                              --ncbi-build ${ncbiBuild}

       rm ${tmpDir}/$(basename ${inputVCF} .vcf).vep.vcf
    }

    output {
        File maf = "${outputFilePrefix}.maf"
    }

    runtime {
        docker: "opengenomics/vcf2maf"
    }
}

workflow run {
    call maf2maf
}
