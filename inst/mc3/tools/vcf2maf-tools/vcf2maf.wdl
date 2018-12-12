task vcf2maf {
    File inputVCF
    File? vepAnnotatedInputVCF
    String tmpDir = "."
    File vepOfflineCacheDir
    File refFasta
    File refFastaFai
    File? remapChain
    String ncbiBuild = "GRCh37"
    String species = "homo_sapiens"
    Array[String]? retainInfoCols
    String? tumorId
    String? normalId
    String? vcfTumorId
    String? vcfNormalId
    File? customEnst
    String? mafCenter
    Float? minHomVaf
    
    String outputFilePrefix

    command {
        if [ -n "${vepAnnotatedInputVCF}" ]; then
            ln -s ${vepAnnotatedInputVCF} ${tmpDir}/$(basename ${inputVCF} .vcf).vep.vcf
        fi

        # retain extra INFO cols
        if [ -n "${sep="," retainInfoCols}" ]; then
           INFOCOLS="--retain-info  ${sep ="," retainInfoCols}"
        else
           INFOCOLS=""
        fi

        perl /opt/vcf2maf.pl --input-vcf ${inputVCF} \
                              --output-maf ${outputFilePrefix}.maf \
                              --vep-path /root/vep
                              --vep-data ${vepOfflineCacheDir} \
                              --ref-fasta ${refFasta} \
                              --species ${species} \
                              --ncbi-build ${ncbiBuild} \
                              --tmp-dir ${tmpDir} \
                              ${"--remap-chain " + remapChain} \
                              ${"--maf-center " + mafCenter} \
                              ${"--tumor-id " + tumorId} \
                              ${"--normal-id " + normalId} \
                              ${"--vcf-tumor-id " + vcfTumorId} \
                              ${"--vcf-normal-id " + vcfNormalId} \
                              ${"--custom-enst " + customEnst} \
                              ${"--min-hom-vaf " + minHomVaf} \
                              $INFOCOLS

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
    call vcf2maf
}
