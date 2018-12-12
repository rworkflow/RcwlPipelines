task vcf2vcf {
    File inputVCF
    String? tumorId
    String? normalId
    String? vcfTumorId
    String? vcfNormalId
    String outputFilePrefix

    command {
        perl /opt/vcf2vcf.pl --input-vcf ${inputVCF} \
                              --output-vcf ${outputFilePrefix}.vcf \
                              ${"--tumor-id " + tumorId} \
                              ${"--normal-id " + normalId} \
                              ${"--vcf-tumor-id " + vcfTumorId} \
                              ${"--vcf-normal-id " + vcfNormalId}
    }

    output {
        File maf = "${outputFilePrefix}.vcf"
    }

    runtime {
        docker: "opengenomics/vcf2maf"
    }
}

workflow run {
    call vcf2vcf
}
