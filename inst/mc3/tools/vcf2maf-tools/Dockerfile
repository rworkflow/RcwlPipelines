#################################################################
# Dockerfile
#
# Software:         vcf2maf
# Software Version: 1.6.13
# Description:      Convert a VCF into a MAF, where each variant is annotated 
#                   to only one of all possible gene isoforms
# Website:          https://github.com/OpenGenomics/vcft2maf-tools
# Base Image:       opengenomics/variant-effect-predictor-tool
# Run Cmd:          docker run opengenomics/vcf2maf perl vcf2maf.pl --man
#################################################################
FROM opengenomics/variant-effect-predictor-tool

MAINTAINER Adam Struck <strucka@ohsu.edu>

ENV PATH $VEP_PATH/htslib:$PATH
ENV PERL5LIB $VEP_PATH:/opt/lib/perl5:$PERL5LIB

WORKDIR /tmp

# install samtools
RUN apt-get update && \
    apt-get install --yes \
    libncurses5-dev vcftools && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN curl -L -o tmp.tar.gz https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
    mkdir samtools && \
    tar -C samtools --strip-components 1 -jxf tmp.tar.gz && \
    cd samtools && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf *

RUN curl -L -o tmp2.tar.gz https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 && \
    mkdir bcftools && \
    tar -C bcftools --strip-components 1 -jxf tmp2.tar.gz && \
    cd bcftools && \
    make && \
    make install && \
    cd .. && \
    rm -rf *

# install liftOver
RUN curl -L http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver > /usr/local/bin/liftOver && \
    chmod a+x /usr/local/bin/liftOver

# install vcf2maf
WORKDIR /opt/

RUN curl -ksSL -o tmp.tar.gz https://github.com/mskcc/vcf2maf/archive/v1.6.13.tar.gz && \
    tar --strip-components 1 -zxf tmp.tar.gz && \
    rm tmp.tar.gz && \
    chmod +x *.pl

CMD ["perl", "vcf2maf.pl", "--man"]
