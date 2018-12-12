FROM java

RUN   apt-get update && apt-get install -y wget python python-pip python-dev zlib1g-dev libbz2-dev liblzma-dev unzip libncurses-dev
WORKDIR /opt
RUN   wget -O samtools-0.1.18.tar.gz https://github.com/samtools/samtools/archive/0.1.18.tar.gz && tar -xvzf samtools-0.1.18.tar.gz && rm -f samtools-0.1.18.tar.gz
RUN   cd /opt/samtools-0.1.18 && make && cp /opt/samtools-0.1.18/samtools /usr/bin/

RUN   pip install pysam
RUN   pip install PyVCF


#Download and unpack radia
RUN   wget https://github.com/ucscCancer/radia/archive/v1.1.5.tar.gz && tar xvzf v1.1.5.tar.gz 

# blat
RUN   wget -O /usr/bin/blat http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v287/blat/blat && chmod +x /usr/bin/blat

# Download and unpack snpEff
RUN   wget http://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip && unzip snpEff_latest_core.zip

# cleanup
RUN   rm -rf v1.1.5.tar.gz snpEff_latest_core.zip snpEff/galaxy snpEff/examples 

COPY radia.py /opt/
COPY radia_filter.py /opt/
RUN chmod +x /opt/radia*
