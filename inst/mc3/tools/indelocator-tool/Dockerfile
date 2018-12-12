FROM    java:7

RUN     apt-get update
RUN     apt-get install -y wget build-essential git-core zlib1g-dev libncurses-dev python python-pip
RUN     pip install PyVCF

WORKDIR /opt
ENV     JAVA_JAR_PATH   /opt
RUN     wget https://software.broadinstitute.org/cancer/cga/sites/default/files/data/tools/indelocator/IndelGenotyper.36.3336-GenomeAnalysisTK.jar && mv IndelGenotyper.36.3336-GenomeAnalysisTK.jar IndelGenotyper.jar
ADD reformat_vcf.sh /opt/reformat_vcf.sh
