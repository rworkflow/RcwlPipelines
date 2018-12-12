FROM ubuntu:17.10

RUN apt-get update 
RUN apt-get install -y git python python-pip zlib1g-dev
RUN pip install PyYAML

WORKDIR /opt
ADD tcga-vcf-reheader.py /opt/tcga-vcf-reheader.py
ADD reheader_wrapper.sh /opt/reheader_wrapper.sh

