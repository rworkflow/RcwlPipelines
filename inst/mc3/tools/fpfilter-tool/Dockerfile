FROM ubuntu:12.04

RUN apt-get update
RUN apt-get install -y wget build-essential git-core cmake zlib1g-dev libncurses-dev

WORKDIR /opt
RUN wget https://github.com/genome/bam-readcount/archive/v0.7.4.tar.gz && tar xvzf v0.7.4.tar.gz && rm -f v0.7.4.tar.gz
RUN wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2 && tar -xvjf samtools-1.2.tar.bz2 && rm -f samtools-1.2.tar.bz2

RUN cd /opt/bam-readcount-0.7.4 && mkdir build && cd build && cmake ../ && make deps && make -j && make install
RUN cd /opt/samtools-1.2 && make -j && make install
ADD fpfilter.pl /opt/fpfilter.pl