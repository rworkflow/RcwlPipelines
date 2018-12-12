FROM ubuntu


RUN     apt-get update
RUN     apt-get install -y git wget make gcc zlib1g-dev ncurses-dev g++ python python-pip
RUN     pip install PyVCF subprocess32

WORKDIR /opt

RUN		mkdir /opt/bin
RUN     wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2
RUN     tar xvjf samtools-1.2.tar.bz2
RUN     cd /opt/samtools-1.2 && make && make install

RUN     git clone https://github.com/genome/pindel.git
RUN     cd pindel && git checkout v0.2.5b8
RUN     cd pindel && ./INSTALL /opt/samtools-1.2/htslib-1.2.1
RUN     cp /opt/pindel/pindel* /usr/local/bin/
ADD     pindel.py /opt/
