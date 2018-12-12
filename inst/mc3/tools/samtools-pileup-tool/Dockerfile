FROM    ubuntu

RUN     apt-get update
RUN     apt-get install -y python wget python-pip python-dev zlib1g-dev libncurses-dev
RUN     pip install bx-python
WORKDIR /opt
RUN     wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2 && tar -xvjf samtools-1.2.tar.bz2 && rm -f samtools-1.2.tar.bz2
RUN     cd /opt/samtools-1.2 && make -j && make install
