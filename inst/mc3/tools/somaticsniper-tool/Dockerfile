FROM ubuntu:16.04

RUN apt-get update && apt-get install -y \
    wget \
    build-essential \
    git-core \
    cmake \
    zlib1g-dev \
    libbz2-dev \
    libncurses-dev \
    python \
    python-dev \
    python-pip

RUN	pip install pysam

WORKDIR /opt
RUN wget https://github.com/genome/somatic-sniper/archive/v1.0.5.0.tar.gz && tar xvzf v1.0.5.0.tar.gz && rm v1.0.5.0.tar.gz

RUN cd /opt/somatic-sniper-1.0.5.0 && mkdir build && cd build && cmake ../ && make deps && make -j && make install
RUN cp /opt/somatic-sniper-1.0.5.0/build/vendor/samtools/samtools /usr/bin/

COPY SomaticSniper.py /opt/ 
RUN chmod +x /opt/SomaticSniper.py
