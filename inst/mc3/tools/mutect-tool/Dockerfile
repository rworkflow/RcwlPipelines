
FROM java:7

RUN apt-get update
RUN apt-get install -y zip wget
#
# Install samtools and the python vcf libraries for wrapper code
#
RUN apt-get install -y samtools python-pip
RUN pip install PyVCF

# We'll be working in /opt from now on
WORKDIR /opt

#
# Download and unpack Mutect
#
RUN wget https://github.com/broadinstitute/mutect/releases/download/1.1.5/muTect-1.1.5-bin.zip
RUN wget https://github.com/broadinstitute/picard/releases/download/1.122/picard-tools-1.122.zip

RUN unzip muTect-1.1.5-bin.zip
RUN unzip picard-tools-1.122.zip
COPY mutect.py /opt/

# Link the picard tools to /opt/picard
RUN ln -s picard-tools-1.122 picard
