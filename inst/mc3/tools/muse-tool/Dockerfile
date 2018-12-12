FROM    ubuntu

RUN     apt-get update && apt-get install -y python wget samtools tabix python-pip python-dev zlib1g-dev
RUN     pip install pysam

RUN     mkdir /opt/bin
RUN     wget -O /opt/bin/MuSEv0.9.9.5 http://bioinformatics.mdanderson.org/Software/MuSE/MuSEv0.9.9.5
RUN     wget -O /opt/bin/MuSEv1.0rc http://bioinformatics.mdanderson.org/Software/MuSE/MuSEv1.0rc_submission_c039ffa
RUN     chmod +x /opt/bin/MuSE*
ADD     ./muse.py /opt/bin/
ADD     ./vcf_reformat.py /opt/bin/vcf_reformat.py
ADD     hgsc_vcf /opt/bin/hgsc_vcf
ADD     filter_muse.py /opt/bin/filter_muse.py
ENV     PATH /usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/bin
