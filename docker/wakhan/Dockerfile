FROM ubuntu:20.04
MAINTAINER Mikhail Kolmogorov, mikolmogorov@gmail.com
RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
	apt-get install -y cmake git make gcc g++ autoconf bzip2 wget tabix libz-dev libncurses5-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev && \
    apt-get install -y python3-pip python3-virtualenv virtualenv python3-dev && \
	apt-get clean

### samtools
# 1.15
WORKDIR /opt/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2 && \
    tar xvf samtools-1.15.tar.bz2 && \
	rm -r /opt/samtools/samtools-1.15.tar.bz2 && \
	cd samtools-1.15/ && \
	autoheader && \
	autoconf -Wno-header && \
	./configure && \
	make && \
	cp samtools /usr/bin/samtools

#bcftools
ARG bcftoolsVer="1.19"
WORKDIR /opt/bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/${bcftoolsVer}/bcftools-${bcftoolsVer}.tar.bz2 && \
    tar -vxjf bcftools-${bcftoolsVer}.tar.bz2 && \
    rm bcftools-${bcftoolsVer}.tar.bz2 && \
    cd bcftools-${bcftoolsVer} && \
    make && \
    make install

RUN pip3 install pysam pyfaidx numpy pandas plotly ruptures scikit-learn==1.2.0 scipy==1.9.2 vcf_parser kaleido

WORKDIR /opt/wakhan
RUN git clone https://github.com/KolmogorovLab/Wakhan && \
    cd Wakhan && git checkout 02060ae && \
    echo '#!/bin/bash\npython3 /opt/wakhan/Wakhan/wakhan.py "$@"' > /usr/bin/wakhan && \
    chmod +x /usr/bin/wakhan
