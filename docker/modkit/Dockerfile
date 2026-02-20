FROM ubuntu:20.04
MAINTAINER Mikhail Kolmogorov, mikolmogorov@gmail.com
RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
	apt-get install -y cmake git make gcc g++ autoconf bzip2 wget tabix libz-dev libncurses5-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev && \
    apt-get install -y python3-pip python3-virtualenv virtualenv python3-dev && \
	apt-get clean
RUN apt-get install -y bedtools

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

WORKDIR /opt/modkit
RUN wget https://github.com/nanoporetech/modkit/releases/download/v0.4.1/modkit_v0.4.1_u16_x86_64.tar.gz && \
    tar -xvf modkit_v0.4.1_u16_x86_64.tar.gz && \
    cp dist_modkit_v0.4.1_cec0a0b/modkit /usr/bin/modkit
