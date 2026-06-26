FROM ubuntu:22.04
MAINTAINER Mikhail Kolmogorov, mikolmogorov@gmail.com
RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
	apt-get install -y cmake git make gcc g++ autoconf bzip2 wget tabix libz-dev libncurses5-dev libbz2-dev liblzma-dev && \
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

RUN git clone https://github.com/KolmogorovLab/Severus && \
    cd Severus && \
    git checkout 163146a && \
    pip3 install -e .

RUN pip3 install "networkx>=2.6" pydot matplotlib numpy biopython pysam plotly
