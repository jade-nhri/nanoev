FROM docker.io/nvidia/cuda:12.1.0-devel-ubuntu20.04
ENV http_proxy=http://proxy.nhri.org.tw:3128
ENV https_proxy=http://proxy.nhri.org.tw:3128
ENV TZ=Asia/Taipei
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update && apt-get install -y \
    wget dos2unix \
    python3 python3-pip \
    python3-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    zlib1g-dev \
    libcurl4-gnutls-dev libssl-dev \
    libcudnn8 libcudnn8-dev \
    cmake unzip git wget libz-dev vim autoconf curl
#Install pyspoa
WORKDIR /opt
RUN wget https://github.com/nanoporetech/pyspoa/releases/download/v0.0.3/pyspoa-0.0.3.tar.gz
RUN tar -xzvf pyspoa-0.0.3.tar.gz
WORKDIR /opt/pyspoa
RUN make
RUN python3 -m pip install --upgrade pip
RUN pip install tensorflow==2.5.3
RUN pip3 install biopython pynput argparse pandas lxml six ete3 medaka==1.5.0
#Install blast2.13.0
WORKDIR /opt
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz
RUN tar zxvf ncbi-blast-2.13.0+-x64-linux.tar.gz
#Install minimap2
WORKDIR /opt
RUN  git clone https://github.com/lh3/minimap2
WORKDIR /opt/minimap2
RUN make
#Install htslib
WORKDIR /opt
RUN git clone https://github.com/samtools/htslib.git
WORKDIR /opt/htslib
RUN git submodule update --init --recursive
RUN autoreconf -i
RUN ./configure
RUN make
RUN make install
#Install bcftools
WORKDIR /opt
RUN git clone https://github.com/samtools/bcftools.git
WORKDIR /opt/bcftools
RUN make
#samtools 1.13
ADD https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 /opt
RUN apt-get update && apt-get install -y \
    libncurses-dev \
    apt-file \
    liblzma-dev \
    libz-dev \
    libbz2-dev \
    vim parallel
WORKDIR /opt
RUN tar -xjf /opt/samtools-1.13.tar.bz2
WORKDIR /opt/samtools-1.13
RUN make && make install
#Install seqkit
WORKDIR /opt
RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.2.0/seqkit_linux_amd64.tar.gz
RUN tar -xvf seqkit_linux_amd64.tar.gz
RUN cp seqkit /usr/local/bin
#Install nanoev
WORKDIR /opt
Run git clone https://github.com/jade-nhri/nanoev.git
RUN gunzip nanoev/EV.fasta.concor_20210303.fasta.gz
RUN chmod +x nanoev/*.py
WORKDIR /

ENV PATH $PATH:/opt:/opt/nanoev:/opt/ncbi-blast-2.13.0+/bin/:/opt/minimap2:/opt/htslib:/opt/bcftools:/opt/samtools-1.13

# Install miniconda
WORKDIR /opt
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
     /bin/bash ~/miniconda.sh -b -p /opt/conda
ENV PATH=$PATH:/opt/conda/bin
#Install homopolish
WORKDIR /opt
RUN git clone  https://github.com/ythuang0522/homopolish.git
WORKDIR /opt/homopolish
RUN chmod +x homopolish.py
RUN conda env create -f environment.yml
WORKDIR /
