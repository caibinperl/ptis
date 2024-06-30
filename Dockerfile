FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive TZ=Asia/Shanghai

RUN apt-get update -qq \
    && apt-get -y --no-install-recommends install \
    build-essential cmake autoconf wget

RUN apt-get -y --no-install-recommends install \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libbz2-dev \
    libncurses5-dev \
    liblapack-dev \
    libblas-dev \
    libgfortran-10-dev \
    gfortran \
    zlib1g-dev \
    libcurl4-gnutls-dev \
    software-properties-common \
    dirmngr \
    tzdata

# Python
RUN apt-get -y --no-install-recommends install \
    python3-dev \
    python3-pip

RUN pip3 config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple
RUN python3 -m pip install --upgrade pip
RUN pip3 install pandas 
RUN pip3 install pysam 
RUN pip3 install snakemake
# resolve "AttributeError: module 'pulp' has no attribute 'list_solvers'"
RUN pip3 install pulp==2.7.0

# R
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | \
    tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && add-apt-repository -y "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
    && apt-get -y --no-install-recommends install r-base

RUN echo 'options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))' > ~/.Rprofile
RUN echo 'options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")' >> ~/.Rprofile

# Tidyverse
RUN Rscript -e 'install.packages("tidyverse")'

# fastp
RUN wget -O fastp http://opengene.org/fastp/fastp.0.23.4 \
    && chmod +x fastp \
    && mv fastp /usr/local/bin/

# bwa, samtools
RUN apt -y --no-install-recommends install bwa samtools

WORKDIR /ptis
COPY script ./script
COPY Snakefile ./
