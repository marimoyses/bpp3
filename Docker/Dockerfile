# Dockerfile for SV-Adjudicator
# Copyright (c) 2018 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT License

FROM continuumio/miniconda3
MAINTAINER "Harold Wang <haroldzwang@gmail.com>"

# Create python environment
RUN conda update -n base conda
RUN conda install -y -c bioconda python=3.6 numpy scipy cython natsort boto3 pandas scikit-learn

# Install make and zlib
RUN apt-get -y update && apt-get -y install --fix-missing build-essential zlib1g-dev unzip libpng-dev

ARG INSTALL_COMPONENTS

VOLUME ["/root/.config"]


# Install Ryan's WGD repo
# Install svtk
RUN git clone -b master https://github.com/talkowski-lab/bpp3.git /opt/bpp3 && \
    cd /opt/bpp3 && \
    git checkout c4aa51af0077462fdb834ac6cf349cc503d57335 && \
    pip install -e . && \
    cd Software && \
    tar -xzvf primer3-2.3.7.tar.gz && \
    cd primer3-2.3.7/src && \
    make all && \
    cd - && \
    unzip blatSrc35.zip && \
    cd blatSrc && \
    MACHTYPE=x86_64 && \
    export MACHTYPE && \
    mkdir -p ~/bin/x86_64 && \
    make
ENV PATH="/opt/bpp3/Software/primer3-2.3.7/src:/root/bin/x86_64:${PATH}"
RUN apt-get install -y libbz2-dev liblzma-dev libcurl4-openssl-dev
RUN conda install -y -c conda-forge -c bioconda bzip2


# Launch bash
CMD ["/bin/bash"]

