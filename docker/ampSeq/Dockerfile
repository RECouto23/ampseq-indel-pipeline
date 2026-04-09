FROM continuumio/miniconda3:24.11.1-0

LABEL maintainer="RECouto524@gmail.com"
LABEL description="AmpSeq pipeline: FastQC, Trimmomatic, FLASH, BWA, Samtools, Python analysis tools"
LABEL version="1.0.0"

# Set channels and strict channel priority
RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict

# Install all packages except flash (no arm64 bioconda binary available)
RUN conda install -y \
    python=3.13.12 \
    bwa=0.7.19 \
    samtools=1.23.1 \
    fastqc=0.12.1 \
    trimmomatic=0.40 \
    pysam=0.23.3 \
    pandas=3.0.2 \
    numpy=2.4.3 \
    matplotlib=3.10.8 \
    seaborn=0.13.2 \
    openpyxl=3.1.5 \
    scipy=1.17.1 \
    biopython=1.87 \
    pyfaidx=0.9.0.4 \
    pyvcf3=1.0.3 \
    statsmodels=0.14.6 \
    make \
    zlib \
    && conda clean -afy

# Build FLASH2 from source (no arm64 bioconda binary available for flash 1.2.11)
# FLASH2 is the maintained successor and is functionally equivalent for this pipeline
RUN apt-get update && apt-get install -y gcc curl zlib1g-dev && apt-get clean && \
    curl -L https://github.com/dstreett/FLASH2/archive/refs/tags/2.2.00.tar.gz | tar xz && \
    cd FLASH2-2.2.00 && \
    make && \
    mv flash2 /usr/local/bin/flash && \
    cd .. && rm -rf FLASH2-2.2.00

# Expose Trimmomatic adapter files at a stable, hardcoded path
RUN TRIM_SHARE=$(find /opt/conda -name "TruSeq3-PE.fa" 2>/dev/null | head -1 | xargs dirname) && \
    mkdir -p /usr/local/share/trimmomatic && \
    cp ${TRIM_SHARE}/*.fa /usr/local/share/trimmomatic/

ENV TRIMMOMATIC_ADAPTERS=/usr/local/share/trimmomatic

# Verify key tools are on PATH
RUN fastqc --version && \
    bwa 2>&1 | head -3 && \
    samtools --version | head -1 && \
    flash --version && \
    trimmomatic -version
