FROM continuumio/miniconda3:24.11.1-0

LABEL maintainer="RECouto524@gmail.com"
LABEL description="Indel counting: pysam, biopython, scikit-bio, pandas, scipy, ipython"
LABEL version="1.0.0"

# Set channels and strict channel priority
RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict

# Install all packages with exact versions from tsaiIndels.yml
RUN conda install -y \
    python=3.11.15 \
    pysam=0.23.3 \
    pandas=2.3.3 \
    numpy=1.26.4 \
    biopython=1.87 \
    scipy=1.10.0 \
    pyfaidx=0.9.0.4 \
    pyvcf3=1.0.4 \
    scikit-bio=0.5.9 \
    biom-format=2.1.17 \
    matplotlib-base=3.10.8 \
    statsmodels=0.14.6 \
    natsort=8.4.0 \
    h5py=3.16.0 \
    ipython=9.10.1 \
    && conda clean -afy

# Verify import of core dependencies
RUN python -c "import pysam, pandas, numpy, Bio, scipy, skbio; print('All imports OK')"
