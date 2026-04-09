FROM continuumio/miniconda3:24.11.1-0

LABEL maintainer="RECouto524@gmail.com"
LABEL description="Indel visualization: pysam, matplotlib, seaborn, pandas, openpyxl"
LABEL version="1.0.0"

# Set channels and strict channel priority
RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict

# Install all packages with exact versions from excelManip.yml
# Note: Python 3.14 is used here to match your excelManip environment exactly.
# If you encounter solver conflicts, fall back to python=3.13
RUN conda install -y \
    python=3.13 \
    pysam=0.23.3 \
    pandas=3.0.2 \
    numpy=2.4.3 \
    matplotlib=3.10.8 \
    seaborn=0.13.2 \
    openpyxl=3.1.5 \
    scipy=1.17.1 \
    statsmodels=0.14.6 \
    && conda clean -afy

# Install pip-only packages from excelManip.yml
# argparse is part of the standard library in Python 3 but listed explicitly;
# glob2 adds enhanced glob support used by IndelPlot.py
RUN pip install --no-cache-dir \
    glob2==0.7

# Verify imports
RUN python -c "import pysam, pandas, numpy, matplotlib, seaborn, openpyxl, glob2; print('All imports OK')"
