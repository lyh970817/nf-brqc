FROM rocker/r-ver:4.3.0

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    wget \
    curl \
    git \
    procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Add conda to path
ENV PATH="/opt/conda/bin:${PATH}"

# Copy conda environment file
COPY conda/environment.yml /

# Create conda environment
RUN conda env create -f /environment.yml && conda clean -a

# Add conda environment to path
ENV PATH="/opt/conda/envs/bioresource-qc/bin:${PATH}"

# Install R packages
RUN R -e "install.packages(c('data.table', 'caret', 'pROC', 'verification', 'ggplot2', 'cowplot', 'MLmetrics', 'glmnet', 'optparse', 'tidyverse'), repos='https://cloud.r-project.org/')"

# Set working directory
WORKDIR /data

# Entry point
ENTRYPOINT ["nextflow"]
CMD ["run", "main.nf", "--help"]
