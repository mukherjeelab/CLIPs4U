FROM ubuntu:22.04

LABEL base.image="ubuntu:22.04"
LABEL version="1"
LABEL software="CLIPs4U"
LABEL software.version="0.1.0"
LABEL software.description="CLIPs4U"
LABEL software.website="https://github.com/mukherjeelab/CLIPs4U/tree/main"
LABEL software.documentation="https://github.com/mukherjeelab/CLIPs4U/tree/main"
LABEL software.tags="Transcriptomics, PAR-CLIP, RBPs"
LABEL maintainer="marcin.sajek@cuanschutz.edu"
LABEL maintainer.organisation="University of Colorado School of Medicine"
LABEL maintainer.location="RC1S 10403A, Research Complex 1 South, 12801 E 17th Ave, Aurora, CO 80045, USA"
LABEL maintainer.lab="Mukherjee Lab"

# Update apt and install necessary packages
RUN apt-get update \
    && apt-get install -y wget git git-lfs bzip2 build-essential \
    && wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh

# Install Mambaforge
ENV PATH="/root/mambaforge/bin:$PATH"
RUN mkdir /root/.conda && bash Mambaforge-Linux-x86_64.sh -b -p /root/mambaforge

# Create a working directory
RUN mkdir /workspace

# Add the environment file to the working directory
ADD env.yml /workspace/env.yml

# Set the working directory
WORKDIR /workspace

# Create conda environment with Mamba
RUN mamba init bash \
    && . ~/.bashrc \
    && mamba env create -f env.yml

# Activate environment
RUN echo "conda activate CLIPs4U_1" >> ~/.bashrc

# Clone the repository into the working directory
RUN git lfs install \
    && git clone https://github.com/mukherjeelab/CLIPs4U.git

# Set the working directory to the cloned repository
WORKDIR /workspace
