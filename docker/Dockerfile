FROM continuumio/miniconda3
WORKDIR /app

# biofilter version argument
ARG BIOFILTER_VERSION=2.4.3

RUN apt-get update \    
    # install packages needed to install GWAMA, biofilter, and NEAT-plots
    && apt-get install -y --no-install-recommends libz-dev g++ gcc git wget tar unzip make \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    # install GWAMA
    && wget https://www.geenivaramu.ee/tools/GWAMA_v2.2.2.zip \
    && unzip GWAMA_v2.2.2.zip \
    && make \
    # install python packages needed for pipeline
    && conda install -y -n base -c conda-forge wget libtiff conda-build scipy pandas seaborn matplotlib numpy apsw sqlite \
    && conda clean --all --yes \
    # install NEAT-plots
    && git clone https://github.com/path/to/directory/NEAT-Plots.git \
    && mv NEAT-Plots/manhattan-plot/ /app/ \
    && conda develop /app/manhattan-plot/ \
    # install biofilter
    && wget https://github.com/path/to/directory/biofilter-${BIOFILTER_VERSION}.tar.gz -O biofilter.tar.gz \
    && tar -zxvf biofilter.tar.gz --strip-components=1 -C /app \
    && /opt/conda/bin/python setup.py install \
    # make biofilter executable
    && chmod a+rx /app/biofilter.py \
    # remove biofilter tarball and NEAT-plots directory
    && rm -R biofilter.tar.gz NEAT-Plots

USER root
