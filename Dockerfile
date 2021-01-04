FROM continuumio/miniconda3
# continuumio/miniconda3 is FROM debian:latest

# Make RUN commands use `bash --login` (always source ~/.bashrc on each RUN)
SHELL ["/bin/bash", "--login", "-c"]

# install apt depedencies and update conda
RUN apt-get update && apt-get install git -y \
    && apt-get install -y apt-transport-https ca-certificates wget unzip bzip2 libfontconfig1 \
    && update-ca-certificates \
    && apt-get -qq -y remove curl \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && apt-get install -y dos2unix \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log 

# link global conda to expected sciserver location
#RUN ln -s /opt/conda/ /home/idies/workspace/covid19/miniconda3

ENV PATH /opt/conda/bin:$PATH

RUN conda install -y python=3 \
    && conda update -y conda \
    && conda clean --all --yes


# configure directory structure exactly as it is on SciServer for ease of transition

RUN mkdir -p /root/idies/workspace/covid19/code \
    && chmod g+s /root/idies/workspace/covid19/code

# install openjdk
WORKDIR /root/idies/workspace/covid19/code
RUN wget https://download.java.net/java/GA/jdk14/076bab302c7b4508975440c56f6cc26a/36/GPL/openjdk-14_linux-x64_bin.tar.gz \
    && tar -xzf openjdk-14_linux-x64_bin.tar.gz \
    && rm openjdk-14_linux-x64_bin.tar.gz

# install samtools 1.10
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 \
    && tar -xjf samtools-1.10.tar.bz2 \
    && rm samtools-1.10.tar.bz2

# clone VCF IGV repo and install local IGV
RUN git clone https://github.com/mkirsche/vcfigv \
    && wget https://data.broadinstitute.org/igv/projects/downloads/2.8/IGV_2.8.10.zip -P vcfigv \
    && unzip -d vcfigv vcfigv/IGV_2.8.10.zip \
    && rm vcfigv/IGV_2.8.10.zip

# install ONT guppy tools
RUN wget --no-check-certificate https://americas.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_3.6.1_linux64.tar.gz \
    && tar -xzf ont-guppy-cpu_3.6.1_linux64.tar.gz \
    && rm ont-guppy-cpu_3.6.1_linux64.tar.gz


# Clone Artic and Pangolin repos and create environments to install packages
RUN git clone --recurse-submodules https://github.com/artic-network/artic-ncov2019 \
    && git clone https://github.com/cov-lineages/pangolin.git \
    && conda env create -f artic-ncov2019/environment.yml \
    && sed -i 's/  - python=3.6/  - python=3.7/' pangolin/environment.yml \
    && conda env create -f pangolin/environment.yml \
    && conda activate pangolin \
    && cd pangolin \
    && python setup.py install


# Clone the Variant and CoverageNormalization repos
WORKDIR /root/idies/workspace/covid19/code/ncov/pipeline_scripts
RUN git clone https://github.com/mkirsche/CoverageNormalization.git \
    && git clone https://github.com/mkirsche/VariantValidator.git


# clone cov-lineages
WORKDIR /root/idies/workspace/covid19/ncov_reference
RUN git clone https://github.com/cov-lineages/lineages.git

#RUN find /root/idies/workspace/covid19/ -type f -print0 | xargs -0 dos2unix

# compile java libraries
ENV PATH="/root/idies/workspace/covid19/code/ncov/pipeline_scripts:${PATH}"
ENV PATH="/root/idies/workspace/covid19/code/jdk-14/bin:${PATH}"
ENV PATH="/root/idies/workspace/covid19/code/ont-guppy-cpu/bin:${PATH}"


# install TeX libraries
RUN wget -qO- "https://yihui.org/tinytex/install-bin-unix.sh" | sh \
    && export PATH=/root/.TinyTeX/bin/x86_64-linux:$PATH \
    && tlmgr path add \
    && tlmgr install mnsymbol \
    && tlmgr install multirow \
    && tlmgr install wrapfig \
    && tlmgr install colortbl \
    && tlmgr install pdflscape \
    && tlmgr install tabu \
    && tlmgr install threeparttable \
    && tlmgr install threeparttablex \
    && tlmgr install environ \
    && tlmgr install ulem \
    && tlmgr install makecell 
##################################################################

#Copy just the environment.yml file for quick debugging purposes. Comment this out in production
COPY ./environment.yml /root/idies/workspace/covid19/code/ncov/

#Finally, copy over the local files into the working directory and copy rest of necessary environment over to workspace
#COPY ./ /root/idies/workspace/covid19/code/ncov/


RUN conda env create -f /root/idies/workspace/covid19/code/ncov/environment.yml 

# Re-copy yml file and the rest for quick debugging. Comment this out in production
COPY ./ /root/idies/workspace/covid19/code/ncov/ \
    && cp -r /root/idies/workspace/covid19/code/ncov/covid19 /root/idies/workspace/



#################################################################
# copy 13-gene genome.json over into RAMPART directory (which has a 9-gene file by default)
RUN cp /root/idies/workspace/covid19/ncov_reference/genome.json /root/idies/workspace/covid19/code/artic-ncov2019/rampart/genome.json

RUN ls /root/idies/workspace/covid19/code/ncov/pipeline_scripts
# set up final environment and default working directory
RUN cat /root/idies/workspace/covid19/bashrc >> /root/.bashrc
RUN cat /root/.bashrc
RUN chmod -R 755 /root/idies/workspace/covid19/code/ncov/pipeline_scripts/
WORKDIR /root/idies/workspace/covid19
