# CLIPs4U
CLIPs4U is the snakemeke workflow to analyze PAR-CLIP data. 

## Overview

[Introduction](#introduction)

[Prerequisites](#prerequisites)

[Dependencies](#dependencies)

[Preparing configuration file](#preparing_configuration_file)

[Workflow](#workflow)

[Usage](#usage)

[Output](#output)

[Contributors](#contributors)

[License](#license)

## Introduction
CLIPs4U runs whole analysis starting from gzipped fastq files end ending with annotation tables using single command and yaml configuration file. It is user friendly and a highly customizable tool.


## Prerequisites
The only thing you must know before running CLIPs4U is the location of your fastq.gz files and your adapters sequences. If you want you can also provide your own genome version and annotation, star index etc (see below). 

## Dependencies
Required libraries are stored in the `env.yml` file. All required dependencies can be installed using *conda* by executing the following command:
```
conda env create --name clips4u --file env.yml
```
Where clips4u is your name of the environment.

**NOTE:** creating environment will take some time, a few hundred packages must be downloaded and installed.

The environment then needs to be activated in order to run CLIPs4U:
```
conda activate clips4u
```

If you prefere working with Docker, CLIPS4U docker image can be found [here](https://hub.docker.com/repository/docker/msajek/clips4u/general) 

## Preparing configuration file
DETAILS ABOUT CONFIGURATION FILE CAN BE FOUND IN config/README.md file


## Workflow

The example directed acyclic graph for the workflow is show below.

![DAG](img/dag.png)

Steps:
* downloading genome fasta and gtf files, prepare genome 2bit file, `main_annotation.rds` file and index for selected aligner (performs only once for selected genome, can be omitted if user specify paths to files in `config.yaml`)
* raw reads quality control using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* adapter trimming (one or two rounds - second round OPTIONAL) using [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html)
* reads collapsing ant OPTIONAL UMI removal using [seqtk](https://github.com/lh3/seqtk) and [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
* OPTIONAL removing repetitive elements using [STAR](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf)
* genome alignment using [bowtie](https://bowtie-bio.sourceforge.net/index.shtml) or [STAR](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf)
* running [PARalyzer](https://ohlerlab.mdc-berlin.de/files/duke/PARalyzer/README.txt) to detect enriched clusters 
* calculation of mismatches statistics and conversion specificity for the PARalyzer output
* annotation of clusters
* creating bigWig files for clusters using [deepTools](https://deeptools.readthedocs.io/en/develop/) for visualization in genome viewers
* motif enrichment analysis using [meme](https://meme-suite.org/meme/tools/meme), [dreme](https://meme-suite.org/meme/doc/dreme-tutorial.html) or [streme](https://meme-suite.org/meme/doc/streme-tutorial.html)
* generating final report and plots


## Usage
**First,** install git-lfs, to ensure that test data will be correctly downloaded.
Dependent on your system git-lfs can be installed using the following commands:
### Ubuntu/Debian
```
sudo apt-get update
sudo apt-get install git-lfs
```
### Fedora
```
sudo dnf update
sudo dnf install git-lfs
```
### CentOS7/RHEL7
```
sudo yum install epel-release
sudo yum install git-lfs
```
### CentOS8/RHEL8
```
sudo dnf install epel-release
sudo dnf install git-lfs
```
### openSUSE
```
sudo zypper refresh
sudo zypper install git-lfs
```
### Arch linux
```
sudo pacman -S git-lfs
```

Then initialize git-lfs using:
```
git lfs install
```

**Second,** clone the repository using: 
```
gh repo clone mukherjeelab/CLIPs4U
```
or:
```
git clone https://github.com/mukherjeelab/CLIPs4U.git
```
or just download zipped package and unpack it.

Ensure that `ZFP36.fq.gz` file in your `test_data` directory has 242MB. After navigating to your CLIPs4U directory and typing `ls -lh test_data` you should have output similar to the one presented below:
```
total 242M
-rw-r----- 1 marcin marcin 242M lip 11 20:26 ZFP36.fq.gz
-rw-r----- 1 marcin marcin 3,7K lip 12 12:46 ZFP36.yaml
```
If for some reasons your output looks like one below:
```
total 4,5K
-rw-r----- 1 sajekmar mukherjee  134 07-12 14:09 ZFP36.fq.gz
-rw-r----- 1 sajekmar mukherjee 3,7K 07-12 14:09 ZFP36.yaml
```
and ZFP36.fq.gz has only 134 B download it manually from [repository](https://github.com/mukherjeelab/CLIPs4U/blob/main/test_data/ZFP36.fq.gz) and move to `/path/to/CLIPs4U/test_data`.

**Third,** create directory for your project, e.g.:
```
mkdir my_parclip_dir
```

**Fourth,** prepare your config YAML file. 
It can be located anywhere, but if you put it in your directory it will be automatically detected.
Parameters not specified in your config file will be set to default values using default_config file.

Please note, that parameters that will be shared between all analyses might be put in the (clips4u)/config/default_config.yaml.
File default_config.yaml contains predefined default parameters, and you are free to change them. 

**Fifth,** create and activate conda environment as described above.

**Running analysis**
specify clips4u snakefile using flag "--snakefile",
specify your directory using flag "--directory" if it is not current working directory (if you are not in this directory),
specify your configfile using flag "--configfile" if configfile is located outside your directory or your directory contains multiple YAML files,
specify maximum number of threads using flag "--threads" 

example dry run command:
```
snakemake -n 1 --snakefile clips4u/workflow/Snakefile --directory my_parclip_dir
```
example test data run command:
```
snakemake --snakefile /path/to/CLIPs4U/workflow/Snakefile --configfile /path/to/CLIPs4U/test_data/ZFP36.yaml --threads 10
```
For further options and possibilities please check [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/index.html).

Test data contains fastq.gz file from ZFP36 PAR-CLIP described [here](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-1-r12), which can be also downloaded from [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53185).

## Output
Every run generates multiple output files. The most important one are:
* tabular TSV files with clusters annotation `workdir/annot`
* plots from various steps of analysis in pdf format `workdir/plots`, including annotation, mismatches statistics, motif enrichment, metagene plots
* final report in html format `workdir/final_report.html`, containing summary tables and plots from various steps of analysis in user friendly format, plots from final report in pdf format can be also found in `workdir/plots`
* bigWig files for visualization of the clusters in genomic viewers `workdir/genome_viewer_files`

**NOTE:** There are four bigWig (bw) files for every sample - filtered and unfiltered for positive and negative strand. Unfiltered files contain all clusters, filtered ones clusters with conversion specificity > 0.6 based. This value was set based on our previous [study](https://academic.oup.com/nar/article/47/2/570/5230955).


**NOTE:** For reproducibility purposes final_config.json will be created in your working directory. This file will contain all default and computed parameters alongside user defined parameters.

**NOTE:** Example output files for test data can be found in `test_out`. Please compare your output from test run with these files. Consider you will have more output files. In test_out we are storing only the most important ones.

## Contributors
* Marcin Sajek
* Neelanjan Mukherjee
* Samantha Lisy
* Yelena Prevalova
* Manuel Ascano Jr
* Tomasz Woźniak


## License
MIT license
