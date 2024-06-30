# PTIS: A Pipeline for Fast Identifying T-DNA Integration Sites in Plant Genomes

Whole genome sequencing (WGS) has become a powerful and popular approach for detecting transfer-DNA (T-DNA) integration sites in the *Agrobacterium tumefaciens*-mediated T-DNA transformation system in plant genetics and breeding. 
We present the PTIS pipeline, which detects *Agrobacterium*-mediated T-DNA integration sites in transgenic plants. 
PTIS offers several improvements compared to similar tools: 
- PTIS covers the entire analysis process, starting from preprocessing raw reads to outputting integration sites.
- Specifically, PTIS is packaged as a Docker image, making it easy to install and can be executed on various operating systems with a single command.
- Additionally, the execution time of PTIS is faster than other pipelines.

PTIS is a user-friendly, fast, Docker-based, end-to-end pipeline for efficiently identifying *Agrobacterium*-mediated T-DNA integration sites in transgenic plants. 
The open-source code is available in the GitHub repository: https://github.com/caibinperl/ptis.git. 
The Docker image can be accessed publicly from Docker Hub at https://hub.docker.com/r/caibinperl/ptis.

## Background
In plant genetics and breeding, the *Agrobacterium tumefaciens*-mediated transfer-DNA (T-DNA) transformation system is a powerful tool for generating transgenic plants. 
T-DNA is a fragment of the *A. tumefaciens* tumor-inducing plasmid, which is flanked by two 24-bp borders (left border and right border). 
*A. tumefaciens* can deliver the tumor-inducing plasmid and transfer DNA to susceptible plant cells as single-stranded DNA. 
This DNA is targeted to the cell nucleus, where the T-strand is probably converted to double strands and integrated into the host genome through non-homologous end-joining (NHEJ).
Most likely, NHEJ leads to a small microhomology of up to 5 bp overlapping between the host and vector sequence or incorporation of short random DNA fragments (referred to as fillers).
In rare cases, T-DNA is integrated precisely without creating microhomology or fillers.
The affordable cost and high throughput of Next-Generation Sequencing (NGS) make it possible to detect T-DNA integration sites from ‘split reads’ of whole genome sequencing (WGS) from transgenic plants.
A split read is a scenario where one portion of a read aligns with the plant genome while the remaining part maps to the vector genome.

We have developed a pipeline called PTIS that automates quality control on raw reads, trimming adapters, indexing the reference genome, aligning sequences, and detecting integration sites.
PTIS is an end-to-end pipeline that covers all analysis tasks, enabling it to be initiated from raw reads. 
PTIS presents results within some TSV files, which display the integration site position, type, read amount, etc. 
PTIS allows users to inspect informative read alignments using the Integrative Genomics Viewer (IGV). 
To utilize Docker's methodologies, we offer PTIS as a Docker image, streamlining the installation process by pulling the image and executing it with a single command. 
This is achieved by attaching the data volumes and running the Docker container. 

## Implementation
PTIS is built as a Docker container image and shared through Docker Hub (https://hub.docker.com/r/caibinperl/ptis). 
The container can run on various operating systems, such as Linux, Mac, and Windows, all supported by the Docker engine. 
The pipeline encompasses a series of essential procedures: conducting quality control on the raw reads, removing adapters, establishing an index for the reference genome, aligning sequences, identifying integration sites, and visualizing the resulting output.

### Preprocessing
WGS reads that have low-quality scores could result in false integration sites. 
For this reason, PTIS begins by performing quality control, adapter trimming, and other data-filtering operations using fastp with default options. 
To reduce the size of read files, we mapped the reads to the vector sequence and extracted all informative reads aligned to the vector.
These filtered reads are prepared for further analysis.

### Align reads to vector and host genomes

The plant reference genome and the transformed binary vector genome FASTA files are first indexed with BWA.
The filtered reads file is then aligned to the indexed FASTA file using BWA MEM with the default options. 
The resulting Sequence Alignment Map (SAM) or Binary Alignment Map (BAM) file contains information about reads that are partially or fully mapped to the vector or host genome.

### Load alignment from SAM

In SAM/BAM format, the alignments are represented by Compact Idiosyncratic Gapped Alignment Report (CIGAR) strings. 
The CIGAR is composed of several operations, including M (sequence match or mismatch), I (insertion to the reference), D (deletion from the reference), and S (clipped sequences). 
CIGAR is not an ideal method for determining breakpoints on a read that spans both the host and vector genome. 
For this reason, we transform CIGAR into a coordinate-based alignment. 
The following definitions were used to describe each alignment. 
Here, q_name, q_size, q_start, and q_end refer to reads name, reads sequence size, alignment start and end position in query respectively; t_name, t_size, t_start, and t_end represent target sequence name, target sequence size, alignment start and end position on host or vector genomes, respectively; the column ‘strand’ defined as ‘+’ (forward) or ‘-’ (reverse) for query strand. 
Alignment also includes a global identity score, defined as $q\_size-ed/q\_size$, for quality control on alignments. It also has a local identity score defined as $M+I-ed/M+I$, where ed represents the edit distance from the reference obtained from the tag NM in SAM. 
The alignment information is saved as a tab-separated values (TSV) file for further analysis.

### Filter alignment

The TSV files containing alignment data are loaded into a tibble (data frame) using the `read_csv()` function from the readr library in R version 4.4. 
The dplyr R library filters reads based on local identity, with a minimum requirement of 90%, corresponding to `min_local_identity` in the configuration file.

### Determine split reads from alignments

Within a reliable split read, three informative regions are required: a portion aligned to the plant genome (referred to as the host portion), a portion aligned to the vector genome (referred to as the vector portion), and an overlapping or insertion between the host and vector portions, if microhomology or filler exists. 
Utilizing the position in the alignment, we can filter split reads and determine junction sites from these reads.
To reduce the chance of false positives resulting from split reads with homologous DNA fragments between the vector and the host, we filtered reads with the entire reads mapped to either the plant or the vector genome.

For split reads, the junction site between the plant genome and vector T-DNA is determined by either the distance between the alignment start position on the vector and the end position on the host or the distance between the alignment start position on the host and the end position on the vector. 

If the host portion in a split read is directly adjacent to the vector portion, integration occurs without using microhomology and without creating fillers, resulting in zero. 
A negative distance describes a situation where two alignment parts in a split read overlap. 
Microhomology was detected with a negative distance shorter than the `max_microhomology_len` cutoff predefined in the configuration file. 
A positive distance was assigned if there was a gap between the vector and host portions. 
With a threshold (`max_filler_len`), the gap can be considered a filler.

### Determine integration sites

An informative split read supports an integration site; the above step will produce multiple sites close to each other.
If these sites are located within a `cluster_dist` range (predefined in the configuration file), they will be clustered into one integrated site. 
This integration is supported by multiple split reads.

### View of integration sites
The outputs for each sample are stored in tabular-separated value files, “site.tsv” and “site_aln.tsv”, representing integration sites and alignments. 
Using IGV's features, users can easily view and explore alignment data from indexed BAM files.


## Installation

### Requirements
Running the pipeline requires a recent version of the docker client, found on the docker website.
Follow the instructions to [download and install Docker](https://docs.docker.com/get-docker/).


## Usage

The PTIS Docker container includes the complete environment, scripts, and software (along with their respective versions). 
This eliminates any concerns about compatibility issues when installing bioinformatics tools on a platform different from the one they were developed on. 
The execution of PTIS depends on the Docker engine, which needs to be installed on the host machine. 
The installation of Docker on different operating systems is easy and described on the Docker website (https://docs.docker.com/get-docker). 
With Docker, PTIS can be executed on various operating systems with a simple command: `docker run --mount type=bind,src=data-dir,target=/ptis/data caibinperl/ptis snakemake --cores 4`. 
The `data-dir` refers to the data directory that contains configuration files, metadata, and raw reads. 

**The `data-dir` required by running PTIS is stored on the host's file system**. 

### Data directory (data-dir)
A data directory should be presented as the following structure: raw reads, host and vector references, and configure files.

The following will introduce the configuration and sequence of files using the data directory in the code source repository (https://github.com/caibinperl/ptis.git).

```sh
data
├── config
│   ├── config.yml
│   └── samples.tsv
├── fastq
│   ├── SRR16230416_R1.fastq.gz
│   ├── SRR16230416_R2.fastq.gz
│   ├── SRR16230426_R1.fastq.gz
│   └── SRR16230426_R2.fastq.gz
└── resource
    ├── 46bg.fasta
    ├── 49btm.fasta
    └── osa.fasta
```

### fastq

* `SRR16230416_R1.fastq.gz` and `SRR16230416_R2.fastq.gz` are raw reads from NGS sequencing for sample `SRR16230416`.
`SRR16230426_R1.fastq.g` and `SRR16230426_R2.fastq.gz` are for sample `SRR16230426`.
Input files must be in FASTQ format, with ".fastq.gz" extensions supported.

### resource

`46bg.fasta` and `49btm.fasta` are vector sequence files as specified by the `vector` column in config/samples.tsv

`osa.fasta` is the host reference FASTA file as specified by the `host` column in config/samples.tsv

### config
* `config.yml` contains processing parameters:
    + `samples: samples.tsv` specifies the name of the sample information file in the same directory, e.g., `samples.tsv`. 
    + `max_filler_len: 10` maximum length of small DNA fragments insertion at T-DNA integration sites.
    + `max_microhomology_len: 5` maximum number of microhomology at the integration site.
    + `min_local_identity: 0.9` is the minimum percent identity of matched bases in aligned regions.
    + `min_reads_num: 3` minimum number of supporting informative split reads.
    + `cluster_dist: 5` maximum distance used to cluster sites into the same integration site.
    + `bwa:` 
        - `threads: 8` number of BWA threads.

* `samples.tsv` contains sample metadata:
    ```tsv
    sample  host    vector
    SRR16230416 osa.fasta   46bg.fasta
    SRR16230426 osa.fasta   49btm.fasta
    ```

### Running PTIS
PTIS is built as a Docker container image and shared through Docker Hub (https://hub.docker.com/repository/docker/caibinperl/ptis/general). 
The container can run on various operating systems, such as Linux, Mac, and Windows, all supported by the Docker engine. 

Before running the application, specify the parameter and place all the read files in the data directory, as mentioned above. For example, the data directory is `/home/bc/data`.

```sh
docker run --mount type=bind,src=/home/bc/data,target=/ptis/data caibinperl/ptis snakemake --cores 4
```

The --mount option tells Docker to create a bind mount, where `/home/bc/data` is the directory on your host machine, and `/ptis/data` would appear inside the container.
The `--cores 4` option is used to set threads.
