import os

import pandas as pd


configfile: "config/config.yml"


def prepare_dir():
    result_dir = os.path.join("results")
    temp_dir = os.path.join(result_dir,"temp")
    trim_dir = os.path.join(temp_dir,"trim_reads")
    aln_dir = os.path.join(temp_dir,"aln")
    load_aln_dir = os.path.join(temp_dir,"load_aln")
    filter_aln_dir = os.path.join(temp_dir,"filter_aln")
    filter_fq_dir = os.path.join(temp_dir,"filter_fq")
    align_vector_dir = os.path.join(temp_dir,"align_vector")
    load_vector_aln_dir = os.path.join(temp_dir,"load_vector_aln")
    align_host_dir = os.path.join(temp_dir,"align_host")
    load_host_aln_dir = os.path.join(temp_dir,"load_host_aln")
    bam_dir = os.path.join(temp_dir,"bam")

    if not os.path.exists(result_dir):
        os.mkdir(result_dir)
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
    if not os.path.exists(trim_dir):
        os.mkdir(trim_dir)
    if not os.path.exists(aln_dir):
        os.mkdir(aln_dir)
    if not os.path.exists(load_aln_dir):
        os.mkdir(load_aln_dir)
    if not os.path.exists(filter_aln_dir):
        os.mkdir(filter_aln_dir)
    if not os.path.exists(filter_fq_dir):
        os.mkdir(filter_fq_dir)
    if not os.path.exists(align_vector_dir):
        os.mkdir(align_vector_dir)
    if not os.path.exists(load_vector_aln_dir):
        os.mkdir(load_vector_aln_dir)
    if not os.path.exists(align_host_dir):
        os.mkdir(align_host_dir)
    if not os.path.exists(load_host_aln_dir):
        os.mkdir(load_host_aln_dir)
    if not os.path.exists(bam_dir):
        os.mkdir(bam_dir)


def get_vector(wildcards):
    return os.path.join("resources",samples["vector"][wildcards.sample])


def get_host(wildcards):
    return os.path.join("resources",samples["host"][wildcards.sample])


def index_ref(ref):
    indexed_files = [
        f"resources/{ref}.{ext}" for ext in ["amb", "ann", "bwt", "pac", "sa"]
    ]
    indexed = True
    for index_file in indexed_files:
        if not os.path.exists(index_file):
            indexed = False
            break

    if not indexed:
        os.system(f"bwa index resources/{ref}")


def start_on():
    prepare_dir()
    for sample in samples.index:
        index_ref(samples["host"][sample])
        index_ref(samples["vector"][sample])


samples = pd.read_table(os.path.join("config",config["samples"])).set_index(
    "sample",drop=False
)
start_on()


rule all:
    input:
        expand("results/{sample}_site.tsv",sample=samples.index),
        expand("results/bam/{sample}",sample=samples.index),


rule trim_reads:
    input:
        "fastq/{sample}_R1.fastq.gz",
        "fastq/{sample}_R2.fastq.gz",
    output:
        "results/temp/trim_reads/{sample}_R1.fastq.gz",
        "results/temp/trim_reads/{sample}_R2.fastq.gz",
        "results/temp/trim_reads/{sample}_fastp.html",
        "results/temp/trim_reads/{sample}_fastp.json",
    shell:
        "fastp --dont_eval_duplication -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} "
        "-h {output[2]} -j {output[3]}"


rule align:
    input:
        get_vector,
        "results/temp/trim_reads/{sample}_R1.fastq.gz",
        "results/temp/trim_reads/{sample}_R2.fastq.gz",
    output:
        "results/temp/aln/{sample}_R1.bam",
        "results/temp/aln/{sample}_R2.bam",
    threads: config["bwa"]["threads"]
    shell:
        "bwa mem -t {threads} {input[0]} {input[1]} |samtools view -b |"
        "samtools sort - -o {output[0]} && samtools index {output[0]} && "
        "bwa mem -t {threads} {input[0]} {input[2]} |samtools view -b |"
        "samtools sort - -o {output[1]} && samtools index {output[1]}"


rule load_aln:
    input:
        "results/temp/aln/{sample}_R1.bam",
        "results/temp/aln/{sample}_R2.bam",
    output:
        "results/temp/load_aln/{sample}_R1.aln",
        "results/temp/load_aln/{sample}_R2.aln",
    script:
        "scripts/load_aln.py"


rule filter_aln:
    input:
        "results/temp/load_aln/{sample}_R1.aln",
        "results/temp/load_aln/{sample}_R2.aln",
    output:
        "results/temp/filter_aln/{sample}_read_name.ids",
    script:
        "scripts/filter_vector_aln.R"


rule filter_fq:
    input:
        "results/temp/filter_aln/{sample}_read_name.ids",
        "results/temp/trim_reads/{sample}_R1.fastq.gz",
        "results/temp/trim_reads/{sample}_R2.fastq.gz",
    output:
        "results/temp/filter_fq/{sample}_R1.fastq",
        "results/temp/filter_fq/{sample}_R2.fastq",
    script:
        "scripts/filter_fq.py"


rule align_vector:
    input:
        get_vector,
        "results/temp/filter_fq/{sample}_R1.fastq",
        "results/temp/filter_fq/{sample}_R2.fastq",
    output:
        "results/temp/align_vector/{sample}_R1.sam",
        "results/temp/align_vector/{sample}_R2.sam",
    threads: config["bwa"]["threads"]
    shell:
        "bwa mem -t {threads} -o {output[0]} {input[0]} {input[1]} && "
        "bwa mem -t {threads} -o {output[1]} {input[0]} {input[2]}"


rule load_vector_aln:
    input:
        "results/temp/align_vector/{sample}_R1.sam",
        "results/temp/align_vector/{sample}_R2.sam",
    output:
        "results/temp/load_vector_aln/{sample}_R1.aln",
        "results/temp/load_vector_aln/{sample}_R2.aln",
    script:
        "scripts/load_aln.py"


rule align_host:
    input:
        get_host,
        "results/temp/filter_fq/{sample}_R1.fastq",
        "results/temp/filter_fq/{sample}_R2.fastq",
    output:
        "results/temp/align_host/{sample}_R1.sam",
        "results/temp/align_host/{sample}_R2.sam",
    threads: config["bwa"]["threads"]
    shell:
        "bwa mem -t {threads} -o {output[0]} {input[0]} {input[1]} && "
        "bwa mem -t {threads} -o {output[1]} {input[0]} {input[2]}"


rule load_host_aln:
    input:
        "results/temp/align_host/{sample}_R1.sam",
        "results/temp/align_host/{sample}_R2.sam",
    output:
        "results/temp/load_host_aln/{sample}_R1.aln",
        "results/temp/load_host_aln/{sample}_R2.aln",
    script:
        "scripts/load_aln.py"


rule filter_site_aln:
    input:
        "results/temp/load_vector_aln/{sample}_R1.aln",
        "results/temp/load_vector_aln/{sample}_R2.aln",
        "results/temp/load_host_aln/{sample}_R1.aln",
        "results/temp/load_host_aln/{sample}_R2.aln",
    output:
        "results/{sample}_site_aln.tsv",
    script:
        "scripts/filter_site_aln.R"


rule detect_site:
    input:
        "results/{sample}_site_aln.tsv",
    output:
        "results/{sample}_site.tsv",
    script:
        "scripts/detect_site.py"


rule filter_sam:
    input:
        "results/{sample}_site.tsv",
        "results/{sample}_site_aln.tsv",
        "results/temp/align_vector/{sample}_R1.sam",
        "results/temp/align_vector/{sample}_R2.sam",
        "results/temp/align_host/{sample}_R1.sam",
        "results/temp/align_host/{sample}_R2.sam",
    output:
        directory("results/bam/{sample}"),
    script:
        "scripts/filter_sam.py"
