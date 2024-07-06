import util

read_names_file = snakemake.input[0]
fq1_file = snakemake.input[1]
fq2_file = snakemake.input[2]
anchored_fq1_file = snakemake.output[0]
anchored_fq2_file = snakemake.output[1]

util.filter_fq(read_names_file, fq1_file, anchored_fq1_file)
util.filter_fq(read_names_file, fq2_file, anchored_fq2_file)
