import util

bam1_file = snakemake.input[0]
bam2_file = snakemake.input[1]
aln1_file = snakemake.output[0]
aln2_file = snakemake.output[1]

util.load_aln(bam1_file, aln1_file)
util.load_aln(bam2_file, aln2_file)
