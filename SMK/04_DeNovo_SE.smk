"""
workflow 04_DeNovo.smk
Author: Patricia Agudelo-Romero
email : Patricia.AgudeloRomero@telethonkids.org.au
"""

import os
DIR = os.getcwd()

#configfile: "config/config.yaml"
#SAMPLES, = glob_wildcards(os.path.join(config["input_DIR"],"{sample}_R1.fastq.gz"))

#rule all:
#	input:
#		expand(os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_all_seqs.fasta"), sample=SAMPLES),
#		expand(os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_cluster.tsv"), sample=SAMPLES),
#		expand(os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_rep_seq.fasta"), sample=SAMPLES),

checkpoint SPADES_DeNovo:
	input:
		f  = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_dedup_norm_R1.fastq.gz"),
	output:
		fasta = os.path.join(config["output_DIR"], "EVEREST/SPADES/{sample}/scaffolds.fasta"),
		dir = directory(os.path.join(config["output_DIR"], "EVEREST/SPADES/{sample}")),
	params:
		threads = "-t 10",
		extra = "--only-assembler",
		type = "--careful",
#		type = "--meta",
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R04_S04_SPADES_DeNovo_{sample}.log"),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R04_S04_SPADES_DeNovo_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/spades.yml"),
	message:
		"De novo assembly",
	shell:
		(" spades.py {params.type} \
		--s1 {input.f} \
		-o {output.dir} \
		{params.threads} {params.extra} &> {log} ")
# this doesn't work in --meta		--s1 {input.u1} --s2 {input.u2} \

def get_SPADES_files(wildcards):
	checkpoint_output = checkpoints.SPADES_DeNovo.get(**wildcards).output[1]
	print(checkpoint_output)
	return expand(os.path.join(checkpoint_output,"scaffolds.fasta"))

rule MMSEQ2_eLinclust:
	input:
		get_SPADES_files,
#       fasta = os.path.join(config["output_DIR"], "EVEREST/SPADES/{sample}/scaffolds.fasta"),
	output:
		all  = os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_all_seqs.fasta"),
		clus = os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_cluster.tsv"),
		rep  = os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_rep_seq.fasta"),
	params:
		prefix  = os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}"),
		tmp_dir = os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_tmp"),
		threads = "7",
		type    = "--min-seq-id 0.98 --kmer-per-seq-scale 0.3 --sort-results 1 --alignment-mode 3 --cov-mode 1",
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R05_S01_MMSEQ2_eLinclust_{sample}.log"),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R05_S01_MMSEQ2_eLinclust_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/MMSEQS.yml"),
	message:
		"cluster sequences"
	shell:
		(" mmseqs easy-linclust {input} {params.prefix} {params.tmp_dir} --threads {params.threads} {params.type} 2> {log} ")

