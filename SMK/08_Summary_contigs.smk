"""
workflow 08_Summary_contigs.smk
Author: Patricia Agudelo-Romero
email : Patricia.AgudeloRomero@telethonkids.org.au
"""

import os
DIR = os.getcwd() 

#configfile: "config/config.yaml"
#SAMPLES, = glob_wildcards(os.path.join(config["input_DIR"],"{sample}_R1.fastq.gz"))

#rule all:
#	input:
#		expand(os.path.join(config["output_DIR"], "EVEREST/Summary/{sample}_nt_summary_mmseqs2.txt"), sample=SAMPLES),


rule BBMAP_mapping_contigs:
	input:
		contigs  = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fna"),
		R1 = os.path.join(config["input_DIR"], "{sample}_R1.fastq.gz"),
		R2 = os.path.join(config["input_DIR"], "{sample}_R2.fastq.gz"),
	output:
		sam = os.path.join(config["output_DIR"], "EVEREST/BBMAP_stats/{sample}_contig.sam"),
	params:
		mem = "-Xmx20000m",
		extras = "nodisk slow=t ambiguous=random threads=7"
	log:
		temp(os.path.join(config["output_DIR"], "EVEREST/logs/R08_S01_BBMAP_mapping_contigs_{sample}.log")),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R08_S01_BBMAP_mapping_contigs_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/BBMAP.yml"),
	message:
		"Mapping contigs",
	shell:
		(" bbmap.sh {params.mem} ref={input.contigs} in={input.R1} in2={input.R2} out={output.sam} {params.extras}  2> {log} ")

rule BBMAP_pileup_summary:
	input:
		contigs  = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fna"),
		sam = os.path.join(config["output_DIR"], "EVEREST/BBMAP_stats/{sample}_contig.sam"),
	output:
		rpkm = os.path.join(config["output_DIR"], "EVEREST/BBMAP_stats/{sample}_rpkm.txt"),
		stats = os.path.join(config["output_DIR"], "EVEREST/BBMAP_stats/{sample}_stats.txt"),
		normcov = os.path.join(config["output_DIR"], "EVEREST/BBMAP_stats/{sample}_normcov.txt"),
	params:
		mem = "-Xmx20000m",
		extras = "binsize=500 header=t stdev=t countgc=t threads=7",
	log:
		temp(os.path.join(config["output_DIR"], "EVEREST/logs/R08_S02_BBMAP_pileup_summary_{sample}.log")),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R08_S02_BBMAP_pileup_summary_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/BBMAP.yml"),
	message:
		"BBMAP stats contigs"
	shell:
		(" bbmap.sh {params.mem} ref={input.contigs} in={input.sam} {params.extras} \
		rpkm={out.rpkm} out={out.stats} normcov={out.normcov}  2> {log} ")

