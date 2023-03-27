"""
workflow 05_Cleaning_contigs.smk
Author: Patricia Agudelo-Romero
email : Patricia.AgudeloRomero@telethonkids.org.au
"""

import os
DIR = os.getcwd() 

configfile: "config/config.yaml"
SAMPLES, = glob_wildcards(os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fna"))
#SAMPLES, = glob_wildcards(os.path.join(config["input_DIR"],"{sample}_R1.fastq.gz"))


rule all:
	input:
		expand(os.path.join(config["output_DIR"], "EVEREST/BBMAP_stats/{sample}_contig_rpkm.txt"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/BBMAP_stats/{sample}_covstats.txt"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/BBMAP_stats/{sample}_scafstats.txt"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fasta.bacphlip"), sample = SAMPLES),

rule BBMAP_mapping_contigs:
	input:
		contigs  = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fna"),
		R1 = os.path.join(config["input_DIR"], "{sample}_R1.fastq.gz"),
		R2 = os.path.join(config["input_DIR"], "{sample}_R2.fastq.gz"),
	output:
		sam = os.path.join(config["output_DIR"], "EVEREST/BBMAP_stats/{sample}_contig.sam"),
		rpkm = os.path.join(config["output_DIR"], "EVEREST/BBMAP_stats/{sample}_contig_rpkm.txt"),
		covstats = os.path.join(config["output_DIR"], "EVEREST/BBMAP_stats/{sample}_covstats.txt"), 
		scafstats = os.path.join(config["output_DIR"], "EVEREST/BBMAP_stats/{sample}_scafstats.txt"),
	params:
		mem = "-Xmx20000m",
		extras = "nodisk slow=t ambiguous=random threads=7"
	log:
		temp(os.path.join(config["output_DIR"], "EVEREST/logs/R05_S06_BBMAP_mapping_contigs_{sample}.log")),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R05_S06_BBMAP_mapping_contigs_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/BBMAP.yml"),
	message:
		"Mapping contigs",
	shell:
		(" bbmap.sh {params.mem} ref={input.contigs} in={input.R1} in2={input.R2} out={output.sam} rpkm={output.rpkm} scafstats={output.scafstats} covstats={output.covstats} {params.extras}  2> {log} ")

rule BACPHLIP_life_style:
	input:
		fasta = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fna"),
	output:
		txt = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fasta.bacphlip"),
#		hmm = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fasta.hmmsearch.tsv"),
	log:
		temp(os.path.join(config["output_DIR"], "EVEREST/logs/R05_S08_BACPHLIP_life_style_{sample}.log")),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R05_S08_BACPHLIP_life_style_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/bacphlip.yml"),
	message:
		"life style"
	shell:
		(" bacphlip -i {input.fasta} -f --multi_fasta ")

