"""
workflow 03_Host_removal.smk
Author: Patricia Agudelo-Romero
email : Patricia.AgudeloRomero@telethonkids.org.au
"""

import os
DIR = os.getcwd()

#configfile: "config/config.yaml"
#SAMPLES, = glob_wildcards(os.path.join(config["input_DIR"],"{sample}_R1.fastq.gz"))

#rule all:
#	input:
#		os.path.join(config["output_DIR"],"EVEREST/multiQC_rep/fastq_unmapped_multiqc_report.html"),

rule KALLISTO_index:
	input:
		fasta = config["transcriptome"],
	output:
		idx = os.path.join(config["output_DIR"],"index_kallisto/Homo_sapiens.GRCh38.idx"),
	log:
		os.path.join(config["output_DIR"],"EVEREST/logs/03_01_KALLISTO_index.log"),
	benchmark:
		os.path.join(config["output_DIR"],"EVEREST/benchmarks/03_01_KALLISTO_index.txt")
	conda:
		os.path.join(DIR, "envs/kallisto.yml"),
	message:
		"star index",
	shell:
		(" kallisto index -i {output.idx} {input.fasta} ")

rule KALLISTO_align:
	input:
		fq1 = os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_trimm_R1.fastq.gz"),
		idx = os.path.join(config["output_DIR"],"index_kallisto/Homo_sapiens.GRCh38.idx"),
	output:
		bam = os.path.join(config["output_DIR"],"EVEREST/KALLISTO/{sample}/pseudoalignments.bam"),
		dir = directory(os.path.join(config["output_DIR"],"EVEREST/KALLISTO/{sample}")),
	params:
		threads = "5",
	log:
		os.path.join(config["output_DIR"],"EVEREST/logs/03_02_KALLISTO_align_{sample}.log"),
	benchmark:
		os.path.join(config["output_DIR"],"EVEREST/benchmarks/03_02_KALLISTO_index_{sample}.txt")
	conda:
		os.path.join(DIR, "envs/kallisto.yml"),
	message:
		"kallisto alingment",
	shell:
		(" kallisto quant -i {input.idx} \
		-o {output.dir} --single -l 200 -s 20 {input.fq1} \
		-t {params.threads} --pseudobam 2> {log} ")

rule SAMTOOLS_fastq:
	input:
		bam = os.path.join(config["output_DIR"],"EVEREST/KALLISTO/{sample}/pseudoalignments.bam"),
	output:
		o1 = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_R1.fastq"),
	params:
		threads = "7",
		type_u = "view -f 4 -h",
		type_s = "sort -@ 7",
		type_f = "fastq -NO -@ 7",
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R03_S02_SAMTOOLS_fastq_{sample}.log"),
	benchmark: 
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R03_S02_SAMTOOLS_fastq_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/minimap2.yml"),
	message:
		"obtain unmapped from kallisto and convert to fastq",
	shell:
		(" ( samtools {params.type_u} {input.bam} | samtools {params.type_s} \
		| samtools {params.type_f} - > {output.o1} ) &> {log} ")  

rule PIGZ_fastq:
	input:
		os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_R1.fastq")
	output:
		os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_R1.fastq.gz")
	params:
		"-p 7 -5",
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R03_S05_PIGZ_fastq_{sample}.log"),
	conda:
		os.path.join(DIR, "envs/minimap2.yml"),
	message:
		"compress fastq files",
	shell:
		(" pigz {params} {input} ")

rule BBMAP_dup:
	input:
		f1 = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_R1.fastq.gz"),
	output:
		f  = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_dedup.fastq.gz"),
	params:
		mem = "-Xmx20000m",
		other = "ac=f s=5 e=5 minidentity=95",
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R03_S06_BBMAP_dup_{sample}.log"),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R03_S06_BBMAP_dup_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/BBMAP.yml"),
	message:
		"BBMAP deduplication",
	shell:
		(" dedupe.sh {params.mem} {params.other} in={input.f1} out={output.f} 2> {log} ")

rule BBMAP_duduped_normalisation:
	input:
		f1  = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_dedup.fastq.gz"),
	output:
		f1  = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_dedup_norm_R1.fastq.gz"),
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R03_S08_BBMAP_dup_norm_{sample}.log"),
	benchmark: 
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R03_S08_BBMAP_dup_norm_{sample}.txt"),
	params:
		mem = "-Xmx20000m",
	conda:
		os.path.join(DIR, "envs/BBMAP.yml"),
	message:
		"BBMAP digital normalisation",
	shell:
		("  bbnorm.sh {params.mem} in={input.f1} out={output.f1} 2> {log} ")

rule FASTQC_unmapped:
	input:
		expand([os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_dedup_norm_R1.fastq.gz")], sample=SAMPLES),
	output:
		html = expand([os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_dedup_norm_R1_fastqc.html")], sample=SAMPLES),
	params:
		threads = "7",
		out_dir = os.path.join(config["output_DIR"], "EVEREST/FASTQ"),
	log:
		expand([os.path.join(config["output_DIR"], "EVEREST/logs/R03_S09_FASTQC_before_merge_{sample}.log")], sample=SAMPLES),
	conda:
		os.path.join(DIR, "envs/QC.yml"),
	message:
		"FASTQC - fastq files before merge",
	shell:
		(" fastqc {input} -t {params.threads} --outdir {params.out_dir}  2> {log} ")
	
rule multiQC_unmapped:
	input:
		expand([os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_dedup_norm_R1_fastqc.html")], sample=SAMPLES),	
	output:
		os.path.join(config["output_DIR"],"EVEREST/multiQC_rep/fastq_unmapped_multiqc_report.html"),
	params:
		in_dir = os.path.join(config["output_DIR"],"EVEREST/FASTQ"),
	log:
		os.path.join(config["output_DIR"],"EVEREST/logs/R03_S10_multiQC_before_merge.log"),
	conda:
		os.path.join(DIR, "envs/QC.yml"),
	message:
		"MultiQC report - fastq files before merge",
	shell:
		(" multiqc {params.in_dir} -n {output} -i fastq_before_merge 2> {log} ")
