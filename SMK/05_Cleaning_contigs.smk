"""
workflow 05_Cleaning_contigs.smk
Author: Patricia Agudelo-Romero
email : Patricia.AgudeloRomero@telethonkids.org.au
"""

import os
DIR = os.getcwd() 

configfile: "config/config.yaml"
SAMPLES, = glob_wildcards(os.path.join(config["input_DIR"],"{sample}_R1.fastq.gz"))

rule all:
	input:
		expand(os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_rep_seq_FilterLen.fasta"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/VIRSORTER/{sample}/final-viral-combined.fa"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/VIRSORTER/{sample}/final-viral-score.tsv"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/VIRSORTER/{sample}/final-viral-boundary.tsv"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses.fna"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fasta"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/BBMAP_stats/{sample}_contig_rpkm.txt"), sample = SAMPLES),
#		expand(os.path.join(config["output_DIR"], "EVEREST/BBMAP_stats/{sample}_stats.txt"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/BBMAP_stats/{sample}_covstats.txt"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/BBMAP_stats/{sample}_scafstats.txt"), sample = SAMPLES),
#abricate
		expand(os.path.join(config["output_DIR"], "EVEREST/ABRICATE/{sample}_argannot.txt"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/ABRICATE/{sample}_ecoh.txt"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/ABRICATE/{sample}_ecoli_vf.txt"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/ABRICATE/{sample}_megares.txt"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/ABRICATE/{sample}_ncbi.txt"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/ABRICATE/{sample}_plasmidfinder.txt"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/ABRICATE/{sample}_resfinder.txt"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/ABRICATE/{sample}_vfdb.txt"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fasta.bacphlip"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fasta.hmmsearch.tsv"), sample = SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/PROKKA/{sample}/fastas"), sample = SAMPLES),
		os.path.join(config["output_DIR"], "EVEREST/PROKKA/database.txt"),
		expand(os.path.join(config["output_DIR"], "EVEREST/PROKKA/{sample}/annotation"), sample = SAMPLES),

rule PROKKA_db:
	input:
		db = config["Prokka"],
	output:
		file = os.path.join(config["output_DIR"], "EVEREST/PROKKA/database.txt"),
	log:
		temp(os.path.join(config["output_DIR"], "EVEREST/logs/R05_S10_PROKKA_db.log")),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R05_S10_PROKKA_db.txt"),
	conda:
		os.path.join(DIR, "envs/prokka.yml"),
	message:
		"Prokka PHROGS database",
	shell:

if [ -f $( cat text.txt )"/hmm/all_phrogs.hmm" ]; then      echo "It exists"; fi

		(" echo $(prokka --listd 2>&1 | sed 1q | sed 's:.*\: ::') > {output.file} ; \
		if [ ! -f $(cat {output.file})"/hmm/all_phrogs.hmm" ] \
		then \
			cp {input.db} $(prokka --listd 2>&1 | sed 1q | sed 's:.*\: ::') ; \
			prokka --setupdb 2> {log} ; \
			sed -i '1 i Database ready' {output.file} \
		else \
			print("PHROGS database already copied!")
		fi ")

#path="$(prokka --listd 2>&1 | sed 1q | sed 's:.*\: ::')" ; \
#		(" cp {input.db} $(prokka --listd 2>&1 | sed 1q | sed 's:.*\: ::') ; \
#		prokka --setupdb ; echo "PHROGS added" > {output.file} ")

rule SEQKIT_filter:
	input:
		fasta  = os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_rep_seq.fasta"),
	output:
		filter = os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_rep_seq_FilterLen.fasta"),
	params:
		"seq -m 5000",
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R05_S02_SEQKIT_filter_{sample}.log"),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R05_S02_SEQKIT_filter_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/seqkit.yml"),
	message:
		"remove contigs < 5000 bp",
	shell:
		(" seqkit {params} {input.fasta} -o {output.filter} 2> {log} ")

checkpoint VIRSORTER_detect:
	input:
		fasta = os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_rep_seq_FilterLen.fasta"),
		db = config["VIRSORTER_DB"],
	output:
		dir = directory(os.path.join(config["output_DIR"], "EVEREST/VIRSORTER/{sample}")),
		viral = os.path.join(config["output_DIR"], "EVEREST/VIRSORTER/{sample}/final-viral-combined.fa"),
		score = os.path.join(config["output_DIR"], "EVEREST/VIRSORTER/{sample}/final-viral-score.tsv"),
		boundary = os.path.join(config["output_DIR"], "EVEREST/VIRSORTER/{sample}/final-viral-boundary.tsv"),
	params:
		type = "run --keep-original-seq",
		groups = "--include-groups dsDNAphage,NCLDV,RNA,ssDNA",
		others = "--min-length 5000 --min-score 0.5",
		threads = "-j 7 all",
#		dir  = os.path.join(config["output_DIR"], "EVEREST/VIRSORTER/{sample}/"),
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R05_S03_VIRSORTER_detect_{sample}.log"),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R05_S03_VIRSORTER_detect_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/virsorter2.yml"),
	message:
		"rid off from non-viral contigs",
	shell:
		(" virsorter {params.type} -i {input.fasta} -w {output.dir} --db-dir {input.db} {params.groups} {params.others} {params.threads} 2> {log} ")

def get_VIRSORTER_files(wildcards):
	checkpoint_output = checkpoints.VIRSORTER_detect.get(**wildcards).output[0]
	return expand(os.path.join(checkpoint_output,"final-viral-combined.fa"))

checkpoint CHECKV_viral_seq:
	input:
		fasta = get_VIRSORTER_files,
		DB = config["CHECKV_DB"],
#		fasta = os.path.join(config["output_DIR"], "EVEREST/VIRSORTER/{sample}/final-viral-combined.fa"),
	output:
		dir  = directory(os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}")),
		fasta = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses.fna"),
	params:
		type = "end_to_end",
		threads = "-t 7",
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R05_S04_CHECKV_viral_seq_{sample}.log")
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R05_S04_CHECKV_viral_seq_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/checkv.yml"),
	message:
		"rid off non-viral reads",
	shell:
		(" checkv {params.type} -d {input.DB} {input.fasta} {output.dir} 2> {log} ")

def get_CHECKV_files(wildcards):
	checkpoint_output = checkpoints.CHECKV_viral_seq.get(**wildcards).output[0]
	return expand(os.path.join(checkpoint_output,"viruses.fna"))

rule RENAME_viral_seq:
	input:
		get_CHECKV_files,
#		fasta = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses.fna"),
	output:
		fasta = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fasta"),
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R05_S05_RENAME_viral_seq_{sample}.log"),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R05_S05_RENAME_viral_seq_{sample}.txt"),
	message:
		"rename viral contigs",
	shell:
		(" sed 's/||.*//' {input} > {output.fasta} 2> {log} ")

rule BBMAP_mapping_contigs:
	input:
		contigs  = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fasta"),
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

rule ABRICATE_db:
	input:
		contigs = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fasta"), 
	output:
		argannot = os.path.join(config["output_DIR"], "EVEREST/ABRICATE/{sample}_argannot.txt"),
		card = os.path.join(config["output_DIR"], "EVEREST/ABRICATE/{sample}_card.txt"),
		ecoh = os.path.join(config["output_DIR"], "EVEREST/ABRICATE/{sample}_ecoh.txt"),
		ecoli_vf = os.path.join(config["output_DIR"], "EVEREST/ABRICATE/{sample}_ecoli_vf.txt"),
		megares = os.path.join(config["output_DIR"], "EVEREST/ABRICATE/{sample}_megares.txt"), 
		ncbi = os.path.join(config["output_DIR"], "EVEREST/ABRICATE/{sample}_ncbi.txt"),
		plasmidfinder = os.path.join(config["output_DIR"], "EVEREST/ABRICATE/{sample}_plasmidfinder.txt"),
		resfinder = os.path.join(config["output_DIR"], "EVEREST/ABRICATE/{sample}_resfinder.txt"),
		vfdb = os.path.join(config["output_DIR"], "EVEREST/ABRICATE/{sample}_vfdb.txt"),
	log:
		temp(os.path.join(config["output_DIR"], "EVEREST/logs/R05_S07_ABRICATE_BD_{sample}.log")),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R05_S07_ABRICATE_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/abricate.yml"),
	message:
		"contigs' screening for antimicrobial resistance or virulence genes"
	shell:
		(" abricate --db argannot --quiet {input.contigs} --nopath > {output.argannot} ; \
		abricate --db card --quiet {input.contigs} > {output.card} ; \
		abricate --db ecoh --quiet {input.contigs} > {output.ecoh} ; \
		abricate --db ecoli_vf --quiet {input.contigs} > {output.ecoli_vf} ; \
		abricate --db megares --quiet {input.contigs} > {output.megares} ; \
		abricate --db ncbi --quiet {input.contigs} > {output.ncbi} ; \
		abricate --db plasmidfinder --quiet {input.contigs} > {output.plasmidfinder} ; \
		abricate --db resfinder --quiet {input.contigs} > {output.resfinder} ; \
		abricate --db vfdb --quiet {input.contigs} > {output.vfdb} ")

rule BACPHLIP_life_style:
	input:
		fasta = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fasta"),
	output:
		txt = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fasta.bacphlip"),
		hmm = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fasta.hmmsearch.tsv"),
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

rule SEQKIT_split_multifasta:
	input:
		fasta = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fasta"),
	output:
		dir = directory(os.path.join(config["output_DIR"], "EVEREST/PROKKA/{sample}/fastas")),
	log:
		temp(os.path.join(config["output_DIR"], "EVEREST/logs/R05_S09_SEQKIT_split_multifasta_{sample}.log")),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R05_S09_SEQKIT_split_multifasta_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/seqkit.yml"),
	message:
		"splitting fasta for Prokka annotation",
	shell:
		(" seqkit split2 --by-size 1 {input.fasta} -O {output.dir} ; \
		for f in {output.dir}/*.fasta ; do mv $f {output.dir}/$( seqkit seq --name --only-id $f).fasta ; done ")

rule PROKKA_annotation:
	input:
		dir = os.path.join(config["output_DIR"], "EVEREST/PROKKA/{sample}/fastas"),
	output:
		dir = directory(os.path.join(config["output_DIR"], "EVEREST/PROKKA/{sample}/annotation")),
	params:
		name = "{sample}"
	log:
		temp(os.path.join(config["output_DIR"], "EVEREST/logs/R05_S11_PROKKA_annotation_{sample}.log")),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R05_S11_PROKKA_annotation_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/prokka.yml"),
	message:
		"Prokka annotation",
	shell:
		(" for fasta in {input.dir}/*.fasta ; do prokka --outdir {output.dir} --force --locustag {params.name} --addgenes --addmrna --metagenome --rfam --kingdom viruses --rfam $fasta ; done ")



## Annotation
#        prokka \
#		{input.fasta} \
#		--locustag ${base} \
#		--hmms ${PAPH}/all_phrogs.hmm \
#		--outdir ${OUT} \
#		--kingdom viruses

#sed '/^>/ s/ .*//' ${TEMP}/prokka.faa >${TEMP}/network.faa

#cat ${TEMP}/*.faa >${TEMP}/combined.faa
# R -f ${RUN}/combine-networks.r

#    vcontact2 \
#	--raw-proteins ${FIN}/combined.faa \
#	--proteins-fp ${FIN}/combined.csv \
#	--db 'ProkaryoticViralRefSeq211-Merged' \
#	--output-dir ${FIN}/results

#    # MAFFT
#    mafft --auto ${FIN}/combined.fasta >${FIN}/mafft_alignment.fasta
#	echo "NUCL substitution model"
#	raxmlHPC \
#	-f a \
#	-p $RANDOM \
#	-x $RANDOM \
#	-N 100 \
#	-m GTRGAMMAI \
#	-s ${FIN}/mafft_alignment.fa \
#	-n nucl_tree \
#	-w ${FIN}
#	echo "Finished tree_drawing"

#	echo "AA substitution model"
#	raxmlHPC \
#	-f a \
#	-p $RANDOM \
#	-x $RANDOM \
#	-N 100 \
#	-m PROTGAMMAI \
#	-s ${FIN}/mafft_alignment.fasta \
#	-n prot_tree \
#	-w ${FIN}



