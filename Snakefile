# reference genome
REF = "data/MN908947.3.fasta"
# extract values of wildcard -- forward and reverse direction of the files 
(READS,) = glob_wildcards("data/sample.{read}.paired.fq.gz")

# target rule that tells snakemake what final files to produce 
rule all:
	input:
		"QC/fastqc/",
		"QC/multiqc/",
		"results/alinged_reads/sample.sorted.bam"

rule fastqc:
	input:
		expand("data/sample.{read}.paired.fq.gz", read= READS)
	output:
		"QC/fastqc/"
	log:
		"logs/fastqc.log"
	shell:
		"fastqc {input} -o {output} 2> {log}"

rule multiqc:
	input:
		"QC/fastqc/"
	output:
		"QC/multiqc/"
	log:
		"logs/multiqc.log"
	shell:
		"multiqc {input} --force -o {output} 2> {log}"

rule bwa_aling:
	input:
		ref = REF,
		r1= "data/sample.R1.paired.fq.gz",
		r2 = "data/sample.R2.paired.fq.gz"
	output: 
		"results/alinged_reads/sample.bam"
	log:
		"logs/bwa_aling.log"
	shell: 
		"(bwa index {input.ref} "
		"&& bwa mem {input.ref} {input.r1} {input.r2} | samtools view -Sbh > {output}) 2> {log}"

rule sort:
	input:
		"results/alinged_reads/sample.bam"
	output:
		"results/alinged_reads/sample.sorted.bam"
	log:
		"logs/sort.log"
	shell:
		"samtools sort {input} -o {output} 2> {log}"
