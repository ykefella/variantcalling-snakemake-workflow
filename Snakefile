# reference genome
REF = "data/MN908947.3.fasta"
# extract values of wildcard -- forward and reverse direction of the files 
(READS,) = glob_wildcards("data/sample.{read}.paired.fq.gz")

# target rule that tells snakemake what final files to produce 
rule all:
	input:
		"QC/fastqc/",
		"QC/multiqc/",
		"results/aligned_reads/dup.reads.bam",
		"results/aligned_reads/metrics.txt",
		"results/aligned_reads/dup.reads.bam.bai",
		"results/variantcall/sample.variants.vcf"

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

rule bwa_align:
	input:
		ref = REF,
		r1= "data/sample.R1.paired.fq.gz",
		r2 = "data/sample.R2.paired.fq.gz"
	output: 
		"results/aligned_reads/sample.bam"
	log:
		"logs/bwa_align.log"
	shell: 
		"(bwa index {input.ref} "
		"&& bwa mem {input.ref} {input.r1} {input.r2} | samtools view -Sbh > {output}) 2> {log}"

rule sort:
	input:
		"results/aligned_reads/sample.bam"
	output:
		"results/aligned_reads/sample.sorted.bam"
	log:
		"logs/sort.log"
	shell:
		"samtools sort {input} -o {output} 2> {log}"

rule markDuplication:
	input:
		"results/aligned_reads/sample.sorted.bam"
	output:
		dup ="results/aligned_reads/dup.reads.bam",
		met ="results/aligned_reads/metrics.txt"
	log:
		"logs/markDuplication.log"
	shell:
		"(java -jar /usr/local/Cellar/picard-tools/2.25.7/libexec/picard.jar MarkDuplicates INPUT={input} OUTPUT={output.dup} "
		"METRICS_FILE={output.met}) 2> {log}"

rule index:
	input:
		"results/aligned_reads/dup.reads.bam"
	output:
		"results/aligned_reads/dup.reads.bam.bai"
	log:
		"logs/index.log"
	shell:
		"samtools index {input} 2> {log}"

rule freebayes:
	input:
		ref = REF,
		bam = "results/aligned_reads/dup.reads.bam"
	output:
		"results/variantcall/sample.variants.vcf"
	log:
		"logs/freebayes.log"
	shell:
		"freebayes -f {input.ref} {input.bam} > {output} 2> {log}"