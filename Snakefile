# reference genome
REF = "data/MN908947.3.fasta"
# extract values of wildcard -- forward and reverse direction of the files 
(READS,) = glob_wildcards("data/sample.{read}.paired.fq.gz")

# target rule that tells snakemake what final files to produce 
rule all:
	input:
		"QC/fastqc/",
		"QC/multiqc/"

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