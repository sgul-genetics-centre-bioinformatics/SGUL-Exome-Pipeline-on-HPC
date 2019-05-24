#!/usr/bin/env python

#SGUL GENETICS CENTRE EXOME PIPELINE APRIL 2019
#Alan Pittman & Dionysios Grigoriadis

#Dependencies for the pipeline

#DRIVE
DRIVE = "/storage/root/homes/apittman/resources/"

#resources: Global
BWAsoftware = DRIVE + "bwa/bwa"
BWAindex = DRIVE + "Genome_reference_files/human_g1k_v37.fasta"
samtoolsSoftware = DRIVE + "samtools-1.8/samtools"
java = DRIVE + "java/jre1.8.0_171/bin/java"
picard = DRIVE + "picard-2.815/picard.jar"
gatk = DRIVE + "gatk-4.0.4.0/gatk-package-4.0.4.0-local.jar"
gatk4_1 = DRIVE + "gatk-4.1.1.0/gatk-package-4.1.1.0-local.jar"
refknownsitesSNPS = DRIVE + "Genome_reference_files/common_all.vcf"
refknownsitesINDELS = DRIVE + "Genome_reference_files/Mills_and_1000G_gold_standard.indels.hg19_modified.vcf"

convert2annovar = DRIVE + "annovar/convert2annovar.pl"
tableAnnovar = DRIVE + "annovar/table_annovar.pl"
humandb = DRIVE + "annovar/humandb"

#resources: FastQC and Trimming
fastqc = DRIVE 
multiqc = DRIVE + "multiqc"
Trimmomatic= DRIVE + "Trimmomatic-0.38/trimmomatic-0.38.jar"

#resources: Data Pre-processing
###############################################################################################

externalBAMdir = DRIVE + "/Exomes/UBAMs"

#resources: Variant Calling
###############################################################################################
ExomeTarget = DRIVE + "Genome_reference_files/BroadExACExomeIntervlas.bed"


