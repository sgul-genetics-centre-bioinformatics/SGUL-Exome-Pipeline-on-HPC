#!/usr/bin/env python

#SGUL GENETICS CENTRE EXOME PIPELINE MAY 2019
#SGUL HPC
#Alan Pittman

#standard pipeline steps:
###############################################################################################

# run BWA command
# run samtools sam to bam conversion command
# run samtools bam sort command 
# run samtools index
# run picard mark PCR Duplicates
# run samtools index again
# run GATK BaseRecalibrator #### need to use big VCF file
# run GATK AnalyzeCovariates ### work in progress need some R libraries installed 
# run GATK ApplyBQSR
	
###############################################################################################

import os
import sys
import subprocess
import csv
import time

#software and resources:
###############################################################################################

from dependencies import BWAsoftware
from dependencies import BWAindex
from dependencies import samtoolsSoftware
from dependencies import java
from dependencies import picard
from dependencies import gatk
from dependencies import refknownsitesSNPS
from dependencies import refknownsitesINDELS
from dependencies import ExomeTarget
from dependencies import convert2annovar
from dependencies import tableAnnovar
from dependencies import humandb

##################################################################################################

def display(message):
    print(message)
	
def ScriptWriter(step):
	for command in step:
		jobScript.write(str(command) + " ")
		
	jobScript.write(" \n")
	jobScript.write(" ")
	jobScript.write(" \n")
	

myProject = sys.argv
del myProject[0]

myProject = str(myProject)
myProject = myProject.lstrip('[')
myProject = myProject.lstrip("'")
myProject = myProject.rstrip(']')
myProject = myProject.rstrip("'")    

display("Genetics Research Centre Alignment Pipeline on SGUL HPC\n")

print("your project is :")
print(myProject)

print("\n")

dirpath = os.getcwd()
projectDataDir=dirpath + "/" + "Unaligned" + "/" + myProject + "/"
print(os.listdir(projectDataDir))

#generate list of all the samples in the project directory 
mySamples = os.listdir(projectDataDir)

VCFOUTdir = "/storage/root/homes/apittman/Exomes/VCF/"
gVCFOUTdir = "/storage/root/homes/apittman/Exomes/gVCF/"
AlignedDir = "/storage/root/homes/apittman/Exomes/Aligned/"
Annodir = "/storage/root/homes/apittman/Exomes/Annotated/"
	
#for each sample in the project folder, identify the fastq files (R1 and R2) for BWA alignment
for sample in mySamples:

	sampleDirectory = projectDataDir + sample + "/"

	myinputfiles = os.listdir(sampleDirectory) # assuming we have two fastq files per sample
	
	inputFASTQ1 = myinputfiles[0]
	inputFASTQ2 = myinputfiles[1]
	
	print("your input fastq files are:")
	
	#L1_2.fq.gz # check naming format of fastq files 
	print inputFASTQ1
	print inputFASTQ2
		
	#make analysis output directory for alignment
	
	alignedOUTdir = "/storage/root/homes/apittman/Exomes/Aligned/"
	
	projectFolder = alignedOUTdir + myProject 
	makeDirectoryProject = ['mkdir', projectFolder] 
	subprocess.call(makeDirectoryProject)
	
	sampleFolder = alignedOUTdir + myProject + "/" + sample
	makeDirectory = ['mkdir', sampleFolder] 
	subprocess.call(makeDirectory)
	
	#make analysis output directorys (VCF and GVCF)
	
	projectFolder = VCFOUTdir + myProject 
	makeDirectoryVCFProject = ['mkdir', projectFolder] 
	subprocess.call(makeDirectoryVCFProject)
	
	sampleFolderVCF = VCFOUTdir + myProject + "/" + sample
	makeDirectoryVCFProjectSample = ['mkdir', sampleFolderVCF] 
	subprocess.call(makeDirectoryVCFProjectSample)
	
	projectFolder = gVCFOUTdir + myProject 
	makeDirectoryGVCFProject = ['mkdir', projectFolder] 
	subprocess.call(makeDirectoryGVCFProject)
	
	sampleFolderGVCF = gVCFOUTdir + myProject + "/" + sample
	makeDirectoryGVCFProjectSample = ['mkdir', sampleFolderGVCF] 
	subprocess.call(makeDirectoryGVCFProjectSample)
	
	#make analysis output directorys (Annotated)
	
	projectFolder = Annodir + myProject 
	makeDirectoryAnnoProject = ['mkdir', projectFolder] 
	subprocess.call(makeDirectoryAnnoProject)
	
	sampleFolderAnno = Annodir + myProject + "/" + sample
	makeDirectoryAnnoProjectSample = ['mkdir', sampleFolderAnno] 
	subprocess.call(makeDirectoryAnnoProjectSample)
	
	############# SAMPLE SPECIFIC VARIABLES ###################################

	outputSAM = sampleFolder + "/" + sample + ".sam" # output samfile
	samHeader = "'@RG\\tID:" + sample + "\\tSM:" + sample + "\\tLB:" + sample + "\\tPL:ILLUMINA'"
	path_inputFASTQ1 = sampleDirectory + inputFASTQ1
	path_inputFASTQ2 = sampleDirectory + inputFASTQ2
	outputBAM = sampleFolder + "/" + sample + ".bam" # output bam file
	outputSorstedBAM = sampleFolder + "/" + sample + "_sorted.bam" # output sorted bam file
	SorstedUniqueBAM = sampleFolder + "/" + sample + "_sorted_unique" + ".bam" #output sorted unique bam file
	picardIbam = "I=" + outputSorstedBAM
	picardObam = "O=" + SorstedUniqueBAM
	picardMfile = "M=" + sampleFolder + "/" + sample + "_marked_dup_metrics.txt"
	SorstedUniqueBAMindex = sampleFolder + "/" + sample + "_sorted_unique" + ".bam" + ".bai"
	outputSorstedBAMindex = sampleFolder + "/" + sample + "_sorted.bam.bai"	
	recal_data_table = sampleFolder + "/" + sample + "_recal_data.table"
	AnalyzeCovariates_pdf = sampleFolder + "/" + sample + "_AnalyzeCovariates.pdf"
	SorstedUniqueRecalibratedBAM = sampleFolder + "/" + sample + "_sorted_unique_recalibrated" + ".bam"
	bamout = AlignedDir + "/" + myProject + "/" + sample + "/" + sample + "bamout" + ".bam"
	sampleVCF = sampleFolderVCF + "/" + sample + "_raw.vcf"
	FilteredsampleVCF = sampleFolderVCF + "/" + sample + "_MetricFilters.vcf"
	sampleGVCF = sampleFolderGVCF + "/" + sample + ".g.vcf.gz"
	sampleFolderVCF = VCFOUTdir + myProject + "/" + sample
	FilteredsampleVCF = sampleFolderVCF + "/" + sample + "_MetricFilters.vcf"
	AVinput = sampleFolderAnno + "/" + sample + ".avinput"
	Annovaroutput = sampleFolderAnno + "/" +sample + ".annovar"
	
	########## Job submission at top of script ########################

	interpreter = "#!/bin/bash"
	nodesandcores = "#PBS -l nodes=1:ppn=10"
	walltime = "#PBS -l walltime=36:00:00"
	pbs = "#PBS -S /bin/bash"
	abe = "#PBS -m abe"
	email = "#PBS -M apittman@sgul.ac.uk"
	
	#-l walltime=100:00:00
	
	########## COMMAND LINE ARGUMENTS #################################
	
	BWAcommand = [BWAsoftware, 
	'mem', 
	BWAindex, 
	path_inputFASTQ1, 
	path_inputFASTQ2, 
	'-t', 
	'8', 
	'-R', samHeader,
	'-o', outputSAM]
	
	samtoolsViewCommand = [samtoolsSoftware, 
	'view', 
	'-Sb', 
	outputSAM, 
	'-o', outputBAM]
	
	samtoolsSortCommand = [samtoolsSoftware,
	'sort', 
	'-m',
	'5000000000', 
	outputBAM, 
	'-o', outputSorstedBAM]
	
	samtoolsINDEXCommand = [samtoolsSoftware,
	'index',
	outputSorstedBAM]
	
	picardMARKDUPLICATESCommand = [java, 
	'-jar', 
	picard, 
	'MarkDuplicates',
	picardIbam, 
	picardObam, 
	picardMfile]
	
	samtoolsINDEXCommand2 = [samtoolsSoftware, 
	'index', 
	SorstedUniqueBAM]
	
	gatkBaserecalibratorCommand = [java, 
	'-jar', 
	gatk, 
	'BaseRecalibrator', 
	'-I', SorstedUniqueBAM,
	'-R', BWAindex, 
	'--known-sites', refknownsitesSNPS, 
	'-O', recal_data_table]
	
	gatkAnalyseCovariatesCommand = [java, 
	'-jar', 
	gatk, 
	'AnalyzeCovariates', 
	'-bqsr', recal_data_table,
	'-plots' , AnalyzeCovariates_pdf]
	
	gatkApplyBQSRCommand = [java, 
	'-jar', 
	gatk, 
	'ApplyBQSR', 
	'-R', BWAindex, 
	'-I', SorstedUniqueBAM,
	'--bqsr-recal-file', recal_data_table, 
	'-O', SorstedUniqueRecalibratedBAM]
	
	gatkHaplotypeCaller_VCF_Command = [java, 
	'-jar', 
	gatk, 
	'HaplotypeCaller', 
	'-I', SorstedUniqueRecalibratedBAM,
	'-R', BWAindex, 
	'--intervals', ExomeTarget, 
	'-O', sampleVCF, 
	'-bamout', bamout]
	
	gatkVariantFiltration_VCF_Command = [java, 
	'-jar', 
	gatk, 
	'VariantFiltration', 
	'-V', sampleVCF, 
	'-O', FilteredsampleVCF,
	'-R', BWAindex, 
	'--genotype-filter-expression', '"GQ < 30.0"', '--genotype-filter-name',
	'"LowGQ"', '--filter-expression', '"QD < 1.5"', '--filter-name', '"LowQD"', 
	'--filter-expression', '"DP < 6"', '--filter-name', '"LowCoverage"',
	'--filter-expression', '"SOR > 10.0"', '--filter-name', '"StrandBias"']
	
	gatkHaplotypeCaller_gVCF_Command = [java, 
	'-jar', 
	gatk, 
	'HaplotypeCaller', 
	'-I', SorstedUniqueRecalibratedBAM,
	'-R', BWAindex, 
	'--intervals', ExomeTarget, 
	'-O', sampleGVCF, 
	'-ERC', 'GVCF']
	
	convert2annovarCommand = [convert2annovar, 
	'-format', 
	'vcf4', 
	FilteredsampleVCF, 
	'-allsample', '-withfreq', '-includeinfo', '-outfile', AVinput]
	
	annovarCommand = [tableAnnovar, AVinput, humandb, '-buildver', 'hg19', '-out', Annovaroutput, '-remove',
	'-protocol', 'refGene,ensGene,cytoBand,genomicSuperDups,1000g2015aug_all,esp6500siv2_ea,esp6500siv2_aa,esp6500siv2_all,exac03,kaviar_20150923,gnomad_exome,gnomad_genome,hrcr1,gme,avsnp150,dbnsfp33a,dbnsfp31a_interpro,spidex,revel,mcap,clinvar_20170905',
	'-operation', 'g,g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f',
	'-genericdbfile', 'hg19_spidex.txt',
	'-argument', '"-hgvs,-hgvs,,,,,,,,,,,,,,,,,,,"',
	'-arg', '"-splicing 5,-splicing 5,,,,,,,,,,,,,,,,,,,"',
	'-nastring', '"NA"', '-polish', '-otherinfo']
	
	removeSAM = ['rm', outputSAM]
	
	removeoutputBAM = ['rm', outputBAM]
	
	removeoutputSorstedBAM = ['rm', outputSorstedBAM]
	
	removeSorstedUniqueBAM = ['rm', SorstedUniqueBAM]
	
	################ Now writing the job scripts ######################
	
	ScriptDirectory = dirpath + "/" + "Jobs" + "/" + myProject + "/"
		
	if not os.path.exists(ScriptDirectory):
		os.mkdir(ScriptDirectory)
		print("Directory " , ScriptDirectory ,  " Created ")
	else:    
		print("Directory " , ScriptDirectory ,  " already exists")
	
	jobScript = open( ScriptDirectory + sample + ".sh","wra+") # add path here
	
	
	jobScript.write(str(interpreter))
	jobScript.write(" \n")
	jobScript.write(str(nodesandcores))
	jobScript.write(" \n")
	jobScript.write(str(walltime))
	jobScript.write(" \n")
	jobScript.write(str(pbs))
	jobScript.write(" \n")
	jobScript.write(str(abe))
	jobScript.write(" \n")
	jobScript.write(str(email))
	jobScript.write(" \n")
	jobScript.write(" ")
	jobScript.write(" \n")
	
	ScriptWriter(BWAcommand)
	ScriptWriter(samtoolsViewCommand)
	ScriptWriter(samtoolsSortCommand)
	ScriptWriter(samtoolsINDEXCommand)
	ScriptWriter(picardMARKDUPLICATESCommand)
	ScriptWriter(samtoolsINDEXCommand2)
	ScriptWriter(gatkBaserecalibratorCommand)
	ScriptWriter(gatkAnalyseCovariatesCommand)
	ScriptWriter(gatkApplyBQSRCommand)
	ScriptWriter(gatkHaplotypeCaller_VCF_Command)
	ScriptWriter(gatkVariantFiltration_VCF_Command)
	ScriptWriter(gatkHaplotypeCaller_gVCF_Command)
	ScriptWriter(removeSAM)
	ScriptWriter(removeoutputBAM)
	ScriptWriter(removeoutputSorstedBAM)
	ScriptWriter(removeSorstedUniqueBAM)
	ScriptWriter(convert2annovarCommand)
	ScriptWriter(annovarCommand)
	
	jobScript.write(" \n")
	jobScript.write("exit")
	jobScript.write(" \n")
		
	jobScript.close()	

	################ Submit Jobs to Q #####################################
	
	job = ScriptDirectory + sample + ".sh"
	
	submission = ['qsub', job]
	
	subprocess.call(submission)
	
	time.sleep(3)

