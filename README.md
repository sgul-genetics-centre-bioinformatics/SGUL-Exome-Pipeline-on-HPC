# SGUL-Exome-Pipeline-on-HPC

Scripts for QC, Alignemnt, Variant Calling and Joint Genotyping of Exome Sequence Data. Following the GATK Best Practices guidelines.

Sample fastq data (paired end R1 and R2 files) is stored in individual sample folders in Unaligned as input for pipeline

e.g.

- ./Unaligned

   - ./MyProject1
      - ./MYsample1
      - ./MYsample2
      - ./MYsample3
      


   - Directory/file structure:
      -  /Exomes
          
           - ./Aligned
           - ./Unaligned
           - ./VCF
           - ./Annotated
           - ./gVCF
           - ./Jobs


   - Needs the following resources to run: For more information check dependencies.py
       - FastQC
       - MultiQC
       - python2.7.5 These scripts have been developed and tested using python2.7. Feel free to update them.
       - bwa
       - human_g1k_v37.fasta
       - samtools
       - java
       - picard-2.815
       - gatk-4.0.4.0
       - Genome_reference_files/common_all.vcf
       - Mills_and_1000G_gold_standard.indels.hg19_modified.vcf
       - 1000G_phase1.snps.high_confidence.hg19.sites.vcf
       - BroadExACExomeIntervlas.bed (Exome Target, the this case the Broad definition)
       
       
       ## Usage
       
       ### in the Exome directory:
       
       python pipeline_Exome_SGUL_HPC_v1.1.py <prject directory>

       
       
