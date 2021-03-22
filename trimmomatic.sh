#!/bin/sh

#load the modules
source/opt/asn/etc/asn-bash-profiles-special/modules.sh
module load fastqc/0.10.1
module load trimmomatic

java -jar /mnt/beegfs/home/aubmxa/.conda/envs/BioInfo_Tools/share/trimmomatic-0.39-1/trimmomatic.jar PE -threads 10 -phred33 \
	~/Pinfish_Transcriptome/raw_reads/all_mate_1.fq ~/Pinfish_Transcriptome/raw_reads/all_mate_2.fq \
	~/Pinfish_Transcriptome/trimmomatic_reads/mate1_paired.fq ~/Pinfish_Transcriptome/trimmomatic_reads/mate1_unpaired.fq \
	~/Pinfish_Transcriptome/trimmomatic_reads/mate2_paired.fq ~/Pinfish_Transcriptome/trimmomatic_reads/mate2_unpaired.fq \
	ILLUMINACLIP:Trimmomatic_Primers.fa:2:35:10 LEADING:30 TRAILING:30 MINLEN:25

fastqc ~/Pinfish_Transcriptome/trimmomatic_reads/mate*.fq
