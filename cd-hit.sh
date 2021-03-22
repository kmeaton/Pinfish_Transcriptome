#!/bin/sh

#load the module
source /opt/asn/etc/asn-bash-profiles-special/modules.sh

module load anaconda/3-2020.02

#run cd-hit
#cd-hit-est is used for DNA sequences that you want to cluster
#-i is the name of the input file
#-o is the name of the output file
#-c is the desired sequence similarity (0.95 corresponds to 95%)
#-n is the word size, usually 10 or 11 is good for a -c of 0.95
#-d is how many characters of the fasta header you want to keep
#-M is the maximum memory limit, in MB (20000 MB = 20 GB)
#-T is the number of threads, setting it to 0 will use the maximum
cd-hit-est -i ~/Pinfish_Transcriptome/new_feb_2021_trimmomatic/transdecoder/output_trinity_trimmomatic.Trinity.fasta.transdecoder.cds -o ~/Pinfish_Transcriptome/new_feb_2021_trimmomatic/cdhit/cdhit_corrected_10feb21.fa -c 0.95 -n 10 -d 200 -M 20000 -T 0
