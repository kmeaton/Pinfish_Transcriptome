#!/bin/bash

#load the module
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load transdecoder/5.5.0

#run transdecoder
#for transdecoder.longorfs, the default is to identify only orfs that are at minimum 100 amino acids long
#can be lowered with the -m parameter, but the rate of false positive orf predictions increases drastically wiht shorter minimum length criteria
TransDecoder.LongOrfs -t ~/Pinfish_Transcriptome/new_feb_2021_trimmomatic/trinity/output_trinity_trimmomatic.Trinity.fasta
TransDecoder.Predict -t ~/Pinfish_Transcriptome/new_feb_2021_trimmomatic/trinity/output_trinity_trimmomatic.Trinity.fasta
