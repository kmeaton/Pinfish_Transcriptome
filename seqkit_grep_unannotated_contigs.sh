#!/bin/sh

#load the module
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load seqkit/0.10.1

cd /home/aubkme/Pinfish_Transcriptome/new_feb_2021_trimmomatic/annotation/

#seqkit grep allows you to search an input file and output records with names matching those found in an input file (-f)
#the -v option inverse matches
seqkit grep -f ./intermediate_files_do_not_delete/trinity_id_ensembl_id_gene_name.csv cdhit_trinityids_only_copy.fa -o transcripts_annotated_sparus_aurata.fa
seqkit grep -v -f ./intermediate_files_do_not_delete/trinity_id_ensembl_id_gene_name.csv cdhit_trinityids_only_copy.fa -o unannotated_transcripts.fa
