#!/bin/bash
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load blast+/2.6.0

input_database_file="/home/aubkme/Pinfish_Transcriptome/new_feb_2021_trimmomatic/annotation/gilthead_seabream_biomart_peptides_25feb21.fasta"
database="sparus_aurata_ensembl_prot"
transcriptome="/home/aubkme/Pinfish_Transcriptome/new_feb_2021_trimmomatic/annotation/cdhit_corrected_10feb21.fa"
output="blastx_pin_v_s_aurata_27feb21.blast"
annotation_folder="/home/aubkme/Pinfish_Transcriptome/new_feb_2021_trimmomatic/annotation"

cd $annotation_folder

#only run the following line if it's your first time blasting to this database
#makeblastdb -in $input_database_file -input_type fasta -dbtype prot -out $database

blastx -db $database -query $transcriptome -out $output -evalue 1e-5 -num_threads 16 -outfmt '6 std qcovhsp' -max_target_seqs 5  
