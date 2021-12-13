#!/bin/sh

# KM Eaton, Auburn University, 2021
# Code associated with Eaton et al. 2021 Frontiers in Ecology and Evolution
# This code was run on the Alabama Supercomputer. 

# Source and load the required modules for this script
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load blast+/2.6.0

# Create the blast database
makeblastdb -in uniprot_sprot.fasta -input_type fasta -dbtype prot -out uniprot_sprot_16june21.fasta.DB

blastx -db uniprot_sprot_16june21.fasta.DB -query zebrafish_unannotated_contigs.fa -out uniprot_sprot.OUT -evalue 1e-10 -num_threads 16 -outfmt '6 std qcovhsp' -max_target_seqs 5
