#!/bin/sh

# KM Eaton, Auburn University, 2021
# Code associated with Eaton et al. 2021 Frontiers in Ecology and Evolution
# This code was run on the Alabama Supercomputer. 
# This will take EXTREMELY long if you run it like this, it might be best to break up your query file into several smaller files, parallelize several smaller BLASTs, and then concatenate the results. 

# Source and load the required modules for this script
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load blast+/2.6.0

# Create the blast database
makeblastdb -in uniprot_trembl.fasta -input_type fasta -dbtype prot -out uniprot_trembl_16june21.fasta.DB

blastx -db uniprot_trembl_16june21.fasta.DB -query uniprot_sprot_chordates_unannotated_contigs.fa -out uniprot_trembl.OUT -evalue 1e-10 -num_threads 36 -outfmt '6 std qcovhsp' -max_target_seqs 5
