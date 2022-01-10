#!/bin/sh

# KM Eaton, Auburn University, 2021
# Code associated with Eaton et al. 2022 Frontiers in Ecology and Evolution
# This code was run on the Alabama Supercomputer. 

# You HAVE to run this program from the scratch folder because the files it generates are huge! 
# Do this by cd to /scratch/eaton_trinity, and run the script from a scripts folder in there
# The following parameters worked well:
# Queue: large    Wall time: 100:00:00   Mem: 150GB   Cores: 32

# Load the module
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load trinity/2.9.1

# Run the program
# Options:
# seqType specifies the type of reads (fa or fq)
# left specifies the file containing mate 1 reads
# right specifies the file containing mate 2 reads
# CPU specifies the number of CPUs to use
# max_memory specifies the max memory to use by Trinity where limiting can be enabled, specify in GB, such as 120G
# output specifies the name of the directory for output files - the program will create it if it doesn't already exist. The documentation says to include "trinity" in the name of the directory as a safety precaution. I don't know what that means, but I will do it.
# min_contig_length specifies the minumum assembled contig length to report (will not report anything smaller)
# full_cleanup makes the program only retain the Trinity fasta file, and renames it as ${output_dir}.Trinity.fasta
Trinity --seqType fq --left mate1_trimmomatic_paired.fq --right mate2_trimmomatic_paired.fq --CPU 32 --max_memory 150G --output /scratch/aubkme/trinity_trimmomatic_reads/output_trinity_trimmomatic --min_contig_length 300 --verbose --full_cleanup
