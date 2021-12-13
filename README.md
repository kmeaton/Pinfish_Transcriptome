# Pinfish transcriptome and Tag-Seq read processing
The scripts in this repository were used to assemble and annotate the transcriptome of the pinfish (_Lagodon rhomboides_) and to analyze Tag-Seq expression data from Eaton et al. (2021). 

### Pinfish transcritpome assembly
The transcriptome was assembled from raw Illumina reads, which were processed by running the following scripts (in this order):

`run_script trimmomatic.sh`

`run_script trinity.sh`

`run_script`
