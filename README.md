# Pinfish transcriptome and Tag-Seq read processing
The scripts in this repository were used to assemble and annotate the transcriptome of the pinfish (_Lagodon rhomboides_) and to analyze Tag-Seq expression data from Eaton et al. (2021). 

### Pinfish transcriptome assembly
The transcriptome was assembled from raw Illumina reads, which were processed by running the following scripts on the Alabama Supercomputer as such (in this order):

`run_script trimmomatic.sh`

`run_script trinity.sh`

`run_script cd-hit.sh`

This generates our assembled transcriptome. Despite the fact that Transdecoder is often used in transcriptome assembly pipelines, we did not choose to include it in our analysis, because Transdecoder specifically identifies ORFs and coding sequence and removes 5' and 3' UTRs and other noncoding sequence. For downstream analyses, it was necessary to have 3' UTRs included, as Tag-Seq preferentially sequences the extreme 3' end of transcripts. 

### Pinfish transcriptome annotation

Prior to annotation, we removed extraneous info from sequence headers using the following script:

`run_script rename.sh`
