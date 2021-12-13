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

Then, we iteratively searched our transcript sequences via BLASTn against databases of known fish transcripts from (1) _Sparus aurata_, (2) _Larimichthys crocea_, (3) _Oreochromis niloticus_, and (4) _Danio rerio_, filtering the results after each search to only include high-quality matches, using the following scripts:

`run_script sparus_blastn.sh`

`run_script filter_blastn_results_sparus.sh`

`run_script croaker_blastn.sh`

`run_script filter_blastn_results_croaker.sh`

`run_script tilapia_blastn.sh`

`run_script filter_blastn_results_tilapia.sh`

`run_script zebrafish_blastn.sh`

`run_script filter_blastn_results_zebrafish.sh`

Finally, to annotate any remaining unannotated transcripts, we searched our transcript sequences via BLASTx against the UniProt-SwissProt and TrEMBL databases (release 2021_01), running the following scripts in this order:

`run_script uniprot_sprot_blastx.sh`

`run_script filter_blastx_results_uniprot_sprot.sh`

`run_script uniprot_trembl_blastx.sh`

`run_script filter_blastx_results_uniprot_trembl.sh`
