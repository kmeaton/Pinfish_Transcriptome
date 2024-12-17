# Pinfish transcriptome and Tag-Seq read processing
The scripts in this repository were used to assemble and annotate the transcriptome of the pinfish (_Lagodon rhomboides_) and to analyze Tag-Seq expression data from Eaton et al. (2022). 

## Updated December 2024:
I've now re-annotated the pinfish transcriptome based newly available genomic resources for close relatives in family Sparidae. Scripts used for this updated annotation, and new updated files, are in the directory called updated_2024. 

### Pinfish transcriptome assembly
The transcriptome was assembled from raw Illumina reads, which were processed by running the following scripts on the Alabama Supercomputer as such (in this order):

`run_script trimmomatic.sh`

`run_script trinity.sh`

`run_script cd-hit.sh`

This generates our assembled transcriptome. Despite the fact that Transdecoder is often used in transcriptome assembly pipelines, we did not choose to include it in our analysis, because Transdecoder specifically identifies ORFs and coding sequence and removes 5' and 3' UTRs and other noncoding sequence. For downstream analyses, it was necessary to have 3' UTRs included, as Tag-Seq preferentially sequences the extreme 3' end of transcripts. 

### Pinfish transcriptome annotation

Prior to annotation, we removed extraneous info from sequence headers using the following script:

`run_script rename.sh`

Then, we iteratively searched our transcript sequences via BLASTn against databases of known fish transcripts from (1) _Sparus aurata_, (2) _Larimichthys crocea_, (3) _Oreochromis niloticus_, and (4) _Danio rerio_ (each downloaded from ENSEMBL's BioMart), filtering the results after each search to only include high-quality matches, using the following scripts. To add relevant annotation information, we used a set of files containing ENSEMBL IDs (linked to the known transcript databases), as well as gene names, GO terms, and gene descriptions associated with these sequence IDs. These files (also downloaded from ENSEMBL's BioMart) are located in the annotation_information folder of this repository, and are named \*\_gene_GO.csv.  

`run_script sparus_blastn.sh`

`run_script filter_blastn_results_sparus.sh`

`run_script croaker_blastn.sh`

`run_script filter_blastn_results_croaker.sh`

`run_script tilapia_blastn.sh`

`run_script filter_blastn_results_tilapia.sh`

`run_script zebrafish_blastn.sh`

`run_script filter_blastn_results_zebrafish.sh`

Finally, to annotate any remaining unannotated transcripts, we searched our transcript sequences via BLASTx against the UniProt-SwissProt and TrEMBL databases (release 2021_01), running the following scripts. Again, to add relevant annotation information, we used a set of files containing UniProt IDs linked to gene names and GO terms associated with these sequence IDs. These files were downloaded from UniProt and are also located in the annotation_information folder of this repository. 

`run_script uniprot_sprot_blastx.sh`

`run_script filter_blastx_results_uniprot_sprot.sh`

`run_script uniprot_trembl_blastx.sh`

`run_script filter_blastx_results_uniprot_trembl.sh`

### Final assembly and annotation table

The final assembly is available on NCBI's Transcriptome Shotgun Assembly site, under BioProject Number PRJNA776622.

The final annotation table was reformatted to contain only one line per gene with relevant functional information, and this complete annotation table is uploaded as transcriptome_annotation.csv. This file is the same as Supplementary File 1 from Eaton et al. (2022). The detailed results from the ENSEMBL and UniProt annotations (prior to their concatenation into a single, more concise file) are also uploaded, as ENSEMBL_annotations.csv and UNIPROT_annotations.csv. 

### Tag-Seq read processing

Reads were processed following the pipeline available [here](https://github.com/z0on/tag-based_RNAseq).

Differentially expressed genes were identified for both juveniles and adults with the program DESeq2 in R, using the scripts DESeq2_juveniles.R and DESeq2_adults.R, respectively. 
