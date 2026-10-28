# The contribution of transposable element insertions to genetic diversity in Aedes aegypti populations

This repository hosts the scripts used to run the pipeline to generate the data and make the figures for the paper “The contribution of transposable element insertions in Aedes aegypti populations”.

The Reference Genome used for these analyses → VectorBase-54_AaegyptiLVP_AGWG_Genome

The consensus TE file added to this rep was from the paper Daron et al. 2024. 

The Repetitive regions annotations came from the paper Love et al. 2023.

## General description of directories within this repo:
- IR_genes: contains the list of genes found to be important for insecticide resistance in Aedes aegypti. These lists were made in Love et at. 2023.
- annotation_genome: vectorbase annotation of the Aedes aegypti genome used for analysis of gene and exon overlaps.
- consensus_TEs: TE library obtained from the paper Daron et al. 2024.
- final_VCF_TEs: final VCF files of reference and nonreference TEs made in this paper after all filtering and post-processing steps.
- nonreference: contains all files used to generate dataset of nonreference TEs. This includes the scripts used to run McClintock, post-processing steps to overlap RetroSeq and TEforest and make the final VCF for non-reference TEs.
- onecodetofindthemall: results obtained from running the onecodetofindthemall pipeline using RepeatMasker results from reference genome and consensus sequences.
- reference: contains all files used to generate dataset of reference TEs. Scripts used to run McClintock, make the location file from onecodetofindthemall results to run mcclintock, make final VCF, so that this one could be use in the genotyping step, and finally all the scripts used to genotype these reference TEs. 
- reference_genome: reference genome file of Aedes aegypti used in this project.
- running_analysis: all scripts made to run final analysis using VCFs for reference and non-reference TEs. These were used to generate all the figures present in the paper. 
