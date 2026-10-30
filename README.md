# The Contribution of Transposable Element Insertions to Genetic Diversity in *Aedes aegypti* Populations

This repository contains the scripts and files used to generate the datasets, perform the analyses, and produce the figures presented in the paper, "The Contribution of Transposable Element Insertions to Genetic Diversity in *Aedes aegypti* Populations."

**Reference genome:** VectorBase-54_AaegyptiLVP_AGWG_Genome

**TE consensus library:** Obtained from Daron et al. (2024).

**Repetitive region annotations:** Obtained from Love et al. (2023).

## General Description of Directories

The repository is organized into the following directories:

- **`IR_genes`** – Contains the lists of genes known to be associated with insecticide resistance in *Aedes aegypti*. These lists were obtained from Love et al. (2023).

- **`annotation_genome`** – Contains the VectorBase annotations for the *Aedes aegypti* reference genome, used to analyze TE overlaps with genes and exons.

- **`consensus_TEs`** – Contains the TE consensus library obtained from Daron et al. (2024).

- **`final_VCF_TEs`** – Contains the final VCF files for reference and non-reference TEs generated in this study after all filtering and post-processing steps.

- **`nonreference`** – Contains all files and scripts used to generate the dataset of non-reference TEs. This includes the commands used to run McClintock and the post-processing scripts used to combine and overlap RetroSeq and TEforest results and generate the final VCF.

- **`onecodetofindthemall`** – Contains the results obtained from running the One Code to Find Them All pipeline using RepeatMasker results from the reference genome and TE consensus sequences.

- **`reference`** – Contains all files and scripts used to generate the dataset of reference TEs. This includes the commands used to run McClintock, the scripts used to generate the location files from One Code to Find Them All results required by McClintock, the scripts used to generate the initial VCF, and the scripts used to genotype reference TEs.

- **`reference_genome`** – Contains the *Aedes aegypti* reference genome used in this project.

- **`running_analysis`** – Contains all scripts used to perform the final analyses using the reference and non-reference TE VCFs. These scripts were used to generate the figures and results presented in the paper.
