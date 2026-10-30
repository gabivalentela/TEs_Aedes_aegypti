## Population Genetic Analyses Using the Final VCFs of Reference and Non-Reference TEs

This directory contains all scripts used to perform the population genetic analyses and generate the figures and results presented in the paper, using the final VCFs of reference and non-reference TEs.

The directory includes the following analyses:

- **`characterization_analysis`** – Generating PCA plots and calculating TE counts across countries.
- **`evolutionary_constraints`** – Identifying TEs with different frequencies across countries and determining which of these overlap immune-related (IR) genes and have evidence of selective sweeps based on SweepFinder2 results.
- **`fst_analysis`** – Calculating FST values for TEs and genome-wide FST values using the full TE dataset, allowing comparisons with SNP-based results.
- **`genome_contribution_lengths_reference_TEs`** – Calculating the proportion of the genome occupied by TEs and assessing differences in genome composition.
- **`high_frequency_TEs`** – Identifying TEs present at high frequencies in the populations.
- **`joint_SFS`** – Generating joint site frequency spectrum (SFS) plots for each pairwise combination of countries.
- **`landscape_analysis`** – Plotting the genomic distribution of TEs across the genome.
- **`ld_analysis`** – Calculating linkage disequilibrium (LD) between TEs and variants known to be associated with insecticide resistance in *Aedes aegypti*.
- **`overlapping_exons`** – Calculating the mean frequency of TEs overlapping and not overlapping exons, as well as examining TE frequencies in repetitive regions.
- **`site_frequency_spectrum`** – Calculating and plotting the site frequency spectrum (SFS) for each country.
