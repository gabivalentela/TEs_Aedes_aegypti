## Scripts Used to Generate the Final VCF of Non-Reference TEs

The `running_mcclintock` directory contains the commands used to run RetroSeq on all samples, as well as the command used to run TEforest.

The `post-processing_nonreference_TEs` directory contains all scripts used to generate the final VCF of non-reference TEs, including:

- Generating the initial VCF from RetroSeq results.
- Generating the initial VCF from TEforest results.
- Combining and overlapping the RetroSeq and TEforest results.
- Removing fixed TEs and singletons from the final dataset.
- Performing variation analyses to determine which variation value to use.
