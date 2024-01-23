# MDS_MOFA
Scripts and data used in the Multi modal MDS analysis paper

<img width="394" alt="image" src="https://github.com/Karimi-Lab/MDS_MOFA/assets/98902126/24ef4743-5691-4e2f-b248-a61b989e0245">

## How to use the scripts
[matrices_preparation_for_MOFA](Scripts/matrices_preparation_for_MOFA.R) : Create list of matrices of same samples' different modalities for both cohort separately
[BMMNC_MOFA_04012023.Rmd](Scripts/BMMNC_MOFA_04012023.Rmd) : Downstream analysis of latent factors from MOFA, including survival analysis for BMMNC cohort
[CD34RNAseq_MOFA_04012023.Rmd](Scripts/CD34RNAseq_MOFA_04012023.Rmd) : Downstream analysis of latent factors from MOFA, including survival analysis for CD34+ RNASeq cohort
[Characterization_of_mutations_relation_with_geneSets.R](Scripts/Characterization_of_mutations_relation_with_geneSets.R) : Characterization of SF3B1 & SRSF2 mutation association with gene sets scores, including survival analysis
