# MDS_MOFA
Scripts and data used in the Multi modal MDS analysis paper

<img width="394" alt="image" src="https://github.com/Karimi-Lab/MDS_MOFA/assets/98902126/24ef4743-5691-4e2f-b248-a61b989e0245">

## How to use the scripts
[matrices_preparation_for_MOFA](Scripts/matrices_preparation_for_MOFA.R): Create list of matrices of same samples' different modalities for both cohort separately 
[BMMNC_MOFA_04012023.Rmd](Scripts/BMMNC_MOFA_04012023.Rmd): Downstream analysis of latent factors from MOFA, including survival analysis for BMMNC cohort 
[CD34RNAseq_MOFA_04012023.Rmd](Scripts/CD34RNAseq_MOFA_04012023.Rmd): Downstream analysis of latent factors from MOFA, including survival analysis for CD34+ RNASeq cohort 
[Characterization_of_mutations_relation_with_geneSets.R](Scripts/Characterization_of_mutations_relation_with_geneSets.R): Characterization of SF3B1 & SRSF2 mutation association with gene sets scores, including survival analysis 

## Example figure generation
### Heatmap of important features within views per factors (_Fig 2 a-b_)

<img width="717" alt="image" src="https://github.com/Karimi-Lab/MDS_MOFA/assets/98902126/bf464642-0538-498b-bcaa-415f975c48b7">

Run the "Important features per factors" chunk in _BMMNC_MOFA_04012023.Rmd_. Important features were selected by having absolute weights over 0.5. 

## Heatmap of important features within Factor 1 (_Fig 2 c-d_)

<img width="1214" alt="image" src="https://github.com/Karimi-Lab/MDS_MOFA/assets/98902126/5997d49b-5157-4e5e-a25e-856e558d7752">

Run the "Important features for Factor 1" chunk in _BMMNC_MOFA_04012023.Rmd_. Important features were selected by having absolute scaled weights over 0.5 for Factor 1.

## Correlogram of mutations and gene sets (_Fig 5a & 6a_)

<img width="401" alt="image" src="https://github.com/Karimi-Lab/MDS_MOFA/assets/98902126/98537c86-84b3-4b3b-85e5-01dc892e9250">

The "_corrTabs_" function takes gene sets scores and mutations and calcutes their correlation coeffecients in addition to p-values from Wilcox rank-sum test (_Characterization_of_mutations_relation_with_geneSets.R_). Heatmap then generated for each mutation separately. Stars indicating p-values manually added.

## Boxplots of significant gene sets and mutation status (_Fig 5b-d & 6b-d_)

<img width="407" alt="image" src="https://github.com/Karimi-Lab/MDS_MOFA/assets/98902126/c2e00132-e1e7-41b9-8130-bb0275162209">

The "_my_Plot_class_" function use significant gene sets names as an argument and create boxplots for each of them with mutation status. For each mutation function run separately (_Characterization_of_mutations_relation_with_geneSets.R_). 

## Kaplan-Meier plots for gene set based on mutation status

<img width="512" alt="image" src="https://github.com/Karimi-Lab/MDS_MOFA/assets/98902126/2727563d-4cfc-46a5-b40d-b809c6d4a964">

The "BMMNC/CD34+ Clinical Outcomes for Aging signatures of patients based on SF3B1 status" parts in the _Characterization_of_mutations_relation_with_geneSets.R_ create Kaplan-Meier plots for specific Inflammation/Aging gene sets (Inflammator chemokines and Inflammatory cytokines, first splitting them based on SF3B1 status (WT and mutated) then splitting mutated ones into "high" and "low". Pairwise p-values and legend titles manually added. 

