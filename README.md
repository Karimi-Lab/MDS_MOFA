# MDS_MOFA
Scripts and data used in the Multi modal MDS analysis paper

## Summary of the paper


<img width="354" alt="image" src="https://github.com/Karimi-Lab/MDS_MOFA/assets/98902126/8c546da9-765f-4e69-85e7-2d1285e3d2b9">



## How to use the scripts
  Create list of matrices of different biological views to run MOFA, and fit MOFA model
  Downstream analysis of MOFA to investigate heterogeneity and relationships of different modalities with latent factors
  Script to analyse association between mutation (SF3B1 and SRSF2) and gene set signatures and their clinical outcome
  Boxplots of gene signatures that significantly associated with mutations based on mutation status

## Example figure generation 

### Generate heatmap of important features with high weights for each views per factor

<img width="262" alt="image" src="https://github.com/Karimi-Lab/MDS_MOFA/assets/98902126/62b5b662-5895-4b83-a0a2-a08f372daae4">

Run the "Important weights per factors" chunk in _BMMNC_MOFA.Rmd_. Significant features were determined based on absolute weight over 0.5. 

### Generate heatpmap of important features for each factor 

<img width="249" alt="image" src="https://github.com/Karimi-Lab/MDS_MOFA/assets/98902126/7aeed141-277c-4ffd-acd1-d59234d3af16">

Run the "Important features for factor 1" chunk in _BMMNC_MOFA.Rmd_. Significant features were determined based on scaled absolute weight over 0.5 for each factor. Weights sorted based on factor values.
