##############################################################################################################
### Project: Multi-modal data analysis reveals new risk and protective factors for Myelodysplastic Syndromes 
### Script purpose: Preparation of views for MOFA model  
### Author: Sila Gerlevik
### Date: 13 June, 2023
##############################################################################################################

# Tutorial resource: https://biofam.github.io/MOFA2/tutorials.html

# Packages -------------
library(MOFA2)
library(ggplot2)
library(data.table)
library(tidyverse)
library(stringr)
library(caret)
library(dplyr)
library(survival)
library(survminer)

## -- ##
# nested list of matrices, where the first index refers to the view and the second index refers to the group.
# samples are stored in the rows and features are stored in the columns.
# Missing values must be filled with NAs, including samples missing an entire view


# Data preparation ----------------------
load("data/MDS_BMMNC_data/MDS_BMMNC_clinical_df_newCat.Rda")
load("data/MDS_BMMNC_data/MDS_BMMNC_genotype_df_newCat.Rda")
load("data/MDS_BMMNC_data/MDS_BMMNC_singscores_df_newCat.Rda")
load("data/MDS_BMMNC_data/MDS_BMMNC_TE_df_newCat.Rda")


# clinical features divided into numerical and categorical data for defining distributions of data for MOFA ----
clinical_category_bmmnc <- clinical_bmmnc %>% select("Sex", "AML")
clinical_num_bmmnc <- clinical_bmmnc %>%  select("Age", "Blast", "Hb", "ANC", "Plt")

clinical_category_bmmnc[] <- lapply(clinical_category_bmmnc, as.numeric)
clinical_num_bmmnc[] <- lapply(clinical_num_bmmnc, as.numeric)
clinical_num_bmmnc$Plt <- as.numeric(clinical_num_bmmnc$Plt)
clinical_num_bmmnc$ANC <- as.numeric(clinical_num_bmmnc$ANC)
clinical_num_bmmnc$Hb <- as.numeric(clinical_num_bmmnc$Hb)

# gene-signature data ------------------------------
# check NA columns
ind <- rowSums(is.na(singscore_bmmnc)) == ncol(singscore_bmmnc)
singscore_bmmnc <- singscore_bmmnc[!ind,]

# get gene sets 
geneSignature <- singscore_bmmnc
gene_sets <- gsub("_.*$", "", rownames(geneSignature))

# split each gene sets
list_of_genesig <- split(geneSignature, f=gene_sets)
list_of_genesig <- list_of_genesig[c("Aging", "Cellular", "Immunology")]

# make them list of matrices for MOFA object
list_of_genesig$Aging <- as.matrix(list_of_genesig$Aging)
list_of_genesig$Cellular <- as.matrix(list_of_genesig$Cellular)
list_of_genesig$Immunology <- as.matrix(list_of_genesig$Immunology)


# genotype data -----------------------
# decided to use 0 and 1 as an indication of not having/having mutation in any type on that gene/loci
table(bmmnc_genotype > 1)
bmmnc_genotype[bmmnc_genotype > 1] <- 1

# regulate gene/loci names
who <- colnames(bmmnc_genotype)[grepl("WHO", colnames(bmmnc_genotype))]
colnames(bmmnc_genotype) <- gsub("WHO_subtype_", "", colnames(bmmnc_genotype))

# getting only gene/loci names
genotypes <- gsub("_.*$", "", colnames(bmmnc_genotype))
genotypes <- unique(genotypes)

# making empty matrix with gene/loci names 
gen_df <- as.data.frame(matrix(0, nrow = 94, ncol = 112))
colnames(gen_df) <- genotypes
rownames(gen_df) <- rownames(bmmnc_genotype)

# combine all mutation types for each gene/loci into 0 or 1 to indicate overall mutation status for that gene/loci
for (i in genotypes) {
  tmp <- as.data.frame(bmmnc_genotype[, which(str_extract(colnames(bmmnc_genotype), "[^_]+") == i)])
  
  if(ncol(tmp) < 2) {
    gen_df[, which(colnames(gen_df) == i)] <- tmp[, 1]
  } else {
    gen_df[, which(colnames(gen_df) == i)] <- ifelse(rowSums(tmp) >= 1, 1, 0)
  }
  
}

# make sample names in same order with other matrices
rownames(gen_df) %in% colnames(geneSignature)
rownames(gen_df) == colnames(geneSignature)
gen_df <- gen_df[colnames(singscore_bmmnc),]

# check zero variables first and remove near zero variables
set.seed(111)
zeros <- nearZeroVar(gen_df)
colnames(gen_df[zeros])

gen_df <- gen_df[,-zeros]

gen_df <- t(gen_df)

# remove unknown TE ---------------------------------

bmmnc_te_df <- bmmnc_te_df[-grep("Unknown", rownames(bmmnc_te_df)),]

# remove the family/class

bmmnc_te_df <- bmmnc_te_df[grepl(":", rownames(bmmnc_te_df)),]

# create list of matrices with each views for MOFA model
matrices <- list(TE = as.matrix(bmmnc_te_df),
                 mutation = as.matrix(gen_df),
                 clinical_category = as.matrix(t(clinical_category_bmmnc)),
                 clinical_numeric = as.matrix(t(clinical_num_bmmnc)))

matrices_new <- append(matrices, list_of_genesig)

save(matrices_new, file = "data/list_of_matrices_with_seven_views_for_MOFA.Rda")

# check matrices dimensions
lapply(matrices_new, dim)

# MOFA object and train the model  ----------------------------
set.seed(111)
MOFAobject <- create_mofa(matrices_new)
data_opts <- get_default_data_options(MOFAobject) # getting default options to make some changes or applied as it is
data_opts$scale_views <- T # scale views to have the same unit variance
data_opts

model_opts <- get_default_model_options(MOFAobject)
# setting likelihoods for mutation and categorical clinical data
model_opts$likelihoods[3] <- "bernoulli" 
model_opts$likelihoods[2] <- "bernoulli"
model_opts$num_factors <- 15 # initiating model with 15 factors
head(model_opts)

train_opts <- get_default_training_options(MOFAobject)
train_opts$freqELBO <- 1 # to set default ELBO
head(train_opts)

# prepare MOFA object to run with adjusted options
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# run the model ----
MOFAobject.trained <- run_mofa(MOFAobject, "data/1102023_MOFA_without_antigens_scale_view_T.hdf5", use_basilisk=TRUE)

# add metadata to model ------
# The number of rows must match the total number of samples in the model 
sum(MOFAobject.trained@dimensions$N)
bmmnc_df <- clinical_bmmnc # clinical df backup

# it needs a column named "sample" with patients were used in MOFA model
bmmnc_df$sample <- rownames(bmmnc_df)

# add metadata to MOFA 
samples_metadata(MOFAobject.trained) <- bmmnc_df
head(MOFAobject.trained@samples_metadata, n=3)
