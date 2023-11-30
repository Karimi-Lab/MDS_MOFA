# Script to merge MDS data
library(readxl)
setwd("~/OneDrive/chemo_pred/Scripts/")
# Load gene sets 
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

clinical_df <- as.data.frame(readxl::read_xlsx("../Data/MDS_Ogawa_Clinical_info.xlsx", sheet = "Survival"))
clinical_df$ID <- paste0("PV", clinical_df$ID)
rownames(clinical_df) <- clinical_df$ID
clinical_df$ID <- NULL

singscores_list <- read_excel_allsheets("../Results/MDS_sign_scores_230530.xlsx")
singscore_df <- do.call(args = singscores_list,
                        what = rbind)
sing_names <- c()
for(x in names(singscores_list)){
  if (dim(singscores_list[[x]])[1] > 0){
    sing_names <- c(sing_names, paste0(x, "_", singscores_list[[x]][,1]))
  }
}
rownames(singscore_df) <- sing_names

# Preprocess genotype df
genotype_df <- read.table("../Data/MDS_Ogawa_genotype.txt", header = F, sep = "\t")
colnames(genotype_df) <- c("Gene", "Type", "ID")
genotype_df$Type <- gsub(" ", ".", genotype_df$Type)
genotype_df$Mutation <- paste(genotype_df$Gene, genotype_df$Type, sep = "_")

# Convert genotype df to wide format
# Keep number of mutations
genotype_df <- genotype_df %>%
  dplyr::group_by(ID, Mutation) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") 
genotype_wide <- as.data.frame(tidyr::pivot_wider(genotype_df, id_cols = "ID", names_from = "Mutation", values_from = n))
rownames(genotype_wide) <- genotype_wide$ID
genotype_wide$ID <- NULL
# Replace NAs with 0
genotype_wide[is.na(genotype_wide)] <- 0
# Order columns aplhabetically
# Good for keeping SVs, WHO subtypes etc. grouped together
genotype_wide <- genotype_wide %>% select(order(colnames(genotype_wide)))


# BMMNC #####
load("../Data/MDS_BMMNC_antigens.Rda")

load("../Data/TE_class_rpm.BMMNC.Rda")
colnames(TE_class_rpm)[1] <- "TE_name"
load("../Data/TE_family_rpm.BMMNC.Rda")
colnames(TE_family_rpm)[1] <- "TE_name"
load("../Data/TE_subfamily_rpm.BMMNC.Rda")
colnames(TE_subfamily_rpm)[1] <- "TE_name"
bmmnc_te_df <- rbind(TE_class_rpm, TE_family_rpm, TE_subfamily_rpm)
# Remove question marked TEs
bmmnc_te_df <- bmmnc_te_df[bmmnc_te_df$TE_name[!grepl("\\?", bmmnc_te_df$TE_name)],]
rownames(bmmnc_te_df) <- bmmnc_te_df$TE_name
bmmnc_te_df$TE_name <- NULL
# remove subfamilies
bmmnc_te_df <- bmmnc_te_df[!grepl(pattern = ".*:.*:.*", rownames(bmmnc_te_df)),]

singscore_bmmnc <- singscore_df[,colnames(singscore_df)[grepl("BMMNC", colnames(singscore_df))]]
colnames(singscore_bmmnc) <- unlist(lapply(colnames(singscore_bmmnc), FUN = function(x){strsplit(x, split = "_")[[1]][2]}))
# Filter for samples in TE matrices
singscore_bmmnc <- singscore_bmmnc[,colnames(bmmnc_te_df)]

clinical_bmmnc <- clinical_df[colnames(bmmnc_te_df),]
clinical_bmmnc <- clinical_bmmnc[complete.cases(clinical_bmmnc),]

bmmnc_te_df <- bmmnc_te_df[,rownames(clinical_bmmnc)]
singscore_bmmnc <- singscore_bmmnc[,rownames(clinical_bmmnc)]
antigens_bmmnc <- as.data.frame(fpkm_MDS_tumor_filt[,rownames(clinical_bmmnc)])

# Merge all data modalities
bmmnc_genotype <- genotype_wide[which(rownames(genotype_wide) %in% colnames(singscore_bmmnc)),]
merged_bmmnc <- plyr::rbind.fill(as.data.frame(t(clinical_bmmnc)), as.data.frame(t(bmmnc_genotype)), singscore_bmmnc, bmmnc_te_df, antigens_bmmnc)
rownames(merged_bmmnc) <- c(colnames(clinical_bmmnc), colnames(bmmnc_genotype), rownames(singscore_bmmnc), rownames(bmmnc_te_df), rownames(antigens_bmmnc))

# CD34 #####
load("../Data/MDS_CD34_antigens.Rda")

load("../Data/TE_class_rpm.CD34.Rda")
colnames(TE_class_rpm)[1] <- "TE_name"
load("../Data/TE_family_rpm.CD34.Rda")
colnames(TE_family_rpm)[1] <- "TE_name"
load("../Data/TE_subfamily_rpm.CD34.Rda")
colnames(TE_subfamily_rpm)[1] <- "TE_name"
cd34_te_df <- rbind(TE_class_rpm, TE_family_rpm, TE_subfamily_rpm)
# Remove question marked TEs
cd34_te_df <- cd34_te_df[cd34_te_df$TE_name[!grepl("\\?", cd34_te_df$TE_name)],]
rownames(cd34_te_df) <- cd34_te_df$TE_name
cd34_te_df$TE_name <- NULL
# remove subfamilies
cd34_te_df <- cd34_te_df[!grepl(pattern = ".*:.*:.*", rownames(cd34_te_df)),]

singscore_cd34 <- singscore_df[,colnames(singscore_df)[grepl("CD34", colnames(singscore_df))]]
colnames(singscore_cd34) <- unlist(lapply(colnames(singscore_cd34), FUN = function(x){strsplit(x, split = "_")[[1]][2]}))
# Filter for samples in TE matrices
singscore_cd34 <- singscore_cd34[,colnames(cd34_te_df)]

clinical_cd34 <- clinical_df[colnames(cd34_te_df),]
clinical_cd34 <- clinical_cd34[complete.cases(clinical_cd34),]

cd34_te_df <- cd34_te_df[,rownames(clinical_cd34)]
singscore_cd34 <- singscore_cd34[,rownames(clinical_cd34)]
antigens_cd34 <- as.data.frame(fpkm_MDS_tumor_filt[,rownames(clinical_cd34)])

# Merge all data modalities
cd34_genotype <- genotype_wide[which(rownames(genotype_wide) %in% colnames(singscore_cd34)),]
merged_cd34 <- plyr::rbind.fill(as.data.frame(t(clinical_cd34)), as.data.frame(t(cd34_genotype)), singscore_cd34, cd34_te_df, antigens_cd34)
rownames(merged_cd34) <- c(colnames(clinical_cd34), colnames(cd34_genotype), rownames(singscore_cd34), rownames(cd34_te_df), rownames(antigens_cd34))