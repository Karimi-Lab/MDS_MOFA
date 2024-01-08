## Correlation analysis between mutations and gene sets singscore
## Plot the correlations of all features in terms of absent/present of SF3B1 & SRS2F
## 28 Sep, 2023 
## Sila Gerlevik


# Packages -----
library(ComplexHeatmap)
library(ggsci)
library(readxl)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(dplyr)
library(PerformanceAnalytics)
library(pals)
library(tidyr)
library(grid)
library(gridExtra)
library(plyr)
library(GGally)


# Function for creating correlation matrices along with p-vals
corrTabs <- function(mat1, mat2){
  inf_df <- cbind(mat1, mat2)
  #inf_df <- na.omit(inf_df)
  #inf_df <- mutate_all(inf_df, function(x)as.numeric(x))
  corMat <- list()
  corMat$P <- data.frame()
  corMat$r <- data.frame()
  mat1_feats <- colnames(mat1)
  mat2_feats <- colnames(mat2)
  
  
  inf_df <- t(inf_df)
  
  for (feat1 in mat1_feats){
    for (feat2 in mat2_feats) {
      corMat$P[feat1, feat2] <- tryCatch({
        round(wilcox.test(as.numeric(inf_df[feat1,]) ~ as.factor(inf_df[feat2,]))$p.value, 4)
      },
      error = function(err) {return(NA)}
      )
      
      corMat$r[feat1, feat2] <- tryCatch({
        round(cor.test(as.numeric(inf_df[feat1,]), as.numeric(inf_df[feat2,]))$estimate, 2)
      },
      error = function(err) {return(NA)}
      )
      
    }
  }
  # Convert NA's to insignificant which would show white in the plot
  # If there is one column, complete.cases removes column name
  #corMat$P[corMat$P < 0.0001] <- "<"
  corMat$r[is.na(corMat$r)] <- 0
  
  return(corMat)
}


# Data ------
load("data/merge_df_for_correlogram.Rda")
load("data/mutations_for_correlogram.Rda")

mutations <- as.data.frame(mutations)
merge_df <- as.data.frame(merge_df)


colnames(merge_df)[55] <- "Immunology_Cytolytic score"
merge_df <- merge_df[,-3] # remove inflammaging (chronic inflammation)

# interest of mutations only
mutations$SF3B1 <- as.factor(mutations$SF3B1)
mutations$SRSF2 <- as.factor(mutations$SRSF2)

# get views for heatmap from gene set signatures
merged_features <- as.data.frame(t(merge_df))

merged_features$views <-  str_extract(rownames(merged_features), "[^_]+")


# get correlation results and wilcox significance ------
ttest_disc <- corrTabs(merge_df, mutations)

# add tissue name to differentiation 
colnames(ttest_disc$P) <- c("SF3B1\nBMMNC", "SRSF2\nBMMNC")
colnames(ttest_disc$r) <- c("SF3B1\nBMMNC", "SRSF2\nBMMNC")


# CD34+ cohort -------
# apply same data regulations as before
load("data/merge_validation_for_corr.Rda")
load("data/mutation_validation_for_corr.Rda")

merge_val <- as.data.frame(merge_val)
mutations_val <- as.data.frame(mutations_val)


colnames(merge_val)[55] <- "Immunology_Cytolytic score"
merge_val <- merge_val[,-3]

mutations_val$SF3B1 <- as.factor(mutations_val$SF3B1)
mutations_val$SRSF2 <- as.factor(mutations_val$SRSF2)

merged_features_val <- as.data.frame(t(merge_val))
merged_features_val$views <- str_extract(rownames(merged_features_val), "[^_]+")


ttest_val <- corrTabs(merge_val, mutations_val)

colnames(ttest_val$P) <- c("SF3B1\nCD34", "SRSF2\nCD34")
colnames(ttest_val$r) <- c("SF3B1\nCD34", "SRSF2\nCD34")


# select cellular and aging for both cohort ----
columnsfromboth <- cbind(ttest_disc$r[!grepl("Immunology", rownames(ttest_disc$r)),], 
                         ttest_val$r[!grepl("Immunology", rownames(ttest_val$r)),] )

merged_features <- merged_features[merged_features$views != "Immunology",]


colh = list(views = c("Aging" = "#F39200FF", 
                     "Cell-type" = "#EFD500FF"))
                     
merged_features$views <- gsub("Cellular", "Cell-type", merged_features$views)

# Create the heatmap annotation
haH <- HeatmapAnnotation(
  views = merged_features$views, col = colh, which = "row", show_legend = T, 
  show_annotation_name = F, 
  annotation_legend_param =  list(title = "Views", title_gp = gpar(fontsize = 5, family= "Helvetica"),
                                  labels = c("Aging", "Cell-type"),
                                  labels_gp = gpar(fontsize = 5, family= "Helvetica"), grid_height = unit(2.8, "mm"),
                                  grid_width = unit(2.8, "mm")))

rownames(columnsfromboth) <- gsub(".*\\_", "", rownames(columnsfromboth))
colnames(columnsfromboth)

# columnsfromboth <- columnsfromboth[,c(2,4)]

ht <- Heatmap(as.matrix(columnsfromboth),
              color_space = "RGB",
              # cell_fun = function(j, i, x, y, width, height, fill) {
              #   if(pval_4columns[i, j] < 0.05) {
              #     grid.text("*", x, y)
              #   }
              # },
              row_names_gp = gpar(fontsize = 5, fontfamily = "Helvetica"),
              column_names_gp = gpar(fontsize = 5, fontfamily = "Helvetica", lineheight = 0.78),
              column_split = colnames(columnsfromboth),
              column_title = NULL,
              show_row_dend = F,
              show_column_dend = F,
              right_annotation = haH,
              split = merged_features$views,
              row_title = NULL,
              column_order = c("SF3B1\nBMMNC" ,"SF3B1\nCD34", "SRSF2\nBMMNC", "SRSF2\nCD34"),
              show_column_names = T,
              heatmap_legend_param = list(
                title = "Correlation\nCoefficient", title_gp = gpar(fontsize = 5, family= "Helvetica"),
                labels_gp = gpar(fontsize = 5, family= "Helvetica"), grid_height = unit(2.8, "mm"),
                grid_width = unit(2.8, "mm")
                #legend_height = 2, legend_width = 0.5
              ),
              column_names_rot = 0,
              column_names_centered = T
              # width = ncol(as.matrix(columnsfromboth))*unit(1.8, "cm"),
              # height = nrow(as.matrix(columnsfromboth))*unit(0.3, "cm")
              
) 


pdf("figures/New_figures_2024/allMutations_Heatmaps_Disc_val_agingCellular.pdf", width = 2.84331, height = 2.75591, onefile = T)
draw(ht)
dev.off()

# for immunology ---------------
immunology <- cbind(ttest_disc$r[grepl("Immunology", rownames(ttest_disc$r)),], 
                    ttest_val$r[grepl("Immunology", rownames(ttest_val$r)),] )

rownames(immunology) <- gsub(".*\\_", "", rownames(immunology))

removed_imm <- c("Endothelial cells",	"Exhausted CD8+ T cells", "CD56dim NK cells", "APM",	"APM loss", "Hypoxia",
                 "Immunoproteasome", "JAK-STAT loss", "MAGEs", "TIS",	"Single-gene signatures",	"Negative ICP (manually curated)",	
                 "Positive ICP (manually curated)",	"Tfh","Stemness",	"Exhausted progenitors", "CD8-CD28",	"CD8-CX3CR1",	
                 "CD8-GZMK",	"CD8-LAYN","CD8-ZNF683",	"CD8-SLC4A10",
                 "CART dysfunction",	"IED172 (RNA-seq)")


immunology <- immunology[!rownames(immunology) %in% removed_imm,]

im_ht <- Heatmap(as.matrix(immunology),
                 # cell_fun = function(j, i, x, y, width, height, fill) {
                 #   if(immunology_pval[i, j] < 0.05) {
                 #     grid.text("*", x, y)
                 #   }
                 # },
                 row_order = rownames(immunology),
                 color_space = "RGB",
                 row_names_gp = gpar(fontsize = 5, fontfamily = "Helvetica"),
                 column_names_gp = gpar(fontsize = 5, fontfamily = "Helvetica", lineheight = 0.78),
                 column_title = NULL,
                 show_row_dend = F,
                 show_column_dend = F,
                 column_split = colnames(immunology),
                 row_title = NULL,
                 column_order = c("SF3B1\nBMMNC" ,"SF3B1\nCD34", "SRSF2\nBMMNC", "SRSF2\nCD34"),
                 show_column_names = T,
                 heatmap_legend_param = list(
                   title = "Correlation\nCoefficient", title_gp = gpar(fontsize = 5, family= "Helvetica"),
                   labels_gp = gpar(fontsize = 5, family= "Helvetica"), grid_height = unit(2.8, "mm"),
                   grid_width = unit(2.8, "mm")
                   #legend_height = 2, legend_width = 0.5
                 ),
                 column_names_rot = 0,
                 column_names_centered = T
                 
) 


pdf("figures/New_figures_2024/allMutations_Heatmaps_Disc_Immunology.pdf", width = 1.7346, height = 3.73701, onefile = T)
draw(im_ht)
dev.off()

# get important features ------
signif_SF3B1_cells <- rownames(ttest_disc$P)[which(ttest_disc$P$SF3B1 < 0.05)]
signif_SF3B1_cells_val <- rownames(ttest_val$P)[which(ttest_val$P$SF3B1 < 0.05)]

signif_SRSF2_cells <- rownames(ttest_disc$P)[which(ttest_disc$P$SRSF2 < 0.05)]
signif_SRSF2_cells_val <- rownames(ttest_val$P)[which(ttest_val$P$SRSF2 < 0.05)]


ttest_disc$P < 0.05



# CD34+ Clinical Outcomes for Aging signatures of patients that have SF3B1 and SRSF2 mutation ####
sf3b1mut_val <- mutations_val %>% filter(SF3B1 == 1) %>% rownames()
srsf2mut_val <- mutations_val %>% filter(SRSF2 == 1) %>% rownames()

aging_sf3b1mut_val <- merge_val[sf3b1mut_val, grep("Aging", colnames(merge_val))]
aging_srsf2mut_val <- merge_val[srsf2mut_val, grep("Aging", colnames(merge_val))]

load("data/GSE114922/GSE114922_clinical_df_newCat.Rda")

clinical_df <- as.data.frame(t(clinical_df))
clinical_df$Status <- ifelse(clinical_df$Status == "DEAD", 1, 0)
clinical_df$Status <- as.numeric(clinical_df$Status)
clinical_df$Days <- as.numeric(clinical_df$Days)

sf3b1mut_clinical_val <- clinical_df[sf3b1mut_val,]
#singscore_bmmnc <- as.data.frame(t(singscore_bmmnc))

### SF3B1 mutated ######
mds_sing_val <- aging_sf3b1mut_val
colnames(mds_sing_val) <- gsub(".*\\_", "", colnames(mds_sing_val))


plot_list_v <- list()
i <- 1
lowHighP <- list()

for(sc in colnames(mds_sing_val)){
  tmp <- as.data.frame(mds_sing_val) %>% mutate(quartile = ntile(get(sc), 4))
  tmp$quartile <- plyr::mapvalues(tmp$quartile, from = c(1,2,3,4), to = c("Low", "Med", "Med", "High"))
  
  plot_df <- data.frame(Quartile = tmp$quartile, 
                        Survival = sf3b1mut_clinical_val$Days,
                        Event =  sf3b1mut_clinical_val$Status)
  colnames(plot_df) <- c("Quartile","Survival","Event")
  #plot_df$Event <- ifelse(plot_df$Event == "ALIVE", 0, 1)
  
  plot_df$Quartile <- factor(plot_df$Quartile, levels = c("Low", "Med", "High"), ordered = T)
  
  plot_df <- plot_df %>%
    filter(Quartile != "Med")
  
  km_infChem <- survfit(formula = Surv(time = as.numeric(plot_df$Survival), event = as.numeric(plot_df$Event)) ~ plot_df$Quartile,
                        data = plot_df)
  
  lowHighP[[sc]] <- survminer::surv_pvalue(km_infChem)$pval
  
  g <- ggsurvplot(km_infChem, data = plot_df,
                  conf.int = F, pval = TRUE,
                  fun = function(y) y * 100,
                  legend = "top", xscale = "d_y", break.time.by=365.25*3,
                  xlab = "Time to event", ylab="Survival probability (%)", title= sc, size = 2, censor.size = 15
  )$plot
  
  plot_list_v[[i]] <- g
  i <- i + 1
}

plot_list_v <- list()
i <- 1
for(sc in colnames(mds_sing_val)){
  tmp <- as.data.frame(mds_sing_val) %>% mutate(quartile = ntile(get(sc), 4))
  tmp$quartile <- plyr::mapvalues(tmp$quartile, from = c(1,2,3,4), to = c("Low", "Med", "Med", "High"))
  
  plot_df <- data.frame(Quartile = tmp$quartile, 
                        Survival = sf3b1mut_clinical_val$Days,
                        Event =  sf3b1mut_clinical_val$Status)
  #colnames(plot_df) <- c("Quartile","Survival","Event")
  #plot_df$Event <- ifelse(plot_df$Event == "ALIVE", 0, 1)
  
  plot_df$Quartile <- factor(plot_df$Quartile, levels = c("Low", "Med", "High"), ordered = T)
  
  km_infChem <- survfit(formula = Surv(time = as.numeric(plot_df$Survival), event = plot_df$Event) ~ plot_df$Quartile,
                        data = plot_df)
  
  g <- ggsurvplot(km_infChem, data = plot_df,
                  conf.int = F, pval = paste0("p-value:",round(lowHighP[[sc]],3)),
                  fun = function(y) y * 100,
                  legend = "top", xscale = "d_y", break.time.by=365.25*3,
                  legend.title = " ", size = 0.3, censor.size = 1.5,
                  palette = "npg", pval.size = 1.8, pval.coord = c(0.5, 10),
                  legend.labs = c(paste("High"), paste("Low"), paste("Med")), 
                  surv.plot.height = 0.3,
                  xlab = "Years", ylab="Survival probability (%)")$plot + 
    ggtitle(sc) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          #axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.text.y = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.title.x = element_text(size = 5, family = "Helvetica", color = "black"),
          legend.position = "none", aspect.ratio = 1,
          legend.title = element_text(size = 5, family = "Helvetica", color = "black"),
          legend.text = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.title.y = element_text(size = 5, family = "Helvetica", color = "black"),
          plot.title = element_text(size = 5, family = "Helvetica", color = "black", vjust = -2))
  
  ggsave(paste0("figures/OS_val_SF3B1_mutated_",sc,".pdf"), g, width = 44, height = 44 , units = "mm", dpi = 300)
  
  plot_list_v[[i]] <- g
  i <- i + 1
}

pdf(file = "figures/OS_val_SF3B1mut_Agings_KM.pdf",
    width = 3.46457,
    height = 4.72441, onefile = T, bg = "white")
ggarrange(plotlist = plot_list_v, nrow = 3, ncol = 2, align = "hv", common.legend = T)

dev.off()



### SRSF2 mutated in CD34+ -------

srsf2mut_clinical_val <- clinical_df[srsf2mut_val,]

mds_sing_srsf2_v <- aging_srsf2mut_val
colnames(mds_sing_srsf2_v) <- gsub(".*\\_", "", colnames(mds_sing_srsf2_v))


plot_list_srsf2_v <- list()
i <- 1
lowHighP <- list()

for(sc in colnames(mds_sing_srsf2_v)){
  tmp <- as.data.frame(mds_sing_srsf2_v) %>% mutate(quartile = ntile(get(sc), 4))
  tmp$quartile <- plyr::mapvalues(tmp$quartile, from = c(1,2,3,4), to = c("Low", "Med", "Med", "High"))
  
  plot_df <- data.frame(Quartile = tmp$quartile, 
                        Survival = srsf2mut_clinical_val$Days,
                        Event =  srsf2mut_clinical_val$Status)
  colnames(plot_df) <- c("Quartile","Survival","Event")
  #plot_df$Event <- ifelse(plot_df$Event == "ALIVE", 0, 1)
  
  plot_df$Quartile <- factor(plot_df$Quartile, levels = c("Low", "Med", "High"), ordered = T)
  
  plot_df <- plot_df %>%
    filter(Quartile != "Med")
  
  km_infChem <- survfit(formula = Surv(time = as.numeric(plot_df$Survival), event = as.numeric(plot_df$Event)) ~ plot_df$Quartile,
                        data = plot_df)
  
  lowHighP[[sc]] <- survminer::surv_pvalue(km_infChem)$pval
  
  g <- ggsurvplot(km_infChem, data = plot_df,
                  conf.int = F, pval = TRUE,
                  fun = function(y) y * 100,
                  legend = "top", xscale = "d_y", break.time.by=365.25*3,
                  xlab = "Time to event", ylab="Survival probability (%)", title= sc, size = 2, censor.size = 15
  )$plot
  
  plot_list_srsf2_v[[i]] <- g
  i <- i + 1
}

plot_list_srsf2_v <- list()
i <- 1
for(sc in colnames(mds_sing_srsf2_v)){
  tmp <- as.data.frame(mds_sing_srsf2_v) %>% mutate(quartile = ntile(get(sc), 4))
  tmp$quartile <- plyr::mapvalues(tmp$quartile, from = c(1,2,3,4), to = c("Low", "Med", "Med", "High"))
  
  plot_df <- data.frame(Quartile = tmp$quartile, 
                        Survival = srsf2mut_clinical_val$Days,
                        Event =  srsf2mut_clinical_val$Status)
  #colnames(plot_df) <- c("Quartile","Survival","Event")
  #plot_df$Event <- ifelse(plot_df$Event == "ALIVE", 0, 1)
  
  plot_df$Quartile <- factor(plot_df$Quartile, levels = c("Low", "Med", "High"), ordered = T)
  
  km_infChem <- survfit(formula = Surv(time = as.numeric(plot_df$Survival), event = plot_df$Event) ~ plot_df$Quartile,
                        data = plot_df)
  
  g <- ggsurvplot(km_infChem, data = plot_df,
                  conf.int = F, pval = paste0("p-value:",round(lowHighP[[sc]],3)),
                  fun = function(y) y * 100,
                  legend = "top", xscale = "d_y", break.time.by=365.25*3,
                  legend.title = " ", size = 0.3, censor.size = 1.5,
                  palette = "npg", pval.size = 1.8, pval.coord = c(0.5, 10),
                  legend.labs = c(paste("High"), paste("Low"), paste("Med")), 
                  surv.plot.height = 0.3,
                  xlab = "Years", ylab="Survival probability (%)")$plot + 
    ggtitle(sc) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          #axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.text.y = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.title.x = element_text(size = 5, family = "Helvetica", color = "black"),
          legend.position = "none", aspect.ratio = 1,
          legend.title = element_text(size = 5, family = "Helvetica", color = "black"),
          legend.text = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.title.y = element_text(size = 5, family = "Helvetica", color = "black"),
          plot.title = element_text(size = 5, family = "Helvetica", color = "black", vjust = -2))
  
  ggsave(paste0("figures/OS_val_srsf2_mutated_",sc,".pdf"), g, width = 44, height = 44 , units = "mm", dpi = 300)
  
  plot_list_srsf2_v[[i]] <- g
  i <- i + 1
}

pdf(file = "figures/OS_val_SRSF2mut_Agings_KM.pdf",
    width = 3.46457,
    height = 4.72441, onefile = T, bg = "white")
ggarrange(plotlist = plot_list_srsf2_v, nrow = 3, ncol = 2, align = "hv", common.legend = T)

dev.off()



# BMMNC Clinical Outcomes for Aging signatures of patients that have SF3B1 and SRSF2 mutation ####
sf3b1mut <- mutations %>% filter(SF3B1 == 1) %>% rownames()
srsf2mut <- mutations %>% filter(SRSF2 == 1) %>% rownames()

aging_sf3b1mut <- merge_df[sf3b1mut, grep("Aging", colnames(merge_df))]
aging_srsf2mut <- merge_df[srsf2mut, grep("Aging", colnames(merge_df))]



### SF3B1 mutated in BMMNC-------
load("data/MDS_BMMNC_data/MDS_BMMNC_clinical_df_newCat.Rda")

sf3b1mut_clinical <- clinical_bmmnc[sf3b1mut,]
#singscore_bmmnc <- as.data.frame(t(singscore_bmmnc))

mds_sing <- aging_sf3b1mut
colnames(mds_sing) <- gsub(".*\\_", "", colnames(mds_sing))


plot_list <- list()
i <- 1
lowHighP <- list()

for(sc in colnames(mds_sing)){
  tmp <- as.data.frame(mds_sing) %>% mutate(quartile = ntile(get(sc), 4))
  tmp$quartile <- plyr::mapvalues(tmp$quartile, from = c(1,2,3,4), to = c("Low", "Med", "Med", "High"))
  
  plot_df <- data.frame(Quartile = tmp$quartile, 
                        Survival = sf3b1mut_clinical$OS,
                        Event =  sf3b1mut_clinical$Event)
  colnames(plot_df) <- c("Quartile","Survival","Event")
  #plot_df$Event <- ifelse(plot_df$Event == "ALIVE", 0, 1)
  
  plot_df$Quartile <- factor(plot_df$Quartile, levels = c("Low", "Med", "High"), ordered = T)
  
  plot_df <- plot_df %>%
    filter(Quartile != "Med")
  
  km_infChem <- survfit(formula = Surv(time = as.numeric(plot_df$Survival), event = as.numeric(plot_df$Event)) ~ plot_df$Quartile,
                        data = plot_df)
  
  lowHighP[[sc]] <- survminer::surv_pvalue(km_infChem)$pval
  
  g <- ggsurvplot(km_infChem, data = plot_df,
                  conf.int = F, pval = TRUE,
                  fun = function(y) y * 100,
                  legend = "top", xscale = "d_y", break.time.by=365.25*3,
                  xlab = "Time to event", ylab="Survival probability (%)", title= sc, size = 2, censor.size = 15
  )$plot
  
  plot_list[[i]] <- g
  i <- i + 1
}

plot_list <- list()
i <- 1
for(sc in colnames(mds_sing)){
  tmp <- as.data.frame(mds_sing) %>% mutate(quartile = ntile(get(sc), 4))
  tmp$quartile <- plyr::mapvalues(tmp$quartile, from = c(1,2,3,4), to = c("Low", "Med", "Med", "High"))
  
  plot_df <- data.frame(Quartile = tmp$quartile, 
                        Survival = sf3b1mut_clinical$OS,
                        Event =  sf3b1mut_clinical$Event)
  #colnames(plot_df) <- c("Quartile","Survival","Event")
  #plot_df$Event <- ifelse(plot_df$Event == "ALIVE", 0, 1)
  
  plot_df$Quartile <- factor(plot_df$Quartile, levels = c("Low", "Med", "High"), ordered = T)
  
  km_infChem <- survfit(formula = Surv(time = as.numeric(plot_df$Survival), event = plot_df$Event) ~ plot_df$Quartile,
                        data = plot_df)
  
  g <- ggsurvplot(km_infChem, data = plot_df,
                  conf.int = F, pval = paste0("p-value:",round(lowHighP[[sc]],3)),
                  fun = function(y) y * 100,
                  legend = "top", xscale = "d_y", break.time.by=365.25*3,
                  legend.title = " ", size = 0.3, censor.size = 1.5,
                  palette = "npg", pval.size = 1.8, pval.coord = c(0.5, 10),
                  legend.labs = c(paste("High"), paste("Low"), paste("Med")), 
                  surv.plot.height = 0.3,
                  xlab = "Years", ylab="Survival probability (%)")$plot + 
    ggtitle(sc) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          #axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.text.y = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.title.x = element_text(size = 5, family = "Helvetica", color = "black"),
          legend.position = "none", aspect.ratio = 1,
          legend.title = element_text(size = 5, family = "Helvetica", color = "black"),
          legend.text = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.title.y = element_text(size = 5, family = "Helvetica", color = "black"),
          plot.title = element_text(size = 5, family = "Helvetica", color = "black", vjust = -2))
  
  ggsave(paste0("figures/OS_disc_SF3B1_mutated_",sc,".pdf"), g, width = 44, height = 44 , units = "mm", dpi = 300)
  
  plot_list[[i]] <- g
  i <- i + 1
}

pdf(file = "figures/OS_disc_SF3B1mut_Agings_KM.pdf",
    width = 3.46457,
    height = 4.72441, onefile = T, bg = "white")
ggarrange(plotlist = plot_list, nrow = 3, ncol = 2, align = "hv", common.legend = T)

dev.off()



### SRSF2 mutated in BMMNC-------

srsf2mut_clinical <- clinical_bmmnc[srsf2mut,]
#singscore_bmmnc <- as.data.frame(t(singscore_bmmnc))

mds_sing_srsf2 <- aging_srsf2mut
colnames(mds_sing_srsf2) <- gsub(".*\\_", "", colnames(mds_sing_srsf2))


plot_list_srsf2 <- list()
i <- 1
lowHighP <- list()

for(sc in colnames(mds_sing_srsf2)){
  tmp <- as.data.frame(mds_sing_srsf2) %>% mutate(quartile = ntile(get(sc), 4))
  tmp$quartile <- plyr::mapvalues(tmp$quartile, from = c(1,2,3,4), to = c("Low", "Med", "Med", "High"))
  
  plot_df <- data.frame(Quartile = tmp$quartile, 
                        Survival = srsf2mut_clinical$OS,
                        Event =  srsf2mut_clinical$Event)
  colnames(plot_df) <- c("Quartile","Survival","Event")
  #plot_df$Event <- ifelse(plot_df$Event == "ALIVE", 0, 1)
  
  plot_df$Quartile <- factor(plot_df$Quartile, levels = c("Low", "Med", "High"), ordered = T)
  
  plot_df <- plot_df %>%
    filter(Quartile != "Med")
  
  km_infChem <- survfit(formula = Surv(time = as.numeric(plot_df$Survival), event = as.numeric(plot_df$Event)) ~ plot_df$Quartile,
                        data = plot_df)
  
  lowHighP[[sc]] <- survminer::surv_pvalue(km_infChem)$pval
  
  g <- ggsurvplot(km_infChem, data = plot_df,
                  conf.int = F, pval = TRUE,
                  fun = function(y) y * 100,
                  legend = "top", xscale = "d_y", break.time.by=365.25*3,
                  xlab = "Time to event", ylab="Survival probability (%)", title= sc, size = 2, censor.size = 15
  )$plot
  
  plot_list_srsf2[[i]] <- g
  i <- i + 1
}

plot_list_srsf2 <- list()
i <- 1
for(sc in colnames(mds_sing_srsf2)){
  tmp <- as.data.frame(mds_sing_srsf2) %>% mutate(quartile = ntile(get(sc), 4))
  tmp$quartile <- plyr::mapvalues(tmp$quartile, from = c(1,2,3,4), to = c("Low", "Med", "Med", "High"))
  
  plot_df <- data.frame(Quartile = tmp$quartile, 
                        Survival = srsf2mut_clinical$OS,
                        Event =  srsf2mut_clinical$Event)
  #colnames(plot_df) <- c("Quartile","Survival","Event")
  #plot_df$Event <- ifelse(plot_df$Event == "ALIVE", 0, 1)
  
  plot_df$Quartile <- factor(plot_df$Quartile, levels = c("Low", "Med", "High"), ordered = T)
  
  km_infChem <- survfit(formula = Surv(time = as.numeric(plot_df$Survival), event = plot_df$Event) ~ plot_df$Quartile,
                        data = plot_df)
  
  g <- ggsurvplot(km_infChem, data = plot_df,
                  conf.int = F, pval = paste0("p-value:",round(lowHighP[[sc]],3)),
                  fun = function(y) y * 100,
                  legend = "top", xscale = "d_y", break.time.by=365.25*3,
                  legend.title = " ", size = 0.3, censor.size = 1.5,
                  palette = "npg", pval.size = 1.8, pval.coord = c(0.5, 10),
                  legend.labs = c(paste("High"), paste("Low"), paste("Med")), 
                  surv.plot.height = 0.3,
                  xlab = "Years", ylab="Survival probability (%)")$plot + 
    ggtitle(sc) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          #axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.text.y = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.title.x = element_text(size = 5, family = "Helvetica", color = "black"),
          legend.position = "none", aspect.ratio = 1,
          legend.title = element_text(size = 5, family = "Helvetica", color = "black"),
          legend.text = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.title.y = element_text(size = 5, family = "Helvetica", color = "black"),
          plot.title = element_text(size = 5, family = "Helvetica", color = "black", vjust = -2))
  
  ggsave(paste0("figures/OS_disc_srsf2_mutated_",sc,".pdf"), g, width = 44, height = 44 , units = "mm", dpi = 300)
  
  plot_list_srsf2[[i]] <- g
  i <- i + 1
}

pdf(file = "figures/OS_disc_SRSF2mut_Agings_KM.pdf",
    width = 3.46457,
    height = 4.72441, onefile = T, bg = "white")
ggarrange(plotlist = plot_list_srsf2, nrow = 3, ncol = 2, align = "hv", common.legend = T)

dev.off()





# BMMNC Clinical Outcomes for Aging signatures of patients based on SF3B1 status ####
load("data/merge_df_for_correlogram.Rda")
load("data/mutations_for_correlogram.Rda")

merge_df <- as.data.frame(merge_df)

mutations <- as.data.frame(mutations)
mutations$SF3B1 <- as.factor(mutations$SF3B1)
mutations$SRSF2 <- as.factor(mutations$SRSF2)

sf3b1mut <- mutations %>% filter(SF3B1 == 1) %>% rownames()
sf3b1WT <- mutations %>% filter(SF3B1 == 0) %>% rownames()

# 
colnames(merge_df) <- gsub(".*\\_", "", colnames(merge_df))
agings <- c("Inflammatory chemokines", "Inflammatory cytokines")
# 
# aging_sf3b1mut <- merge_df[sf3b1mut, agings]
# aging_sf3b1WT <- merge_df[sf3b1WT, agings]
# aging_sf3b1WT$quartile <- rep("SF3B1_WT", nrow(aging_sf3b1WT))

aging3 <- merge_df[, agings]
aging3$SF3B1 <- mutations$SF3B1[match(rownames(aging3), rownames(mutations))]


#### OS of singscores for Aging SF3B1 status in BMMNC-------
load("data/MDS_BMMNC_data/MDS_BMMNC_clinical_df_newCat.Rda")

clinical_bmmnc$PatientId <- rownames(clinical_bmmnc)
sf3b1mut_clinical <- clinical_bmmnc[sf3b1mut,]
sf3b1WT_clinical <- clinical_bmmnc[sf3b1WT,]
sf3b1WT_clinical$Quartile <- rep("SF3B1_WT", nrow(sf3b1WT_clinical))
sf3b1WT_clinical <- sf3b1WT_clinical %>% select(Quartile, OS, Event, PatientId)

mds_sing <- aging3

plot_list <- list()
i <- 1
lowHighP <- list()

for(sc in agings){
  tmp <- as.data.frame(mds_sing) %>% filter(SF3B1 == 1) %>%  mutate(quartile = ntile(get(sc), 2))
  tmp$quartile <- plyr::mapvalues(tmp$quartile, from = c(1, 2), to = c("SF3B1_Mut_Low", "SF3B1_Mut_High"))
  
  
  plot_df <- data.frame(Quartile = tmp$quartile, 
                        OS = sf3b1mut_clinical$OS,
                        Event =  sf3b1mut_clinical$Event,
                        PatientId = sf3b1mut_clinical$PatientId)
  plot_df <- rbind(plot_df, sf3b1WT_clinical)
  
  lowHighP[[sc]] <- pairwise_survdiff(Surv(OS, Event) ~ Quartile,
                                      data = plot_df, p.adjust.method = "none")
  km_infChem <- survfit(formula = Surv(time = as.numeric(plot_df$OS), event = as.numeric(plot_df$Event)) ~ plot_df$Quartile)
  
  g <- ggsurvplot(km_infChem, data = plot_df,
                  conf.int = F, pval = F,
                  fun = function(y) y * 100,
                  legend = "top", xscale = "d_y", break.time.by=365.25*3,
                  legend.title = " ", size = 0.3, censor.size = 1.5,
                  palette = "npg", pval.size = 1.8, pval.coord = c(0.5, 10),
                  legend.labs = c(paste0("SF3B1-Mut,"," High"), paste0("SF3B1-Mut,"," Low"), paste("SF3B1-WT")), 
                  surv.plot.height = 0.3,
                  xlab = "Overall Survival (Years)", ylab="Survival probability (%)")$plot + 
    ggtitle(paste0(sc," signature")) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          #axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.text.y = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.title.x = element_text(size = 5, family = "Helvetica", color = "black"),
          legend.position = "none", aspect.ratio = 1,
          legend.title = element_text(size = 5, family = "Helvetica", color = "black"),
          legend.text = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.title.y = element_text(size = 5, family = "Helvetica", color = "black"),
          plot.title = element_text(size = 5, family = "Helvetica", color = "black", vjust = -2))
  
  ggsave(paste0("figures/OS_disc_SF3B1_Aging_",sc,".pdf"), g, width = 58, height = 44 , units = "mm", dpi = 300)
  
  plot_list[[i]] <- g
  i <- i + 1
}

pdf(file = "figures/OS_disc_KM_by_SF3B1_status_and_3Agings.pdf",
    width = 5.46457,
    height = 6.72441, onefile = T, bg = "white")
ggarrange(plotlist = plot_list, nrow = 3,  align = "hv", common.legend = T)

dev.off()




#### EFS of singscores for Aging SF3B1 status in BMMNC-------
load("data/OneDrive_1_20-09-2023/MDS_BMMNC_clinical_df_newCat.Rda")

clinical_bmmnc$PatientId <- rownames(clinical_bmmnc)
sf3b1mut_clinical <- clinical_bmmnc[sf3b1mut,]
sf3b1WT_clinical <- clinical_bmmnc[sf3b1WT,]
sf3b1WT_clinical$Quartile <- rep("SF3B1_WT", nrow(sf3b1WT_clinical))
sf3b1WT_clinical <- sf3b1WT_clinical %>% select(Quartile, EFS, Event, PatientId)

plot_list_efs <- list()
i <- 1
lowHighP_efs <- list()
for(sc in agings){
  tmp <- as.data.frame(mds_sing) %>% filter(SF3B1 == 1) %>%  mutate(quartile = ntile(get(sc), 2))
  tmp$quartile <- plyr::mapvalues(tmp$quartile, from = c(1, 2), to = c("SF3B1_Mut_Low", "SF3B1_Mut_High"))
  
  
  plot_df <- data.frame(Quartile = tmp$quartile, 
                        EFS = sf3b1mut_clinical$EFS,
                        Event =  sf3b1mut_clinical$Event,
                        PatientId = sf3b1mut_clinical$PatientId)
  plot_df <- rbind(plot_df, sf3b1WT_clinical)
  
  lowHighP_efs[[sc]] <- pairwise_survdiff(Surv(EFS, Event) ~ Quartile,
                                      data = plot_df, p.adjust.method = "none")
  km_infChem <- survfit(formula = Surv(time = as.numeric(plot_df$EFS), event = as.numeric(plot_df$Event)) ~ plot_df$Quartile)
  
  g <- ggsurvplot(km_infChem, data = plot_df,
                  conf.int = F, pval = F,
                  fun = function(y) y * 100,
                  legend = "top", xscale = "d_y", break.time.by=365.25*3,
                  legend.title = " ", size = 0.3, censor.size = 1.5,
                  palette = "npg", pval.size = 1.8, pval.coord = c(0.5, 10),
                  legend.labs = c(paste("SF3B1_Mut_High"), paste("SF3B1_Mut_Low"), paste("SF3B1_WT")), 
                  surv.plot.height = 0.3,
                  xlab = "Event Free Survival (Years)", ylab="Survival probability (%)")$plot + 
    ggtitle(paste0(sc," signature")) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          #axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.text.y = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.title.x = element_text(size = 5, family = "Helvetica", color = "black"),
          legend.position = "none", aspect.ratio = 1,
          legend.title = element_text(size = 5, family = "Helvetica", color = "black"),
          legend.text = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.title.y = element_text(size = 5, family = "Helvetica", color = "black"),
          plot.title = element_text(size = 5, family = "Helvetica", color = "black", vjust = -2))
  
  ggsave(paste0("figures/EFS_disc_SF3B1_Aging_",sc,".pdf"), g, width = 58, height = 44 , units = "mm", dpi = 300)
  
  plot_list_efs[[i]] <- g
  i <- i + 1
}

pdf(file = "figures/EFS_disc_KM_by_SF3B1_status_and_3Agings.pdf",
    width = 3.46457,
    height = 4.72441, onefile = T, bg = "white")
ggarrange(plotlist = plot_list_efs, nrow = 3,  align = "hv", common.legend = T)

dev.off()

# sanity check 
sanit <- merge_df[,agings]
mut <- sanit[sf3b1mut,]
mut$inf_group <- ifelse(mut$Inflammaging > quantile(mut$Inflammaging)[3], "High", "Low")
wt <- sanit[sf3b1WT,]
wt$inf_group <- rep("WT", nrow(wt))

all <- rbind(mut, wt)
all$PatientId <- rownames(all)
clinical_bmmnc$infgroup <- all$inf_group[match(clinical_bmmnc$PatientId, rownames(all))]

san <- clinical_bmmnc %>% filter(infgroup != "High")
fit <- surv_fit(Surv(san$EFS, san$Event) ~ san$infgroup, san)
ggsurvplot(fit, data = san,
           conf.int = F, pval = T,
           fun = function(y) y * 100,
           legend = "top", xscale = "d_y", break.time.by=365.25*3,
           legend.title = " ", size = 0.3, censor.size = 1.5,
           palette = "npg", pval.size = 1.8, pval.coord = c(0.5, 10),
           #legend.labs = c(paste("SF3B1_Mut_High"), paste("SF3B1_Mut_Low"), paste("SF3B1_WT")), 
           surv.plot.height = 0.3,
           xlab = "Years", ylab="Survival probability (%)")$plot 

# CD34 Clinical Outcomes of singscores for Aging based on SF3B1 status in -------
load("data/GSE114922/GSE114922_clinical_df_newCat.Rda")

load("data/merge_validation_for_corr.Rda")
load("data/mutation_validation_for_corr.Rda")

merge_val <- as.data.frame(merge_val)

mutations_val <- as.data.frame(mutations_val)
mutations_val$SF3B1 <- as.factor(mutations_val$SF3B1)
mutations_val$SRSF2 <- as.factor(mutations_val$SRSF2)

sf3b1mut <- mutations_val %>% filter(SF3B1 == 1) %>% rownames()
sf3b1WT <- mutations_val %>% filter(SF3B1 == 0) %>% rownames()

clinical_df <- as.data.frame(t(clinical_df))
clinical_df$Status <- ifelse(clinical_df$Status == "DEAD", 1, 0)
clinical_df$Status <- as.numeric(clinical_df$Status)
clinical_df$Days <- as.numeric(clinical_df$Days)

clinical_df$PatientId <- rownames(clinical_df)
sf3b1mut_clinical <- clinical_df[sf3b1mut,]
sf3b1WT_clinical <- clinical_df[sf3b1WT,]
sf3b1WT_clinical$Quartile <- rep("SF3B1_WT", nrow(sf3b1WT_clinical))
sf3b1WT_clinical <- sf3b1WT_clinical %>% select(Quartile, Days, Status, PatientId)

agings <- c("Aging_Inflammatory chemokines", "Aging_Inflammatory cytokines" )
aging_val <- merge_val[, agings]
aging_val$SF3B1 <- mutations_val$SF3B1[match(rownames(aging_val), rownames(mutations_val))]

colnames(aging_val) <- gsub(".*\\_", "", colnames(aging_val))

mds_sing <- aging_val

plot_list <- list()
i <- 1
lowHighP <- list()

for(sc in colnames(aging_val)[1:2]){
  tmp <- as.data.frame(mds_sing) %>% filter(SF3B1 == 1) %>%  mutate(quartile = ntile(get(sc), 2))
  tmp$quartile <- plyr::mapvalues(tmp$quartile, from = c(1, 2), to = c("SF3B1_Mut_Low", "SF3B1_Mut_High"))
  
  
  plot_df <- data.frame(Quartile = tmp$quartile, 
                        Days = sf3b1mut_clinical$Days,
                        Status =  sf3b1mut_clinical$Status,
                        PatientId = sf3b1mut_clinical$PatientId)
  plot_df <- rbind(plot_df, sf3b1WT_clinical)
  
  lowHighP[[sc]] <- pairwise_survdiff(Surv(Days, Status) ~ Quartile,
                                      data = plot_df, p.adjust.method = "none")
  km_infChem <- survfit(formula = Surv(time = as.numeric(plot_df$Days), event = as.numeric(plot_df$Status)) ~ plot_df$Quartile)
  
  g <- ggsurvplot(km_infChem, data = plot_df,
                  conf.int = F, pval = F,
                  fun = function(y) y * 100,
                  legend = "top", xscale = "d_y", break.time.by=365.25*3,
                  legend.title = " ", size = 0.3, censor.size = 1.5,
                  palette = "npg", pval.size = 1.8, pval.coord = c(0.5, 10),
                  legend.labs = c(paste0("SF3B1-Mut,"," High"), paste0("SF3B1-Mut,"," Low"), paste("SF3B1-WT")), 
                  surv.plot.height = 0.3,
                  xlab = "Overall Survival (Years)", ylab="Survival probability (%)")$plot + 
    ggtitle(paste0(sc," signature")) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          #axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.text.y = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.title.x = element_text(size = 5, family = "Helvetica", color = "black"),
          legend.position = "none", aspect.ratio = 1,
          legend.title = element_text(size = 5, family = "Helvetica", color = "black"),
          legend.text = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.title.y = element_text(size = 5, family = "Helvetica", color = "black"),
          plot.title = element_text(size = 5, family = "Helvetica", color = "black", vjust = -2))
  
  ggsave(paste0("figures/OS_val_SF3B1_Aging_",sc,".pdf"), g, width = 58, height = 44 , units = "mm", dpi = 300)
  
  plot_list[[i]] <- g
  i <- i + 1
}

pdf(file = "figures/OS_val_KM_by_SF3B1_status_and_2Agings.pdf",
    width = 5.46457,
    height = 6.02441, onefile = T, bg = "white")
ggarrange(plotlist = plot_list, nrow = 2,  align = "hv", common.legend = T)

dev.off()




