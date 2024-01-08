###############
## Boxplot script from Nogay
## Plot the distributions absent/present of SF3B1 & SRS2F
## 22 Sep, 2023 
## Sila Gerlevik
################

cell_mutations <- cbind(merge_df, mutations)
colnames(cell_mutations) <- gsub(".*\\_", "", colnames(cell_mutations))
cell_mutations$SF3B1 <- ifelse(cell_mutations$SF3B1 == 0, "SF3B1\nWT", "SF3B1\nMut")
cell_mutations$SRSF2 <- ifelse(cell_mutations$SRSF2 == 0, "SRSF2\nWT", "SRSF2\nMut")
signif_SF3B1_cells <- gsub(".*\\_", "", signif_SF3B1_cells)
signif_SRSF2_cells <- gsub(".*\\_", "", signif_SRSF2_cells)

cell_mutations_v <- cbind(merge_val, mutations_val)
colnames(cell_mutations_v) <- gsub(".*\\_", "", colnames(cell_mutations_v))
cell_mutations_v$SF3B1 <- ifelse(cell_mutations_v$SF3B1 == 0, "SF3B1\nWT", "SF3B1\nMut")
cell_mutations_v$SRSF2 <- ifelse(cell_mutations_v$SRSF2 == 0, "SRSF2\nWT", "SRSF2\nMut")
signif_SF3B1_cells_val <- gsub(".*\\_", "", signif_SF3B1_cells_val)
signif_SRSF2_cells_val <- gsub(".*\\_", "", signif_SRSF2_cells_val)

# summary plot for SF3B1 -----
my_comparisons <- list(c("SF3B1\nWT", "SF3B1\nMut"))
my_Plot_class <- function(signif_SF3B1_cells){
  plot_it <- cell_mutations
  plot_it$SF3B1 <- factor(plot_it$SF3B1, levels = c("SF3B1\nWT", "SF3B1\nMut"))
  ggplot(plot_it, aes(x=SF3B1,y=plot_it[,signif_SF3B1_cells], fill=SF3B1))+
    stat_boxplot(geom = "errorbar",width=0.3,aes(fill=SF3B1))+
    geom_boxplot(notch = F,width=0.8,outlier.shape = NA)+
    scale_fill_manual(values=c("#4DBBD5FF","#E64B35FF"))+
    labs(x = NULL,
         y = signif_SF3B1_cells)+
    geom_jitter(alpha = 0.6, shape=16, position = position_jitter(0.1),color="black",size=0.5)+
    ylim((min(plot_it[,signif_SF3B1_cells])-0.001), (max(plot_it[,signif_SF3B1_cells])+0.08)) +
    stat_compare_means(comparisons = my_comparisons,
                       size = 2.6,
                       label.y = (max(plot_it[,signif_SF3B1_cells]) + 0.001),
                       label.x = 0.3,
                       aes(label = paste0("", after_stat(p.signif))),
                       tip.length = 0) +
    #ggtitle(signif_SF3B1_cells_val) +
    theme(axis.title.x = element_text(size = 5, family = "Helvetica", color = "black"), 
          axis.title.y = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.text.y = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.text.x =  element_text(size = 5, family = "Helvetica", color = "black"), 
          #axis.ticks.y = element_blank(), axis.ticks.x = element_line(), 
          legend.position = "none", legend.title = element_blank(), 
          #legend.text = element_text(size = 5, color = "black"), 
          #legend.key = element_rect(fill = "transparent"), 
          axis.line = element_line(colour = "black"), aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
}

my_Plot_List_class_sf3_d <- lapply(signif_SF3B1_cells, my_Plot_class)


s3_pm <- plot_grid(plotlist = my_Plot_List_class_sf3_d, align = "hv", 
                nrow = 4,ncol = 5)
ggsave("figures/sf3b1_disc_boxplot.pdf", s3_pm, width = 180, units = "mm", dpi = 300)



# summary plot for SF3B1 validation -----
my_comparisons <- list(c("SF3B1\nWT", "SF3B1\nMut"))
my_Plot_class <- function(signif_SF3B1_cells_val){
  plot_it <- cell_mutations_v
  plot_it$SF3B1 <- factor(plot_it$SF3B1, levels = c("SF3B1\nWT", "SF3B1\nMut"))
  ggplot(plot_it, aes(x=SF3B1,y=plot_it[,signif_SF3B1_cells_val], fill=SF3B1))+
    stat_boxplot(geom = "errorbar",width=0.3,aes(fill=SF3B1))+
    geom_boxplot(notch = F,width=0.8,outlier.shape = NA)+
    scale_fill_manual(values=c("#4DBBD5FF","#E64B35FF"))+
    labs(x = NULL,
         y = signif_SF3B1_cells_val)+
    geom_jitter(alpha = 0.6, shape=16, position = position_jitter(0.1),color="black",size=0.5)+
    ylim((min(plot_it[,signif_SF3B1_cells_val])-0.001), (max(plot_it[,signif_SF3B1_cells_val])+0.08)) +
    stat_compare_means(comparisons = my_comparisons,
                       size = 2.6,
                       label.y = (max(plot_it[,signif_SF3B1_cells_val]) + 0.001),
                       label.x = 0.3,
                       aes(label = paste0("", after_stat(p.signif))),
                       tip.length = 0) +
    #ggtitle(signif_SF3B1_cells_val) +
    theme(axis.title.x = element_text(size = 5, family = "Helvetica", color = "black"), 
          axis.title.y = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.text.y = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.text.x =  element_text(size = 5, family = "Helvetica", color = "black"), 
          #axis.ticks.y = element_blank(), axis.ticks.x = element_line(), 
          legend.position = "none", legend.title = element_blank(), 
          #legend.text = element_text(size = 5, color = "black"), 
          #legend.key = element_rect(fill = "transparent"), 
          axis.line = element_line(colour = "black"), aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
}

my_Plot_List_class_sf3_v <- lapply(signif_SF3B1_cells_val, my_Plot_class)


s3_pm_val <- plot_grid(plotlist = my_Plot_List_class_sf3_v, align = "hv", 
                nrow = 4,ncol = 5)
ggsave("figures/sf3b1_val_boxplot.pdf", s3_pm_val, width = 180, units = "mm", dpi = 300)




# summary plot for SRSF2 ----
my_comparisons <- list(c("SRSF2\nWT", "SRSF2\nMut"))
my_Plot_class <- function(signif_SRSF2_cells){
  plot_it <- cell_mutations
  plot_it$SRSF2 <- factor(plot_it$SRSF2, levels = c("SRSF2\nWT", "SRSF2\nMut"))
  ggplot(plot_it, aes(x=SRSF2,y=plot_it[,signif_SRSF2_cells], fill=SRSF2))+
    stat_boxplot(geom = "errorbar",width=0.3,aes(fill=SRSF2))+
    geom_boxplot(notch = F,width=0.8,outlier.shape = NA)+
    scale_fill_manual(values=c("#4DBBD5FF","#E64B35FF"))+
    labs(x = NULL,
         y = signif_SRSF2_cells) +
    geom_jitter(alpha = 0.6, shape=16, position = position_jitter(0.1),color="black",size=0.5)+
    ylim((min(plot_it[,signif_SRSF2_cells])-0.001), (max(plot_it[,signif_SRSF2_cells])+0.08)) +
    stat_compare_means(comparisons = my_comparisons,
                       size = 2.6,
                       label.y = (max(plot_it[,signif_SRSF2_cells]) + 0.001),
                       label.x = 0.3,
                       aes(label = paste0("", after_stat(p.signif))),
                       tip.length = 0) +
    #scale_y_continuous(expand = expansion(mult = c(0,0.3))) +
    theme(axis.title.x = element_text(size = 5, family = "Helvetica", color = "black"), 
          axis.title.y = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.text.y = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.text.x =  element_text(size = 5, family = "Helvetica", color = "black"), 
          #axis.ticks.y = element_blank(), axis.ticks.x = element_line(), 
          legend.position = "none", legend.title = element_blank(), 
          #legend.text = element_text(size = 5, color = "black"), 
          #legend.key = element_rect(fill = "transparent"), 
          axis.line = element_line(colour = "black"), aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
}


my_Plot_List_class_s2_d <- lapply(signif_SRSF2_cells, my_Plot_class)


s2_pm <- plot_grid(plotlist = my_Plot_List_class_s2_d, align = "hv", 
              nrow = 4, ncol = 5)

ggsave("figures/srsf2_disc_boxplot.pdf", s2_pm, width = 180, units = "mm", dpi = 300)


# summary plot for SRSF2 validation ----
my_comparisons <- list(c("SRSF2\nWT", "SRSF2\nMut"))
my_Plot_class <- function(signif_SRSF2_cells_val){
  plot_it <- cell_mutations_v
  plot_it$SRSF2 <- factor(plot_it$SRSF2, levels = c("SRSF2\nWT", "SRSF2\nMut"))
  ggplot(plot_it, aes(x=SRSF2,y=plot_it[,signif_SRSF2_cells_val], fill=SRSF2))+
    stat_boxplot(geom = "errorbar",width=0.3,aes(fill=SRSF2))+
    geom_boxplot(notch = F,width=0.8,outlier.shape = NA)+
    scale_fill_manual(values=c("#4DBBD5FF","#E64B35FF"))+
    labs(x = NULL,
         y = signif_SRSF2_cells_val) +
    geom_jitter(alpha = 0.6, shape=16, position = position_jitter(0.1),color="black",size=0.5)+
    ylim((min(plot_it[,signif_SRSF2_cells_val])-0.001), (max(plot_it[,signif_SRSF2_cells_val])+0.08)) +
    stat_compare_means(comparisons = my_comparisons,
                       size = 2.6,
                       label.y = (max(plot_it[,signif_SRSF2_cells_val]) + 0.001),
                       label.x = 0.3,
                       aes(label = paste0("", after_stat(p.signif))),
                       tip.length = 0) +
    #scale_y_continuous(expand = expansion(mult = c(0,0.3))) +
    theme(axis.title.x = element_text(size = 5, family = "Helvetica", color = "black"), 
          axis.title.y = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.text.y = element_text(size = 5, family = "Helvetica", color = "black"),
          axis.text.x =  element_text(size = 5, family = "Helvetica", color = "black"), 
          #axis.ticks.y = element_blank(), axis.ticks.x = element_line(), 
          legend.position = "none", legend.title = element_blank(), 
          #legend.text = element_text(size = 5, color = "black"), 
          #legend.key = element_rect(fill = "transparent"), 
          axis.line = element_line(colour = "black"), aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
}


my_Plot_List_class_s2_v <- lapply(signif_SRSF2_cells_val, my_Plot_class)


s2_pm_val <- plot_grid(plotlist = my_Plot_List_class_s2_v, align = "hv", 
                nrow = 2, ncol = 3)

ggsave("figures/srsf2_val_boxplot.pdf", s2_pm_val, width = 110, height = 110, units = "mm", dpi = 300)

