################## ***CRONOTIPO FFAs T0*** ################## 
{library(Hmisc)
  library(ggpubr)
  library(ggplot2)
  library(vegan)
  library(pairwiseAdonis)
  library(mixOmics)
  library(ecodist)
  library(reshape2)
  library(readxl)
  library(dunn.test)
}

options(scipen=100)
#setwd('../Cronotipo/Cronotipo_2025')
FFA <- read_excel("Faecal_FFAs.xlsx", sheet=2) 
colnames(FFA)
FFA$Type_diet
FFA[! colnames(FFA) %in% c("ID","Time","Sample_Name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("ID","Time","Sample_Name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Type_diet<-factor(FFA$Type_diet, levels = c("control","morning","evening"))
FFA<-subset.data.frame(FFA, Time=='T0')

SCFA<-FFA[,1:12] # from acetic to valeric 
MCFA<-FFA[,c(1:5,13:21)] # to dodecanoic
LCFA<-FFA[,c(1:5,22:24)] # to Octadecanoic, before the raw total columns
colnames(LCFA)

######################### NORMALIZATIONS ##############################
SCFA_raw<-as.data.frame(SCFA)
MCFA_raw<-as.data.frame(MCFA)
LCFA_raw<- as.data.frame(LCFA)

# #### normalization through relative abundance in each sample, e.g. https://www.nature.com/articles/s41598-021-81420-3
# #### normalization of FFA respect total SCFA, MCFA or LCFA
SCFA_prop<-SCFA
SCFA_prop[ !colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-SCFA_prop[!colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(SCFA_prop[!colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
SCFA_prop<-as.data.frame(SCFA_prop)
# 
MCFA_prop<-MCFA
MCFA_prop[ !colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-MCFA_prop[!colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(MCFA_prop[!colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
MCFA_prop<-as.data.frame(MCFA_prop)
# 
LCFA_prop<-LCFA
LCFA_prop[ !colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-LCFA_prop[!colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(LCFA_prop[!colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
LCFA_prop<-as.data.frame(LCFA_prop)
# rowSums(LCFA[!colnames(LCFA)%in% c("Subject","Condition")])

####################### KRUSKAL WALLIS TESTs T0 ####################
dir.create("FFA")
library(FSA)
##### FFA raw abundances
list <- list("SCFA" = SCFA_raw, "MCFA" = MCFA_raw, "LCFA" = LCFA_raw)
names(list)

for (y in 1:length(list)) {
  temp <- list[[y]]
  results <- NULL
  posthoc_results <- NULL
  
  for (x in colnames(temp)[!colnames(temp) %in% c("ID","Time","Sample_Name","Group","Type_diet")]) {
    # Kruskal-Wallis Test
    test <- kruskal.test(temp[[x]] ~ temp$Type_diet)
    kruskal <- cbind(paste(x), test$statistic, test$p.value)
    results <- rbind.data.frame(results, kruskal)
    
    # Post-hoc Dunn Test
    dunn <- dunnTest(temp[[x]] ~ temp$Type_diet, method = "bh")
    dunn_res <- cbind(variable = paste(x), dunn$res)
    posthoc_results <- rbind.data.frame(posthoc_results, dunn_res)
  }
  
  colnames(results) <- c(paste(names(list)[y]), "Chi-squared", "p.value")
  results$'p.adj(BH)' <- p.adjust(results$p.value, method = "BH")
  row.names(results) <- NULL
  write.csv2(results, file = paste0("FFA/", names(list)[y], "_raw_abund_KW.csv"), row.names = F)
  write.csv2(posthoc_results, file = paste0("FFA/", names(list)[y], "_raw_abund_Dunn.csv"), row.names = F)
  
  table <- melt(data = temp, id.vars = c("ID","Time","Sample_Name","Group","Type_diet"))
  if (y == 1) {
    table_SCFA <- table 
  }
  
  # Plot
  ggplot(data = table, mapping = aes(x = variable, y = value, fill = Type_diet)) +
    scale_fill_manual(values = c("control" = "#CCCCCC", "morning" = "#FFD700", "evening" = '#4682B4')) +
    theme(axis.text.x = element_text(angle = -25, size = 8, vjust = 1, hjust = 0.05)) +
    geom_boxplot() +
    labs(title = paste("Raw", names(list)[y], "quantities", sep = " "), x = "", y = "Quantity")
  
  #ggsave(filename = paste0("FFA/", names(list)[y], "_raw_abundance.png"), height = 3, width = 6, dpi = 300)
  
  if (y == 2) {
    table_MCFA <- table
  }
  if (y == 3) {
    table_LCFA <- table
  }
}

####################### Kruskal-wallis TESTs ####################
###### loop for each FFA
library(stringr)
### Stool
Stool<-FFA
Stool$Type_diet<-trimws(Stool$Type_diet)
Stool$ID<-trimws(Stool$ID)
head(Stool)
Stool<-as.data.frame(Stool)
Stool<-reshape::melt(Stool)
Stool$Type_diet
Stool$variable<-factor(Stool$variable, levels=unique(Stool$variable))
Stool$Type_diet<-factor(Stool$Type_diet, levels = c("control","morning","evening"))
Stool$value<-Stool$value*100 # percentual

for(i in unique(Stool$variable)){
  target<-Stool[Stool$variable==i , ]
  
  test <- dunn.test::dunn.test(target$value, target$Type_diet, method = "bh",list = T) # pairwise
  test <- cbind.data.frame(test$comparison, p_adj_bh=round(test$P.adjusted,3))
  #test$p_adj_bh[test$p_adj_bh==0.000]<-0.001 # if...
  specular_vector<-test$p_adj_bh # can not modify the numeric vector with characters at each step, otherwise every number will became a character itself!
  test$p_adj_bh[specular_vector <=0.001]<-"***"
  test$p_adj_bh[specular_vector >0.001 & specular_vector <=0.01]<-"**"
  test$p_adj_bh[specular_vector >0.01 & specular_vector <=0.05]<-"*"
  
  if(length(which( grepl("*",test$p_adj_bh, fixed=T) )) > 0 ) { # then if there are results (otherwise errors with empty tables!)
    test<-test[ grepl("*",test$p_adj_bh, fixed=T) , ]  # mainteining only significant rows
    comparisons<-unlist(str_split(test$`test$comparison`, pattern ="-")) # it gives a unique vector, not two columns
    comparisons<-gsub(" ","",comparisons)
    # NB: c(T,F) --> only even indexes,   c(F,T) --> only pair indexes   (see the output of unlist-str_split functions, the even indexes are the first pseudo-column elements)
    test$group1<-comparisons[c(T,F)] # 'group1' and 'group2' HAVE to be the col names according to stat_pvalue_manual function
    test$group2<-comparisons[c(F,T)]
    # test # a 'manual' check --> OK
    test<-test[ seq(length(test$`test$comparison`),1,-1) , ] # inverting the order of comparisons to enhance the plot aesthetics (the last row here is the one with most width), moreover using and inverse seq instead of c(3,2,1) to adapt each loop to the row numbers (=number of signif results)
    test$`test$comparison`<-NULL
    
    kruskal_res<-kruskal.test(target$value~target$Type_diet)
    kruskal_res<-round(kruskal_res$p.value, 3)
    kruskal_res[kruskal_res==0]<-0.001 # if...
    
    ggboxplot(data=target, x="Type_diet", y="value", fill="Type_diet", 
              width=0.6, size= 0.15, alpha= 0.2, outlier.shape = NA) +
      stat_pvalue_manual(test, label="p_adj_bh",
                         y.position = max(target$value),
                         step.increase = 0.065,
                         size = 2.5,
                         bracket.size = 0.15
      ) +
      scale_fill_manual(values=c("control" = "#CCCCCC", "morning" = "#FFD700", "evening" = '#4682B4')) +
      scale_color_manual(values=c("control" = "#CCCCCC", "morning" = "#FFD700", "evening" = '#4682B4')) +
      geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.6, jitter.width = 0.9),
                 aes(color=Type_diet), size= 0.3, alpha= 0.5) +  
      theme_classic2(base_size = 7) + 
      theme(strip.text.x=element_text(size=12,colour="black"),
            axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=7), 
            axis.text.y = element_text(size=4.5),
            plot.margin = margin(5,2,0,2),
            plot.title= element_text(size=10, hjust = 0.5, vjust=1.8) ,
            plot.subtitle = element_text(size=7.4, vjust=1) ,
            legend.key.size=unit(0.8,"cm"), 
            legend.text=element_text(size=7),
            panel.grid.major.y = element_line(size=0.12, color="gray"),
            panel.grid.minor.y = element_line(size=0.03, color="gray"),
            panel.grid.major.x = element_blank()
      ) +
      guides( color="none", fill="none" ) +
      labs(y="Percentage (%)", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("FFA/Serum_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }
  
  if(length(which( grepl("*",test$p_adj_bh, fixed=T) )) == 0 ) { # otherwise, if there are not pairwise results...
    
    kruskal_res<-kruskal.test(target$value~target$Type_diet)
    kruskal_res<-round(kruskal_res$p.value, 3)
    kruskal_res[kruskal_res==0]<-0.001 # if...
    
    ggboxplot(data=target, x="Type_diet", y="value", fill="Type_diet", 
              width=0.6, size= 0.15, alpha= 0.2, outlier.shape = NA) +
      scale_fill_manual(values=c("control" = "#CCCCCC", "morning" = "#FFD700", "evening" = '#4682B4')) +
      scale_color_manual(values=c("control" = "#CCCCCC", "morning" = "#FFD700", "evening" = '#4682B4')) +
      geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.6, jitter.width = 0.9),
                 aes(color=Type_diet), size= 0.3, alpha= 0.5) +  
      theme_classic2(base_size = 7) + 
      theme(strip.text.x=element_text(size=12,colour="black"),
            axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=7), 
            axis.text.y = element_text(size=5),
            plot.margin = margin(5,2,0,2),
            plot.title= element_text(size=10, hjust = 0.5, vjust=1.8) ,
            plot.subtitle = element_text(size=7.4, vjust=1) ,
            legend.key.size=unit(0.8,"cm"), 
            legend.text=element_text(size=7),
            panel.grid.major.y = element_line(size=0.12, color="gray"),
            panel.grid.minor.y = element_line(size=0.03, color="gray"),
            panel.grid.major.x = element_blank()
      ) +
      guides( color="none", fill="none" ) +
      labs(y="Percentage (%)", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("FFA/Serum_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }
  
}

########### PCoA ########### 
list_PCoA <- list("SCFA" = SCFA_prop, "MCFA" = MCFA_prop, "LCFA" = LCFA_prop)
p_values <- NULL
for (x in 1:length(list_PCoA)) {
  rm(data, coord)
  data <- as.data.frame(list_PCoA[[x]])
  row.names(data) <- data$ID
  
  if (x == 3) {
    data <- data[!is.na(data$Octadecanoic),] 
  }
  
  dist <- vegdist(data[, !colnames(data) %in% c("ID","Time","Sample_Name","Group","Type_diet")], method = "bray") # calcolare dissimilarità Bray-Curtis
  coord <- cmdscale(dist, k = 2, eig = TRUE) # PCoA
  
  coord_df <- as.data.frame(coord$points)
  colnames(coord_df) <- c("PC1", "PC2")
  coord_df$ID <- row.names(coord_df)
  coord_df$Type_diet <- data$Type_diet
  
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Type_diet)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("control" = alpha("#CCCCCC", 0.5), "morning" = alpha("#FFD700", 0.5), "evening"=alpha('4682B4',0.5))) +
    stat_ellipse(size = 0.2) +
    geom_text(aes(label = ID), size = 2, color = "black") +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of different type of diet\n with Bray-Curtis dissimilarity")
    )
  
  # again but with no names
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Type_diet)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("control" = alpha("#CCCCCC", 0.5), "morning" = alpha("#FFD700", 0.5), "evening"=alpha('#4682B4',0.5))) +
    stat_ellipse(size = 0.2) +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of different type of diet\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("FFA/PCoA_BRAY_", names(list_PCoA)[x], "T0_noname.png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )                         
  
  p_value_disp <- permutest(betadisper(dist, group = data$Type_diet), permutations = 9999)
  p_value_diver <- adonis(dist ~ Type_diet, data = data, permutations = 9999)
  p_value_disp <- as.data.frame(p_value_disp$tab[1, ])
  colnames(p_value_disp) <- colnames(p_value_diver$aov.tab[1, ])
  new_values <- rbind(p_value_diver$aov.tab[1, ], p_value_disp)
  row.names(new_values) <- c(paste("Beta_diversity", names(list_PCoA)[x], sep = "_"),
                             paste("Beta_dispersion", names(list_PCoA)[x], sep = "_"))

  p_values <- rbind(p_values, new_values)
  
}

p_values$adjusted_BH_p_values<-p.adjust(p_values$`Pr(>F)`, method = "BH")
p_values
write.csv2(p_values, file="FFA/p_values_BRAY_T0.csv", quote = F,row.names = T)


########### total scfa, mcfa, lcfa #############
FFA <- read_excel("Faecal_FFAs.xlsx", sheet=3) 
colnames(FFA)
FFA$Type_diet
FFA[! colnames(FFA) %in% c("ID","Time","Sample_Name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("ID","Time","Sample_Name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Type_diet<-factor(FFA$Type_diet, levels = c("control","morning","evening"))
FFA<-subset.data.frame(FFA, Time=='T0')

SCFA<-FFA[,1:6] 
MCFA<-FFA[,c(1:5,7)] 
LCFA<-FFA[,c(1:5,8)] 
colnames(LCFA)

######################### NORMALIZATIONS ##############################
SCFA_raw<-SCFA
MCFA_raw<-MCFA
LCFA_raw<- LCFA

####################### Kruskal-wallis TESTs ####################
###### loop for each FFA
library(stringr)

### Stool
Stool<-FFA
Stool$Type_diet<-trimws(Stool$Type_diet)
Stool$ID<-trimws(Stool$ID)
head(Stool)
Stool<-as.data.frame(Stool)
Stool<-reshape::melt(Stool)
Stool$Type_diet
Stool$variable<-factor(Stool$variable, levels=unique(Stool$variable))
Stool$Type_diet<-factor(Stool$Type_diet, levels = c("control","morning","evening"))
Stool$value<-Stool$value*100 # percentual

for(i in unique(Stool$variable)){
  target<-Stool[Stool$variable==i , ]
  
  test <- dunn.test::dunn.test(target$value, target$Type_diet, method = "bh",list = T) # pairwise
  test <- cbind.data.frame(test$comparison, p_adj_bh=round(test$P.adjusted,3))
  #test$p_adj_bh[test$p_adj_bh==0.000]<-0.001 # if...
  specular_vector<-test$p_adj_bh # can not modify the numeric vector with characters at each step, otherwise every number will became a character itself!
  test$p_adj_bh[specular_vector <=0.001]<-"***"
  test$p_adj_bh[specular_vector >0.001 & specular_vector <=0.01]<-"**"
  test$p_adj_bh[specular_vector >0.01 & specular_vector <=0.05]<-"*"
  
  if(length(which( grepl("*",test$p_adj_bh, fixed=T) )) > 0 ) { # then if there are results (otherwise errors with empty tables!)
    test<-test[ grepl("*",test$p_adj_bh, fixed=T) , ]  # mainteining only significant rows
    comparisons<-unlist(str_split(test$`test$comparison`, pattern ="-")) # it gives a unique vector, not two columns
    comparisons<-gsub(" ","",comparisons)
    # NB: c(T,F) --> only even indexes,   c(F,T) --> only pair indexes   (see the output of unlist-str_split functions, the even indexes are the first pseudo-column elements)
    test$group1<-comparisons[c(T,F)] # 'group1' and 'group2' HAVE to be the col names according to stat_pvalue_manual function
    test$group2<-comparisons[c(F,T)]
    # test # a 'manual' check --> OK
    test<-test[ seq(length(test$`test$comparison`),1,-1) , ] # inverting the order of comparisons to enhance the plot aesthetics (the last row here is the one with most width), moreover using and inverse seq instead of c(3,2,1) to adapt each loop to the row numbers (=number of signif results)
    test$`test$comparison`<-NULL
    
    kruskal_res<-kruskal.test(target$value~target$Type_diet)
    kruskal_res<-round(kruskal_res$p.value, 3)
    kruskal_res[kruskal_res==0]<-0.001 # if...
    
    ggboxplot(data=target, x="Type_diet", y="value", fill="Type_diet", 
              width=0.6, size= 0.15, alpha= 0.2, outlier.shape = NA) +
      stat_pvalue_manual(test, label="p_adj_bh",
                         y.position = max(target$value),
                         step.increase = 0.065,
                         size = 2.5,
                         bracket.size = 0.15
      ) +
      scale_fill_manual(values=c("control" = "#CCCCCC", "morning" = "#FFD700", "evening" = '#4682B4')) +
      scale_color_manual(values=c("control" = "#CCCCCC", "morning" = "#FFD700", "evening" = '#4682B4')) +
      geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.6, jitter.width = 0.9),
                 aes(color=Type_diet), size= 0.3, alpha= 0.5) +  
      theme_classic2(base_size = 7) + 
      theme(strip.text.x=element_text(size=12,colour="black"),
            axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=7), 
            axis.text.y = element_text(size=4.5),
            plot.margin = margin(5,2,0,2),
            plot.title= element_text(size=10, hjust = 0.5, vjust=1.8) ,
            plot.subtitle = element_text(size=7.4, vjust=1) ,
            legend.key.size=unit(0.8,"cm"), 
            legend.text=element_text(size=7),
            panel.grid.major.y = element_line(size=0.12, color="gray"),
            panel.grid.minor.y = element_line(size=0.03, color="gray"),
            panel.grid.major.x = element_blank()
      ) +
      guides( color="none", fill="none" ) +
      labs(y="Percentage (%)", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("FFA/Serum_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }
  
  if(length(which( grepl("*",test$p_adj_bh, fixed=T) )) == 0 ) { # otherwise, if there are not pairwise results...
    
    kruskal_res<-kruskal.test(target$value~target$Type_diet)
    kruskal_res<-round(kruskal_res$p.value, 3)
    kruskal_res[kruskal_res==0]<-0.001 # if...
    
    ggboxplot(data=target, x="Type_diet", y="value", fill="Type_diet", 
              width=0.6, size= 0.15, alpha= 0.2, outlier.shape = NA) +
      scale_fill_manual(values=c("control" = "#CCCCCC", "morning" = "#FFD700", "evening" = '#4682B4')) +
      scale_color_manual(values=c("control" = "#CCCCCC", "morning" = "#FFD700", "evening" = '#4682B4')) +
      geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.6, jitter.width = 0.9),
                 aes(color=Type_diet), size= 0.3, alpha= 0.5) +  
      theme_classic2(base_size = 7) + 
      theme(strip.text.x=element_text(size=12,colour="black"),
            axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=7), 
            axis.text.y = element_text(size=5),
            plot.margin = margin(5,2,0,2),
            plot.title= element_text(size=10, hjust = 0.5, vjust=1.8) ,
            plot.subtitle = element_text(size=7.4, vjust=1) ,
            legend.key.size=unit(0.8,"cm"), 
            legend.text=element_text(size=7),
            panel.grid.major.y = element_line(size=0.12, color="gray"),
            panel.grid.minor.y = element_line(size=0.03, color="gray"),
            panel.grid.major.x = element_blank()
      ) +
      guides( color="none", fill="none" ) +
      labs(y="Percentage (%)", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("FFA/Faecal_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }
  
}

################## ***CRONOTIPO FFAs T1*** ################## 
options(scipen=100)
#setwd('../Cronotipo/Cronotipo_2025')
FFA <- read_excel("Faecal_FFAs.xlsx", sheet=2) 
colnames(FFA)
FFA$Type_diet
FFA[! colnames(FFA) %in% c("ID","Time","Sample_Name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("ID","Time","Sample_Name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Type_diet<-factor(FFA$Type_diet, levels = c("control","morning","evening"))
FFA<-subset.data.frame(FFA, Time=='T1')

SCFA<-FFA[,1:12] # from acetic to valeric 
MCFA<-FFA[,c(1:5,13:21)] # to dodecanoic
LCFA<-FFA[,c(1:5,22:24)] # to Octadecanoic, before the raw total columns
colnames(LCFA)

######################### NORMALIZATIONS ##############################
SCFA_raw<-as.data.frame(SCFA)
MCFA_raw<-as.data.frame(MCFA)
LCFA_raw<- as.data.frame(LCFA)

# #### normalization through relative abundance in each sample, e.g. https://www.nature.com/articles/s41598-021-81420-3
# #### normalization of FFA respect total SCFA, MCFA or LCFA
SCFA_prop<-SCFA
SCFA_prop[ !colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-SCFA_prop[!colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(SCFA_prop[!colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
SCFA_prop<-as.data.frame(SCFA_prop)
# 
MCFA_prop<-MCFA
MCFA_prop[ !colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-MCFA_prop[!colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(MCFA_prop[!colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
MCFA_prop<-as.data.frame(MCFA_prop)
# 
LCFA_prop<-LCFA
LCFA_prop[ !colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-LCFA_prop[!colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(LCFA_prop[!colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
LCFA_prop<-as.data.frame(LCFA_prop)
# rowSums(LCFA[!colnames(LCFA)%in% c("Subject","Condition")])

####################### KRUSKAL WALLIS TESTs T1 ####################
#dir.create("FFA")
library(FSA)
##### FFA raw abundances
list <- list("SCFA" = SCFA_raw, "MCFA" = MCFA_raw, "LCFA" = LCFA_raw)
names(list)

for (y in 1:length(list)) {
  temp <- list[[y]]
  results <- NULL
  posthoc_results <- NULL
  
  for (x in colnames(temp)[!colnames(temp) %in% c("ID","Time","Sample_Name","Group","Type_diet")]) {
    # Kruskal-Wallis Test
    test <- kruskal.test(temp[[x]] ~ temp$Type_diet)
    kruskal <- cbind(paste(x), test$statistic, test$p.value)
    results <- rbind.data.frame(results, kruskal)
    
    # Post-hoc Dunn Test
    dunn <- dunnTest(temp[[x]] ~ temp$Type_diet, method = "bh")
    dunn_res <- cbind(variable = paste(x), dunn$res)
    posthoc_results <- rbind.data.frame(posthoc_results, dunn_res)
  }
  
  # Nomi delle colonne
  colnames(results) <- c(paste(names(list)[y]), "Chi-squared", "p.value")
  results$'p.adj(BH)' <- p.adjust(results$p.value, method = "BH")
  row.names(results) <- NULL
  
  write.csv2(results, file = paste0("FFA/", names(list)[y], "_raw_abund_KW.csv"), row.names = F)
  write.csv2(posthoc_results, file = paste0("FFA/", names(list)[y], "_raw_abund_Dunn.csv"), row.names = F)
  table <- melt(data = temp, id.vars = c("ID","Time","Sample_Name","Group","Type_diet"))
  if (y == 1) {
    table_SCFA <- table 
  }
  
  # Plot
  ggplot(data = table, mapping = aes(x = variable, y = value, fill = Type_diet)) +
    scale_fill_manual(values = c("control" = "#555555", "morning" = "#FFB347", "evening" = '#0047AB')) +
    theme(axis.text.x = element_text(angle = -25, size = 8, vjust = 1, hjust = 0.05)) +
    geom_boxplot() +
    labs(title = paste("Raw", names(list)[y], "quantities", sep = " "), x = "", y = "Quantity")
  
  #ggsave(filename = paste0("FFA/", names(list)[y], "_raw_abundance.png"), height = 3, width = 6, dpi = 300)
  
  if (y == 2) {
    table_MCFA <- table
  }
  if (y == 3) {
    table_LCFA <- table
  }
}

####################### Kruskal-wallis TESTs ####################
###### loop for each FFA
library(stringr)
### Stool
Stool<-FFA
Stool$Type_diet<-trimws(Stool$Type_diet)
Stool$ID<-trimws(Stool$ID)
head(Stool)
Stool<-as.data.frame(Stool)
Stool<-reshape::melt(Stool)
Stool$Type_diet
Stool$variable<-factor(Stool$variable, levels=unique(Stool$variable))
Stool$Type_diet<-factor(Stool$Type_diet, levels = c("control","morning","evening"))
Stool$value<-Stool$value*100 # percentual

for(i in unique(Stool$variable)){
  target<-Stool[Stool$variable==i , ]
  
  test <- dunn.test::dunn.test(target$value, target$Type_diet, method = "bh",list = T) # pairwise
  test <- cbind.data.frame(test$comparison, p_adj_bh=round(test$P.adjusted,3))
  #test$p_adj_bh[test$p_adj_bh==0.000]<-0.001 # if...
  specular_vector<-test$p_adj_bh # can not modify the numeric vector with characters at each step, otherwise every number will became a character itself!
  test$p_adj_bh[specular_vector <=0.001]<-"***"
  test$p_adj_bh[specular_vector >0.001 & specular_vector <=0.01]<-"**"
  test$p_adj_bh[specular_vector >0.01 & specular_vector <=0.05]<-"*"
  
  if(length(which( grepl("*",test$p_adj_bh, fixed=T) )) > 0 ) { # then if there are results (otherwise errors with empty tables!)
    test<-test[ grepl("*",test$p_adj_bh, fixed=T) , ]  # mainteining only significant rows
    comparisons<-unlist(str_split(test$`test$comparison`, pattern ="-")) # it gives a unique vector, not two columns
    comparisons<-gsub(" ","",comparisons)
    # NB: c(T,F) --> only even indexes,   c(F,T) --> only pair indexes   (see the output of unlist-str_split functions, the even indexes are the first pseudo-column elements)
    test$group1<-comparisons[c(T,F)] # 'group1' and 'group2' HAVE to be the col names according to stat_pvalue_manual function
    test$group2<-comparisons[c(F,T)]
    # test # a 'manual' check --> OK
    test<-test[ seq(length(test$`test$comparison`),1,-1) , ] # inverting the order of comparisons to enhance the plot aesthetics (the last row here is the one with most width), moreover using and inverse seq instead of c(3,2,1) to adapt each loop to the row numbers (=number of signif results)
    test$`test$comparison`<-NULL
    
    kruskal_res<-kruskal.test(target$value~target$Type_diet)
    kruskal_res<-round(kruskal_res$p.value, 3)
    kruskal_res[kruskal_res==0]<-0.001 # if...
    
    ggboxplot(data=target, x="Type_diet", y="value", fill="Type_diet", 
              width=0.6, size= 0.15, alpha= 0.2, outlier.shape = NA) +
      stat_pvalue_manual(test, label="p_adj_bh",
                         y.position = max(target$value),
                         step.increase = 0.065,
                         size = 2.5,
                         bracket.size = 0.15
      ) +
      scale_fill_manual(values= c("control" = "#555555", "morning" = "#FFB347", "evening" = '#0047AB')) +
      scale_color_manual(values= c("control" = "#555555", "morning" = "#FFB347", "evening" = '#0047AB')) +
      geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.6, jitter.width = 0.9),
                 aes(color=Type_diet), size= 0.3, alpha= 0.5) +  
      theme_classic2(base_size = 7) + 
      theme(strip.text.x=element_text(size=12,colour="black"),
            axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=7), 
            axis.text.y = element_text(size=4.5),
            plot.margin = margin(5,2,0,2),
            plot.title= element_text(size=10, hjust = 0.5, vjust=1.8) ,
            plot.subtitle = element_text(size=7.4, vjust=1) ,
            legend.key.size=unit(0.8,"cm"), 
            legend.text=element_text(size=7),
            panel.grid.major.y = element_line(size=0.12, color="gray"),
            panel.grid.minor.y = element_line(size=0.03, color="gray"),
            panel.grid.major.x = element_blank()
      ) +
      guides( color="none", fill="none" ) +
      labs(y="Percentage (%)", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("FFA/Serum_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }
  
  if(length(which( grepl("*",test$p_adj_bh, fixed=T) )) == 0 ) { # otherwise, if there are not pairwise results...
    
    kruskal_res<-kruskal.test(target$value~target$Type_diet)
    kruskal_res<-round(kruskal_res$p.value, 3)
    kruskal_res[kruskal_res==0]<-0.001 # if...
    
    ggboxplot(data=target, x="Type_diet", y="value", fill="Type_diet", 
              width=0.6, size= 0.15, alpha= 0.2, outlier.shape = NA) +
      scale_fill_manual(values= c("control" = "#555555", "morning" = "#FFB347", "evening" = '#0047AB')) +
      scale_color_manual(values= c("control" = "#555555", "morning" = "#FFB347", "evening" = '#0047AB')) +
      geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.6, jitter.width = 0.9),
                 aes(color=Type_diet), size= 0.3, alpha= 0.5) +  
      theme_classic2(base_size = 7) + 
      theme(strip.text.x=element_text(size=12,colour="black"),
            axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=7), 
            axis.text.y = element_text(size=5),
            plot.margin = margin(5,2,0,2),
            plot.title= element_text(size=10, hjust = 0.5, vjust=1.8) ,
            plot.subtitle = element_text(size=7.4, vjust=1) ,
            legend.key.size=unit(0.8,"cm"), 
            legend.text=element_text(size=7),
            panel.grid.major.y = element_line(size=0.12, color="gray"),
            panel.grid.minor.y = element_line(size=0.03, color="gray"),
            panel.grid.major.x = element_blank()
      ) +
      guides( color="none", fill="none" ) +
      labs(y="Percentage (%)", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("FFA/Serum_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }
  
}

########### PCoA ########### 
list_PCoA <- list("SCFA" = SCFA_prop, "MCFA" = MCFA_prop, "LCFA" = LCFA_prop)
p_values <- NULL
for (x in 1:length(list_PCoA)) {
  rm(data, coord)
  data <- as.data.frame(list_PCoA[[x]])
  row.names(data) <- data$ID
  
  if (x == 3) {
    data <- data[!is.na(data$Octadecanoic),] 
  }
  
  dist <- vegdist(data[, !colnames(data) %in% c("ID","Time","Sample_Name","Group","Type_diet")], method = "bray") # calcolare dissimilarità Bray-Curtis
  coord <- cmdscale(dist, k = 2, eig = TRUE) # PCoA
  
  coord_df <- as.data.frame(coord$points)
  colnames(coord_df) <- c("PC1", "PC2")
  coord_df$ID <- row.names(coord_df)
  coord_df$Type_diet <- data$Type_diet
  
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Type_diet)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("control" = alpha("#555555", 0.5), "morning" = alpha("#FFB347", 0.5), "evening"=alpha('#0047AB',0.5))) +
    stat_ellipse(size = 0.2) +
    geom_text(aes(label = ID), size = 2, color = "black") +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of different type of diet\n with Bray-Curtis dissimilarity")
    )
  
  # again but with no names
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Type_diet)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("control" = alpha("#555555", 0.5), "morning" = alpha("#FFB347", 0.5), "evening"=alpha('#0047AB',0.5))) +
    stat_ellipse(size = 0.2) +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of different type of diet\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("FFA/PCoA_BRAY_", names(list_PCoA)[x], "T1_noname.png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )                         
  
  p_value_disp <- permutest(betadisper(dist, group = data$Type_diet), permutations = 9999)
  p_value_diver <- adonis(dist ~ Type_diet, data = data, permutations = 9999)
  p_value_disp <- as.data.frame(p_value_disp$tab[1, ])
  colnames(p_value_disp) <- colnames(p_value_diver$aov.tab[1, ])
  new_values <- rbind(p_value_diver$aov.tab[1, ], p_value_disp)
  row.names(new_values) <- c(paste("Beta_diversity", names(list_PCoA)[x], sep = "_"),
                             paste("Beta_dispersion", names(list_PCoA)[x], sep = "_"))
  
  p_values <- rbind(p_values, new_values)
  
}

p_values$adjusted_BH_p_values<-p.adjust(p_values$`Pr(>F)`, method = "BH")
p_values
write.csv2(p_values, file="FFA/p_values_BRAY_T1.csv", quote = F,row.names = T)


########### total scfa, mcfa, lcfa #############
FFA <- read_excel("Faecal_FFAs.xlsx", sheet=3) 
colnames(FFA)
FFA$Type_diet
FFA[! colnames(FFA) %in% c("ID","Time","Sample_Name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("ID","Time","Sample_Name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Type_diet<-factor(FFA$Type_diet, levels = c("control","morning","evening"))
FFA<-subset.data.frame(FFA, Time=='T1')

SCFA<-FFA[,1:6] 
MCFA<-FFA[,c(1:5,7)] 
LCFA<-FFA[,c(1:5,8)] 
colnames(LCFA)

######################### NORMALIZATIONS ##############################
SCFA_raw<-SCFA
MCFA_raw<-MCFA
LCFA_raw<- LCFA

####################### Kruskal-wallis TESTs ####################
###### loop for each FFA
library(stringr)

### Stool
Stool<-FFA
Stool$Type_diet<-trimws(Stool$Type_diet)
Stool$ID<-trimws(Stool$ID)
head(Stool)
Stool<-as.data.frame(Stool)
Stool<-reshape::melt(Stool)
Stool$Type_diet
Stool$variable<-factor(Stool$variable, levels=unique(Stool$variable))
Stool$Type_diet<-factor(Stool$Type_diet, levels = c("control","morning","evening"))
Stool$value<-Stool$value*100 # percentual

for(i in unique(Stool$variable)){
  target<-Stool[Stool$variable==i , ]
  
  test <- dunn.test::dunn.test(target$value, target$Type_diet, method = "bh",list = T) # pairwise
  test <- cbind.data.frame(test$comparison, p_adj_bh=round(test$P.adjusted,3))
  #test$p_adj_bh[test$p_adj_bh==0.000]<-0.001 # if...
  specular_vector<-test$p_adj_bh # can not modify the numeric vector with characters at each step, otherwise every number will became a character itself!
  test$p_adj_bh[specular_vector <=0.001]<-"***"
  test$p_adj_bh[specular_vector >0.001 & specular_vector <=0.01]<-"**"
  test$p_adj_bh[specular_vector >0.01 & specular_vector <=0.05]<-"*"
  
  if(length(which( grepl("*",test$p_adj_bh, fixed=T) )) > 0 ) { # then if there are results (otherwise errors with empty tables!)
    test<-test[ grepl("*",test$p_adj_bh, fixed=T) , ]  # mainteining only significant rows
    comparisons<-unlist(str_split(test$`test$comparison`, pattern ="-")) # it gives a unique vector, not two columns
    comparisons<-gsub(" ","",comparisons)
    # NB: c(T,F) --> only even indexes,   c(F,T) --> only pair indexes   (see the output of unlist-str_split functions, the even indexes are the first pseudo-column elements)
    test$group1<-comparisons[c(T,F)] # 'group1' and 'group2' HAVE to be the col names according to stat_pvalue_manual function
    test$group2<-comparisons[c(F,T)]
    # test # a 'manual' check --> OK
    test<-test[ seq(length(test$`test$comparison`),1,-1) , ] # inverting the order of comparisons to enhance the plot aesthetics (the last row here is the one with most width), moreover using and inverse seq instead of c(3,2,1) to adapt each loop to the row numbers (=number of signif results)
    test$`test$comparison`<-NULL
    
    kruskal_res<-kruskal.test(target$value~target$Type_diet)
    kruskal_res<-round(kruskal_res$p.value, 3)
    kruskal_res[kruskal_res==0]<-0.001 # if...
    
    ggboxplot(data=target, x="Type_diet", y="value", fill="Type_diet", 
              width=0.6, size= 0.15, alpha= 0.2, outlier.shape = NA) +
      stat_pvalue_manual(test, label="p_adj_bh",
                         y.position = max(target$value),
                         step.increase = 0.065,
                         size = 2.5,
                         bracket.size = 0.15
      ) +
      scale_fill_manual(values=c("control" = "#555555", "morning" = "#FFB347", "evening" = '#0047AB')) +
      scale_color_manual(values=c("control" = "#555555", "morning" = "#FFB347", "evening" = '#0047AB')) +
      geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.6, jitter.width = 0.9),
                 aes(color=Type_diet), size= 0.3, alpha= 0.5) +  
      theme_classic2(base_size = 7) + 
      theme(strip.text.x=element_text(size=12,colour="black"),
            axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=7), 
            axis.text.y = element_text(size=4.5),
            plot.margin = margin(5,2,0,2),
            plot.title= element_text(size=10, hjust = 0.5, vjust=1.8) ,
            plot.subtitle = element_text(size=7.4, vjust=1) ,
            legend.key.size=unit(0.8,"cm"), 
            legend.text=element_text(size=7),
            panel.grid.major.y = element_line(size=0.12, color="gray"),
            panel.grid.minor.y = element_line(size=0.03, color="gray"),
            panel.grid.major.x = element_blank()
      ) +
      guides( color="none", fill="none" ) +
      labs(y="Percentage (%)", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("FFA/Serum_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }
  
  if(length(which( grepl("*",test$p_adj_bh, fixed=T) )) == 0 ) { # otherwise, if there are not pairwise results...
    
    kruskal_res<-kruskal.test(target$value~target$Type_diet)
    kruskal_res<-round(kruskal_res$p.value, 3)
    kruskal_res[kruskal_res==0]<-0.001 # if...
    
    ggboxplot(data=target, x="Type_diet", y="value", fill="Type_diet", 
              width=0.6, size= 0.15, alpha= 0.2, outlier.shape = NA) +
      scale_fill_manual(values=c("control" = "#555555", "morning" = "#FFB347", "evening" = '#0047AB')) +
      scale_color_manual(values=c("control" = "#555555", "morning" = "#FFB347", "evening" = '#0047AB')) +
      geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.6, jitter.width = 0.9),
                 aes(color=Type_diet), size= 0.3, alpha= 0.5) +  
      theme_classic2(base_size = 7) + 
      theme(strip.text.x=element_text(size=12,colour="black"),
            axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=7), 
            axis.text.y = element_text(size=5),
            plot.margin = margin(5,2,0,2),
            plot.title= element_text(size=10, hjust = 0.5, vjust=1.8) ,
            plot.subtitle = element_text(size=7.4, vjust=1) ,
            legend.key.size=unit(0.8,"cm"), 
            legend.text=element_text(size=7),
            panel.grid.major.y = element_line(size=0.12, color="gray"),
            panel.grid.minor.y = element_line(size=0.03, color="gray"),
            panel.grid.major.x = element_blank()
      ) +
      guides( color="none", fill="none" ) +
      labs(y="Percentage (%)", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("FFA/Faecal_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }
  
}

########## ***TYPE OF DIET*** ##########
########## MORNING ##########
options(scipen=100)
FFA <- read_excel("Faecal_FFAs.xlsx", sheet=2) 
colnames(FFA)
FFA$Type_diet
FFA[! colnames(FFA) %in% c("ID","Time","Sample_Name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("ID","Time","Sample_Name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Type_diet<-factor(FFA$Type_diet, levels = c("control","morning","evening"))
FFA<-subset.data.frame(FFA, Type_diet=='morning')
exclude_cols <- c("ID", "Time", "Sample_Name", "Group", "Type_diet")
FFA[!colnames(FFA) %in% exclude_cols] <- 
  apply(FFA[!colnames(FFA) %in% exclude_cols], 2, function(x) as.numeric(x) * 100)

SCFA<-FFA[,1:12] # from acetic to valeric 
MCFA<-FFA[,c(1:5,13:21)] # to dodecanoic
LCFA<-FFA[,c(1:5,22:24)] # to Octadecanoic, before the raw total columns
colnames(LCFA)

######################### NORMALIZATIONS ##############################
SCFA_raw<-as.data.frame(SCFA)
MCFA_raw<-as.data.frame(MCFA)
LCFA_raw<- as.data.frame(LCFA)

# #### normalization through relative abundance in each sample, e.g. https://www.nature.com/articles/s41598-021-81420-3
# #### normalization of FFA respect total SCFA, MCFA or LCFA
SCFA_prop<-SCFA
SCFA_prop[ !colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-SCFA_prop[!colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(SCFA_prop[!colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
SCFA_prop<-as.data.frame(SCFA_prop)
# 
MCFA_prop<-MCFA
MCFA_prop[ !colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-MCFA_prop[!colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(MCFA_prop[!colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
MCFA_prop<-as.data.frame(MCFA_prop)
# 
LCFA_prop<-LCFA
LCFA_prop[ !colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-LCFA_prop[!colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(LCFA_prop[!colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
LCFA_prop<-as.data.frame(LCFA_prop)
# rowSums(LCFA[!colnames(LCFA)%in% c("Subject","Condition")])

####################### MANN WITHNEY TESTs ####################
#dir.create("FFA")

##### FFA raw abundances
list<-list("SCFA"=SCFA_raw, "MCFA"=MCFA_raw, "LCFA"=LCFA_raw)
names(list)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("ID","Time","Sample_Name","Group","Type_diet")]) {
    test<-wilcox.test(temp[[x]]~temp$Time)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  write.csv2(results, file=paste0("FFA/",names(list)[y],"_raw_abund_Mann_Whit.csv"), row.names = F)
  
  table <- melt(data = temp, id.vars = c("ID","Time","Sample_Name","Group","Type_diet"))
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Time)) +
    scale_fill_manual(values =c("T0" = "#FFD700", "T1" = "#FFB347")) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  #ggsave(filename = paste0("FFA/",names(list)[y],"_raw_abundance.png"), height=3, width = 6, dpi=300)
  
  if(y==2){ # if MCFA
    table_MCFA<-table }
  
  if(y==3){ # if LCFA
    table_LCFA<-table }
}


# specific plot for SCFA
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Acetic", "Propionic","Butyric","Valeric"), "", "")
signif_y<-tapply(table_SCFA$value, table_SCFA$variable, max )
table_SCFA$signif_y<-as.numeric(signif_y[table_SCFA$variable]) + 10 # plus 5 to put the sign higher

ggplot(data = table_SCFA, mapping = aes(x = variable, y = value, fill = Time)) +
  scale_fill_manual(values = c("T0" = alpha("#FFD700", 0.5), "T1" = alpha("#FFB347", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label = signif, y = signif_y, group = variable), 
            color = "darkgray", size = 8) +
  theme(axis.text.x = element_text(angle = -20, size = 9.5, vjust = 1, hjust = 0.05),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 11),
        panel.grid.minor.y = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.35),
        title = element_text(size = 10.5)) +
  scale_y_log10(breaks = c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50, seq(100, 500, 100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width = 0.65, size = 0.5, outlier.size = 1) + 
  labs(title = "SCFA concentration", x = "", fill = "", y = "Percentage (%)")
ggsave(filename = "FFA/SCFA_zoom_plot.png", height=4.5, width = 7, dpi=300)


# specific plot for MCFA
head(table_MCFA, n=2)
unique(table_MCFA$variable)
table_MCFA$signif<-ifelse(table_MCFA$variable %in% c("Octanoic","Decanoic", "Dodecanoic"), "", "")   
signif_y<-tapply(table_MCFA$value, table_MCFA$variable, max )
table_MCFA$signif_y<-as.numeric(signif_y[table_MCFA$variable]) + 2 # plus 5 to put the sign higher

ggplot(data=table_MCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("#FFD700", 0.5), "T1" = alpha("#FFB347", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 500,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="MCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA/MCFA_zoom_plot.png", height=4.5, width = 7, dpi=300)


# specific plot for LCFA
head(table_LCFA, n=2)
table_LCFA$signif<-ifelse(table_LCFA$variable %in% c("Hexadecanoic","Octadecanoic"), "", "")   
signif_y<-tapply(table_LCFA$value, table_LCFA$variable, max )
table_LCFA$signif_y<-as.numeric(signif_y[table_LCFA$variable]) + 10 # plus 5 to put the sign higher

# View(table_LCFA)
ggplot(data=table_LCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("#FFD700", 0.5), "T1" = alpha("#FFB347", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 700,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="LCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA/LCFA_zoom_plot.png", height=4.5, width = 7, dpi=300)

########### PCoA ########### 
list_PCoA <- list("SCFA" = SCFA_prop, "MCFA" = MCFA_prop, "LCFA" = LCFA_prop)

p_values <- NULL
for (x in 1:length(list_PCoA)) {
  rm(data, coord)
  data <- as.data.frame(list_PCoA[[x]])
  row.names(data) <- data$ID
  
  if (x == 3) {
    data <- data[!is.na(data$Octadecanoic),] 
  }
  
  dist <- vegdist(data[, !colnames(data) %in% c("ID","Time","Sample_Name","Group","Type_diet")], method = "bray") # calcolare dissimilarità Bray-Curtis
  coord <- cmdscale(dist, k = 2, eig = TRUE) # PCoA
  
  # Convertire le coordinate in un dataframe e aggiungere ID e Time
  coord_df <- as.data.frame(coord$points)
  colnames(coord_df) <- c("PC1", "PC2")
  coord_df$Subject <- row.names(coord_df)
  coord_df$Time <- data$Time
  
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Time)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("T0" = alpha("#FFD700", 0.5), "T1" = alpha("#FFB347", 0.5))) +
    stat_ellipse(size = 0.2) +
    geom_text(aes(label = Subject), size = 2, color = "black") +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of Morning type of diet subjects\n with Bray-Curtis dissimilarity")
    )
  
  # again but with no names
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Time)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("T0" = alpha("#FFD700", 0.5), "T1" = alpha("#FFB347", 0.5))) +
    stat_ellipse(size = 0.2) +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of Morning type of diet\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("FFA/PCoA_BRAY_", names(list_PCoA)[x], "morning_noname.png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )                        
  
  # checking for significant dispersion
  suppressWarnings(rm(p_value_diver, p_value_disp))
  p_value_diver<-adonis(dist ~ data$Time, data = data, permutations = 9999)
  p_value_disp<-permutest(betadisper(dist, data$Time), permutations = 9999)
  p_value_disp<-as.data.frame(p_value_disp$tab[1,])
  colnames(p_value_disp)<-colnames(p_value_diver$aov.tab[1,])
  
  new_values<-rbind(p_value_diver$aov.tab[1,],p_value_disp)
  row.names(new_values)<-c(paste("Beta_diversity",names(list_PCoA)[x], sep = "_"),
                           paste("Beta_dispersion",names(list_PCoA)[x], sep = "_"))
  p_values<-rbind(p_values,new_values)
}  

p_values$adjusted_BH_p_values<-p.adjust(p_values$`Pr(>F)`, method = "BH")
p_values
write.csv2(p_values, file="FFA/p_values_BRAY_morning.csv", quote = F,row.names = T)


########### total scfa, mcfa, lcfa #############
FFA <- read_excel("Faecal_FFAs.xlsx", sheet=3) 
colnames(FFA)
FFA$Type_diet
FFA[! colnames(FFA) %in% c("ID","Time","Sample_Name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("ID","Time","Sample_Name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Type_diet<-factor(FFA$Type_diet, levels = c("control","morning","evening"))
FFA<-subset.data.frame(FFA, Type_diet=='morning')
exclude_cols <- c("ID", "Time", "Sample_Name", "Group", "Type_diet")
FFA[!colnames(FFA) %in% exclude_cols] <- 
  apply(FFA[!colnames(FFA) %in% exclude_cols], 2, function(x) as.numeric(x) * 100)

SCFA<-FFA[,1:6] 
MCFA<-FFA[,c(1:5,7)] 
LCFA<-FFA[,c(1:5,8)] 
colnames(LCFA)

SCFA_raw<-SCFA
MCFA_raw<-MCFA
LCFA_raw<- LCFA

#WILCOXON TEST
list<-list("SCFA"=SCFA_raw, "MCFA"=MCFA_raw, "LCFA"=LCFA_raw)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("ID", "Time", "Sample_Name", "Group", "Type_diet")]) {
    test<-wilcox.test(temp[[x]]~temp$Time)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  results
  write.csv2(results, file=paste0("FFA/",names(list)[y],"_raw_abund_Mann_Whit_TOT.csv"), row.names = F)
  
  table <- melt(data = temp)
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Time)) +
    scale_fill_manual(values =  c("T0" = alpha("#FFD700", 0.5), "T1" = alpha("#FFB347", 0.5))) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  #ggsave(filename = paste0("FFA_tot/",names(list)[y],"_raw_abundance_TOT.png"), height=3, width = 6, dpi=300)
  
  if(y==2){ # if MCFA
    table_MCFA<-table }
  
  if(y==3){ # if MCFA
    table_LCFA<-table }
}

# specific plot for SCFA
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Acetic", "Propionic","Butyric","Valeric"), "", "")
signif_y<-tapply(table_SCFA$value, table_SCFA$variable, max )
table_SCFA$signif_y<-as.numeric(signif_y[table_SCFA$variable]) + 10 # plus 5 to put the sign higher

ggplot(data = table_SCFA, mapping = aes(x = variable, y = value, fill = Time)) +
  scale_fill_manual(values = c("T0" = alpha("#FFD700", 0.5), "T1" = alpha("#FFB347", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label = signif, y = signif_y, group = variable), 
            color = "darkgray", size = 8) +
  theme(axis.text.x = element_text(angle = -20, size = 9.5, vjust = 1, hjust = 0.05),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 11),
        panel.grid.minor.y = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.35),
        title = element_text(size = 10.5)) +
  scale_y_log10(breaks = c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50, seq(100, 500, 100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width = 0.65, size = 0.5, outlier.size = 1) + 
  labs(title = "SCFA concentration", x = "", fill = "", y = "Percentage (%)")
ggsave(filename = "FFA/SCFA_zoom_plot_TOT.png", height=4.5, width = 7, dpi=300)


# specific plot for MCFA
head(table_MCFA, n=2)
unique(table_MCFA$variable)
table_MCFA$signif<-ifelse(table_MCFA$variable %in% c("Octanoic","Decanoic", "Dodecanoic"), "", "")   
signif_y<-tapply(table_MCFA$value, table_MCFA$variable, max )
table_MCFA$signif_y<-as.numeric(signif_y[table_MCFA$variable]) + 2 # plus 5 to put the sign higher

ggplot(data=table_MCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("#FFD700", 0.5), "T1" = alpha("#FFB347", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 500,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="MCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA/MCFA_zoom_plot_TOT.png", height=4.5, width = 7, dpi=300)


# specific plot for LCFA
head(table_LCFA, n=2)
table_LCFA$signif<-ifelse(table_LCFA$variable %in% c("Hexadecanoic","Octadecanoic"), "", "")   
signif_y<-tapply(table_LCFA$value, table_LCFA$variable, max )
table_LCFA$signif_y<-as.numeric(signif_y[table_LCFA$variable]) + 10 # plus 5 to put the sign higher

# View(table_LCFA)
ggplot(data=table_LCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("#FFD700", 0.5), "T1" = alpha("#FFB347", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 700,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="LCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA/LCFA_zoom_plot_TOT.png", height=4.5, width = 7, dpi=300)

########## EVENING ##########
options(scipen=100)
#setwd('SARA/Bollous pemphigoid')
FFA <- read_excel("Faecal_FFAs.xlsx", sheet=2) 
colnames(FFA)
FFA$Type_diet
FFA[! colnames(FFA) %in% c("ID","Time","Sample_Name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("ID","Time","Sample_Name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Type_diet<-factor(FFA$Type_diet, levels = c("control","morning","evening"))
FFA<-subset.data.frame(FFA, Type_diet=='evening')
exclude_cols <- c("ID", "Time", "Sample_Name", "Group", "Type_diet")
FFA[!colnames(FFA) %in% exclude_cols] <- 
  apply(FFA[!colnames(FFA) %in% exclude_cols], 2, function(x) as.numeric(x) * 100)

SCFA<-FFA[,1:12] # from acetic to valeric 
MCFA<-FFA[,c(1:5,13:21)] # to dodecanoic
LCFA<-FFA[,c(1:5,22:24)] # to Octadecanoic, before the raw total columns
colnames(LCFA)

######################### NORMALIZATIONS ##############################
SCFA_raw<-as.data.frame(SCFA)
MCFA_raw<-as.data.frame(MCFA)
LCFA_raw<- as.data.frame(LCFA)

# #### normalization through relative abundance in each sample, e.g. https://www.nature.com/articles/s41598-021-81420-3
# #### normalization of FFA respect total SCFA, MCFA or LCFA
SCFA_prop<-SCFA
SCFA_prop[ !colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-SCFA_prop[!colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(SCFA_prop[!colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
SCFA_prop<-as.data.frame(SCFA_prop)
# 
MCFA_prop<-MCFA
MCFA_prop[ !colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-MCFA_prop[!colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(MCFA_prop[!colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
MCFA_prop<-as.data.frame(MCFA_prop)
# 
LCFA_prop<-LCFA
LCFA_prop[ !colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-LCFA_prop[!colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(LCFA_prop[!colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
LCFA_prop<-as.data.frame(LCFA_prop)
# rowSums(LCFA[!colnames(LCFA)%in% c("Subject","Condition")])

####################### MANN WITHNEY TESTs ####################
#dir.create("FFA")

##### FFA raw abundances
list<-list("SCFA"=SCFA_raw, "MCFA"=MCFA_raw, "LCFA"=LCFA_raw)
names(list)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("ID","Time","Sample_Name","Group","Type_diet")]) {
    test<-wilcox.test(temp[[x]]~temp$Time)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  write.csv2(results, file=paste0("FFA/",names(list)[y],"_raw_abund_Mann_Whit.csv"), row.names = F)
  
  table <- melt(data = temp, id.vars = c("ID","Time","Sample_Name","Group","Type_diet"))
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Time)) +
    scale_fill_manual(values =c("T0" = "#4682B4", "T1" = "#0047AB")) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  #ggsave(filename = paste0("FFA/",names(list)[y],"_raw_abundance.png"), height=3, width = 6, dpi=300)
  
  if(y==2){ # if MCFA
    table_MCFA<-table }
  
  if(y==3){ # if LCFA
    table_LCFA<-table }
}


# specific plot for SCFA
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Acetic", "Propionic","Butyric","Valeric"), "", "")
signif_y<-tapply(table_SCFA$value, table_SCFA$variable, max )
table_SCFA$signif_y<-as.numeric(signif_y[table_SCFA$variable]) + 10 # plus 5 to put the sign higher

ggplot(data = table_SCFA, mapping = aes(x = variable, y = value, fill = Time)) +
  scale_fill_manual(values = c("T0" = alpha("#4682B4", 0.5), "T1" = alpha("#0047AB", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label = signif, y = signif_y, group = variable), 
            color = "darkgray", size = 8) +
  theme(axis.text.x = element_text(angle = -20, size = 9.5, vjust = 1, hjust = 0.05),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 11),
        panel.grid.minor.y = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.35),
        title = element_text(size = 10.5)) +
  scale_y_log10(breaks = c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50, seq(100, 500, 100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width = 0.65, size = 0.5, outlier.size = 1) + 
  labs(title = "SCFA concentration", x = "", fill = "", y = "Percentage (%)")
ggsave(filename = "FFA/SCFA_zoom_plot.png", height=4.5, width = 7, dpi=300)


# specific plot for MCFA
head(table_MCFA, n=2)
unique(table_MCFA$variable)
table_MCFA$signif<-ifelse(table_MCFA$variable %in% c("Octanoic","Decanoic", "Dodecanoic"), "", "")   
signif_y<-tapply(table_MCFA$value, table_MCFA$variable, max )
table_MCFA$signif_y<-as.numeric(signif_y[table_MCFA$variable]) + 2 # plus 5 to put the sign higher

ggplot(data=table_MCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values =  c("T0" = alpha("#4682B4", 0.5), "T1" = alpha("#0047AB", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 500,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="MCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA/MCFA_zoom_plot.png", height=4.5, width = 7, dpi=300)

# specific plot for LCFA
head(table_LCFA, n=2)
table_LCFA$signif<-ifelse(table_LCFA$variable %in% c("Hexadecanoic","Octadecanoic"), "", "")   
signif_y<-tapply(table_LCFA$value, table_LCFA$variable, max )
table_LCFA$signif_y<-as.numeric(signif_y[table_LCFA$variable]) + 10 # plus 5 to put the sign higher

# View(table_LCFA)
ggplot(data=table_LCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values =  c("T0" = alpha("#4682B4", 0.5), "T1" = alpha("#0047AB", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 700,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="LCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA/LCFA_zoom_plot.png", height=4.5, width = 7, dpi=300)

########### PCoA ########### 
list_PCoA <- list("SCFA" = SCFA_prop, "MCFA" = MCFA_prop, "LCFA" = LCFA_prop)

p_values <- NULL
for (x in 1:length(list_PCoA)) {
  rm(data, coord)
  data <- as.data.frame(list_PCoA[[x]])
  row.names(data) <- data$ID
  
  if (x == 3) {
    data <- data[!is.na(data$Octadecanoic),] 
  }
  
  dist <- vegdist(data[, !colnames(data) %in% c("ID","Time","Sample_Name","Group","Type_diet")], method = "bray") # calcolare dissimilarità Bray-Curtis
  coord <- cmdscale(dist, k = 2, eig = TRUE) # PCoA
  
  coord_df <- as.data.frame(coord$points)
  colnames(coord_df) <- c("PC1", "PC2")
  coord_df$Subject <- row.names(coord_df)
  coord_df$Time <- data$Time
  
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Time)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values =  c("T0" = alpha("#4682B4", 0.5), "T1" = alpha("#0047AB", 0.5))) +
    stat_ellipse(size = 0.2) +
    geom_text(aes(label = Subject), size = 2, color = "black") +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of Evening type of diet subjects\n with Bray-Curtis dissimilarity")
    )
  
  # again but with no names
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Time)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values =  c("T0" = alpha("#4682B4", 0.5), "T1" = alpha("#0047AB", 0.5))) +
    stat_ellipse(size = 0.2) +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of Evening type of diet\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("FFA/PCoA_BRAY_", names(list_PCoA)[x], "morning_noname.png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )                        
  
  # checking for significant dispersion
  suppressWarnings(rm(p_value_diver, p_value_disp))
  p_value_diver<-adonis(dist ~ data$Time, data = data, permutations = 9999)
  p_value_disp<-permutest(betadisper(dist, data$Time), permutations = 9999)
  p_value_disp<-as.data.frame(p_value_disp$tab[1,])
  colnames(p_value_disp)<-colnames(p_value_diver$aov.tab[1,])
  
  new_values<-rbind(p_value_diver$aov.tab[1,],p_value_disp)
  row.names(new_values)<-c(paste("Beta_diversity",names(list_PCoA)[x], sep = "_"),
                           paste("Beta_dispersion",names(list_PCoA)[x], sep = "_"))
  p_values<-rbind(p_values,new_values)
}  

p_values$adjusted_BH_p_values<-p.adjust(p_values$`Pr(>F)`, method = "BH")
p_values
write.csv2(p_values, file="FFA/p_values_BRAY_evening.csv", quote = F,row.names = T)

########### total scfa, mcfa, lcfa #############
FFA <- read_excel("Faecal_FFAs.xlsx", sheet=3) 
colnames(FFA)
FFA$Type_diet
FFA[! colnames(FFA) %in% c("ID","Time","Sample_Name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("ID","Time","Sample_Name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Type_diet<-factor(FFA$Type_diet, levels = c("control","morning","evening"))
FFA<-subset.data.frame(FFA, Type_diet=='evening')
exclude_cols <- c("ID", "Time", "Sample_Name", "Group", "Type_diet")
FFA[!colnames(FFA) %in% exclude_cols] <- 
  apply(FFA[!colnames(FFA) %in% exclude_cols], 2, function(x) as.numeric(x) * 100)

SCFA<-FFA[,1:6] 
MCFA<-FFA[,c(1:5,7)] 
LCFA<-FFA[,c(1:5,8)] 
colnames(LCFA)

SCFA_raw<-SCFA
MCFA_raw<-MCFA
LCFA_raw<- LCFA

#WILCOXON TEST
list<-list("SCFA"=SCFA_raw, "MCFA"=MCFA_raw, "LCFA"=LCFA_raw)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("ID", "Time", "Sample_Name", "Group", "Type_diet")]) {
    test<-wilcox.test(temp[[x]]~temp$Time)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  results
  write.csv2(results, file=paste0("FFA/",names(list)[y],"_raw_abund_Mann_Whit_TOT.csv"), row.names = F)
  
  table <- melt(data = temp)
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Time)) +
    scale_fill_manual(values = c("T0" = alpha("#4682B4", 0.5), "T1" = alpha("#0047AB", 0.5))) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  #ggsave(filename = paste0("FFA_tot/",names(list)[y],"_raw_abundance_TOT.png"), height=3, width = 6, dpi=300)
  
  if(y==2){ # if MCFA
    table_MCFA<-table }
  
  if(y==3){ # if MCFA
    table_LCFA<-table }
}

# specific plot for SCFA
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Acetic", "Propionic","Butyric","Valeric"), "", "")
signif_y<-tapply(table_SCFA$value, table_SCFA$variable, max )
table_SCFA$signif_y<-as.numeric(signif_y[table_SCFA$variable]) + 10 # plus 5 to put the sign higher

ggplot(data = table_SCFA, mapping = aes(x = variable, y = value, fill = Time)) +
  scale_fill_manual(values =c("T0" = alpha("#4682B4", 0.5), "T1" = alpha("#0047AB", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label = signif, y = signif_y, group = variable), 
            color = "darkgray", size = 8) +
  theme(axis.text.x = element_text(angle = -20, size = 9.5, vjust = 1, hjust = 0.05),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 11),
        panel.grid.minor.y = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.35),
        title = element_text(size = 10.5)) +
  scale_y_log10(breaks = c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50, seq(100, 500, 100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width = 0.65, size = 0.5, outlier.size = 1) + 
  labs(title = "SCFA concentration", x = "", fill = "", y = "Percentage (%)")
ggsave(filename = "FFA/SCFA_zoom_plot_TOT.png", height=4.5, width = 7, dpi=300)


# specific plot for MCFA
head(table_MCFA, n=2)
unique(table_MCFA$variable)
table_MCFA$signif<-ifelse(table_MCFA$variable %in% c("Octanoic","Decanoic", "Dodecanoic"), "", "")  
signif_y<-tapply(table_MCFA$value, table_MCFA$variable, max )
table_MCFA$signif_y<-as.numeric(signif_y[table_MCFA$variable]) + 2 # plus 5 to put the sign higher

ggplot(data=table_MCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values =c("T0" = alpha("#4682B4", 0.5), "T1" = alpha("#0047AB", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 500,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="MCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA/MCFA_zoom_plot_TOT.png", height=4.5, width = 7, dpi=300)


# specific plot for LCFA
head(table_LCFA, n=2)
table_LCFA$signif<-ifelse(table_LCFA$variable %in% c("Hexadecanoic","Octadecanoic"), "", "")   
signif_y<-tapply(table_LCFA$value, table_LCFA$variable, max )
table_LCFA$signif_y<-as.numeric(signif_y[table_LCFA$variable]) + 10 # plus 5 to put the sign higher

# View(table_LCFA)
ggplot(data=table_LCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("#4682B4", 0.5), "T1" = alpha("#0047AB", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 700,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="LCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA/LCFA_zoom_plot_TOT.png", height=4.5, width = 7, dpi=300)

########## CONTROL ##########
options(scipen=100)
#setwd('SARA/Bollous pemphigoid')
FFA <- read_excel("Faecal_FFAs.xlsx", sheet=2) 
colnames(FFA)
FFA$Type_diet
FFA[! colnames(FFA) %in% c("ID","Time","Sample_Name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("ID","Time","Sample_Name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Type_diet<-factor(FFA$Type_diet, levels = c("control","morning","evening"))
FFA<-subset.data.frame(FFA, Type_diet=='control')
exclude_cols <- c("ID", "Time", "Sample_Name", "Group", "Type_diet")
FFA[!colnames(FFA) %in% exclude_cols] <- 
  apply(FFA[!colnames(FFA) %in% exclude_cols], 2, function(x) as.numeric(x) * 100)

SCFA<-FFA[,1:12] # from acetic to valeric 
MCFA<-FFA[,c(1:5,13:21)] # to dodecanoic
LCFA<-FFA[,c(1:5,22:24)] # to Octadecanoic, before the raw total columns
colnames(LCFA)

######################### NORMALIZATIONS ##############################
SCFA_raw<-as.data.frame(SCFA)
MCFA_raw<-as.data.frame(MCFA)
LCFA_raw<- as.data.frame(LCFA)

# #### normalization through relative abundance in each sample, e.g. https://www.nature.com/articles/s41598-021-81420-3
# #### normalization of FFA respect total SCFA, MCFA or LCFA
SCFA_prop<-SCFA
SCFA_prop[ !colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-SCFA_prop[!colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(SCFA_prop[!colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
SCFA_prop<-as.data.frame(SCFA_prop)
# 
MCFA_prop<-MCFA
MCFA_prop[ !colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-MCFA_prop[!colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(MCFA_prop[!colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
MCFA_prop<-as.data.frame(MCFA_prop)
# 
LCFA_prop<-LCFA
LCFA_prop[ !colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-LCFA_prop[!colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(LCFA_prop[!colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
LCFA_prop<-as.data.frame(LCFA_prop)
# rowSums(LCFA[!colnames(LCFA)%in% c("Subject","Condition")])

####################### MANN WITHNEY TESTs ####################
#dir.create("FFA")

##### FFA raw abundances
list<-list("SCFA"=SCFA_raw, "MCFA"=MCFA_raw, "LCFA"=LCFA_raw)
names(list)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("ID","Time","Sample_Name","Group","Type_diet")]) {
    test<-wilcox.test(temp[[x]]~temp$Time)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  write.csv2(results, file=paste0("FFA/",names(list)[y],"_raw_abund_Mann_Whit.csv"), row.names = F)
  
  table <- melt(data = temp, id.vars = c("ID","Time","Sample_Name","Group","Type_diet"))
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Time)) +
    scale_fill_manual(values =c("T0" = "#CCCCCC", "T1" = "#555555")) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  #ggsave(filename = paste0("FFA/",names(list)[y],"_raw_abundance.png"), height=3, width = 6, dpi=300)
  
  if(y==2){ # if MCFA
    table_MCFA<-table }
  
  if(y==3){ # if LCFA
    table_LCFA<-table }
}


# specific plot for SCFA
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Acetic", "Propionic","Butyric","Valeric"), "", "")
signif_y<-tapply(table_SCFA$value, table_SCFA$variable, max )
table_SCFA$signif_y<-as.numeric(signif_y[table_SCFA$variable]) + 10 # plus 5 to put the sign higher

ggplot(data = table_SCFA, mapping = aes(x = variable, y = value, fill = Time)) +
  scale_fill_manual(values = c("T0" = alpha("#CCCCCC", 0.5), "T1" = alpha("#555555", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label = signif, y = signif_y, group = variable), 
            color = "darkgray", size = 8) +
  theme(axis.text.x = element_text(angle = -20, size = 9.5, vjust = 1, hjust = 0.05),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 11),
        panel.grid.minor.y = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.35),
        title = element_text(size = 10.5)) +
  scale_y_log10(breaks = c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50, seq(100, 500, 100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width = 0.65, size = 0.5, outlier.size = 1) + 
  labs(title = "SCFA concentration", x = "", fill = "", y = "Percentage (%)")
ggsave(filename = "FFA/SCFA_zoom_plot.png", height=4.5, width = 7, dpi=300)


# specific plot for MCFA
head(table_MCFA, n=2)
unique(table_MCFA$variable)
table_MCFA$signif<-ifelse(table_MCFA$variable %in% c("Octanoic","Decanoic", "Dodecanoic"), "", "")   
signif_y<-tapply(table_MCFA$value, table_MCFA$variable, max )
table_MCFA$signif_y<-as.numeric(signif_y[table_MCFA$variable]) + 2 # plus 5 to put the sign higher

ggplot(data=table_MCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values =  c("T0" = alpha("#CCCCCC", 0.5), "T1" = alpha("#555555", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 500,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="MCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA/MCFA_zoom_plot.png", height=4.5, width = 7, dpi=300)

# specific plot for LCFA
head(table_LCFA, n=2)
table_LCFA$signif<-ifelse(table_LCFA$variable %in% c("Hexadecanoic","Octadecanoic"), "", "")   
signif_y<-tapply(table_LCFA$value, table_LCFA$variable, max )
table_LCFA$signif_y<-as.numeric(signif_y[table_LCFA$variable]) + 10 # plus 5 to put the sign higher

# View(table_LCFA)
ggplot(data=table_LCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("#CCCCCC", 0.5), "T1" = alpha("#555555", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 700,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="LCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA/LCFA_zoom_plot.png", height=4.5, width = 7, dpi=300)

########### PCoA ########### 
list_PCoA <- list("SCFA" = SCFA_prop, "MCFA" = MCFA_prop, "LCFA" = LCFA_prop)

p_values <- NULL
for (x in 1:length(list_PCoA)) {
  rm(data, coord)
  data <- as.data.frame(list_PCoA[[x]])
  row.names(data) <- data$ID
  
  if (x == 3) {
    data <- data[!is.na(data$Octadecanoic),] # rimuovere NA in LCFA
  }
  
  dist <- vegdist(data[, !colnames(data) %in% c("ID","Time","Sample_Name","Group","Type_diet")], method = "bray") # calcolare dissimilarità Bray-Curtis
  coord <- cmdscale(dist, k = 2, eig = TRUE) # PCoA
  
  # Convertire le coordinate in un dataframe e aggiungere ID e Time
  coord_df <- as.data.frame(coord$points)
  colnames(coord_df) <- c("PC1", "PC2")
  coord_df$Subject <- row.names(coord_df)
  coord_df$Time <- data$Time
  
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Time)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values =  c("T0" = alpha("#CCCCCC", 0.5), "T1" = alpha("#555555", 0.5))) +
    stat_ellipse(size = 0.2) +
    geom_text(aes(label = Subject), size = 2, color = "black") +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of Control diet subjects\n with Bray-Curtis dissimilarity")
    )
  
  # again but with no names
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Time)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values =  c("T0" = alpha("#CCCCCC", 0.5), "T1" = alpha("#555555", 0.5))) +
    stat_ellipse(size = 0.2) +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of Control diet\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("FFA/PCoA_BRAY_", names(list_PCoA)[x], "morning_noname.png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )                        
  
  # checking for significant dispersion
  suppressWarnings(rm(p_value_diver, p_value_disp))
  p_value_diver<-adonis(dist ~ data$Time, data = data, permutations = 9999)
  p_value_disp<-permutest(betadisper(dist, data$Time), permutations = 9999)
  p_value_disp<-as.data.frame(p_value_disp$tab[1,])
  colnames(p_value_disp)<-colnames(p_value_diver$aov.tab[1,])
  
  new_values<-rbind(p_value_diver$aov.tab[1,],p_value_disp)
  row.names(new_values)<-c(paste("Beta_diversity",names(list_PCoA)[x], sep = "_"),
                           paste("Beta_dispersion",names(list_PCoA)[x], sep = "_"))
  p_values<-rbind(p_values,new_values)
}  

p_values$adjusted_BH_p_values<-p.adjust(p_values$`Pr(>F)`, method = "BH")
p_values
write.csv2(p_values, file="FFA/p_values_BRAY_control.csv", quote = F,row.names = T)

########### total scfa, mcfa, lcfa #############
FFA <- read_excel("Faecal_FFAs.xlsx", sheet=3) 
colnames(FFA)
FFA$Type_diet
FFA[! colnames(FFA) %in% c("ID","Time","Sample_Name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("ID","Time","Sample_Name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Type_diet<-factor(FFA$Type_diet, levels = c("control","morning","evening"))
FFA<-subset.data.frame(FFA, Type_diet=='control')
exclude_cols <- c("ID", "Time", "Sample_Name", "Group", "Type_diet")
FFA[!colnames(FFA) %in% exclude_cols] <- 
  apply(FFA[!colnames(FFA) %in% exclude_cols], 2, function(x) as.numeric(x) * 100)

SCFA<-FFA[,1:6] 
MCFA<-FFA[,c(1:5,7)] 
LCFA<-FFA[,c(1:5,8)] 
colnames(LCFA)

SCFA_raw<-SCFA
MCFA_raw<-MCFA
LCFA_raw<- LCFA

#WILCOXON TEST
list<-list("SCFA"=SCFA_raw, "MCFA"=MCFA_raw, "LCFA"=LCFA_raw)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("ID", "Time", "Sample_Name", "Group", "Type_diet")]) {
    test<-wilcox.test(temp[[x]]~temp$Time)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  results
  write.csv2(results, file=paste0("FFA/",names(list)[y],"_raw_abund_Mann_Whit_TOT.csv"), row.names = F)
  
  table <- melt(data = temp)
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Time)) +
    scale_fill_manual(values = c("T0" = alpha("#CCCCCC", 0.5), "T1" = alpha("#555555", 0.5))) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  #ggsave(filename = paste0("FFA_tot/",names(list)[y],"_raw_abundance_TOT.png"), height=3, width = 6, dpi=300)
  
  if(y==2){ # if MCFA
    table_MCFA<-table }
  
  if(y==3){ # if MCFA
    table_LCFA<-table }
}

# specific plot for SCFA
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Acetic", "Propionic","Butyric","Valeric"), "", "")
signif_y<-tapply(table_SCFA$value, table_SCFA$variable, max )
table_SCFA$signif_y<-as.numeric(signif_y[table_SCFA$variable]) + 10 # plus 5 to put the sign higher

ggplot(data = table_SCFA, mapping = aes(x = variable, y = value, fill = Time)) +
  scale_fill_manual(values =c("T0" = alpha("#CCCCCC", 0.5), "T1" = alpha("#555555", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label = signif, y = signif_y, group = variable), 
            color = "darkgray", size = 8) +
  theme(axis.text.x = element_text(angle = -20, size = 9.5, vjust = 1, hjust = 0.05),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 11),
        panel.grid.minor.y = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.35),
        title = element_text(size = 10.5)) +
  scale_y_log10(breaks = c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50, seq(100, 500, 100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width = 0.65, size = 0.5, outlier.size = 1) + 
  labs(title = "SCFA concentration", x = "", fill = "", y = "Percentage (%)")
ggsave(filename = "FFA/SCFA_zoom_plot_TOT.png", height=4.5, width = 7, dpi=300)


# specific plot for MCFA
head(table_MCFA, n=2)
unique(table_MCFA$variable)
table_MCFA$signif<-ifelse(table_MCFA$variable %in% c("Octanoic","Decanoic", "Dodecanoic"), "", "")  
signif_y<-tapply(table_MCFA$value, table_MCFA$variable, max )
table_MCFA$signif_y<-as.numeric(signif_y[table_MCFA$variable]) + 2 # plus 5 to put the sign higher

ggplot(data=table_MCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values =c("T0" = alpha("#CCCCCC", 0.5), "T1" = alpha("#555555", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 500,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="MCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA/MCFA_zoom_plot_TOT.png", height=4.5, width = 7, dpi=300)


# specific plot for LCFA
head(table_LCFA, n=2)
table_LCFA$signif<-ifelse(table_LCFA$variable %in% c("Hexadecanoic","Octadecanoic"), "", "")  
signif_y<-tapply(table_LCFA$value, table_LCFA$variable, max )
table_LCFA$signif_y<-as.numeric(signif_y[table_LCFA$variable]) + 10 # plus 5 to put the sign higher

# View(table_LCFA)
ggplot(data=table_LCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("#CCCCCC", 0.5), "T1" = alpha("#555555", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 700,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="LCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA/LCFA_zoom_plot_TOT.png", height=4.5, width = 7, dpi=300)

########## INTERVENTION ##########
options(scipen=100)
#setwd('SARA/Bollous pemphigoid')
FFA <- read_excel("Faecal_FFAs.xlsx", sheet=2) 
colnames(FFA)
FFA$Type_diet
FFA[! colnames(FFA) %in% c("ID","Time","Sample_Name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("ID","Time","Sample_Name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Type_diet<-factor(FFA$Type_diet, levels = c("control","morning","evening"))
FFA<-subset.data.frame(FFA, Group=='intervention')
exclude_cols <- c("ID", "Time", "Sample_Name", "Group", "Type_diet")
FFA[!colnames(FFA) %in% exclude_cols] <- 
  apply(FFA[!colnames(FFA) %in% exclude_cols], 2, function(x) as.numeric(x) * 100)

SCFA<-FFA[,1:12] # from acetic to valeric 
MCFA<-FFA[,c(1:5,13:17)] # to dodecanoic
LCFA<-FFA[,c(1:5,18:20)] # to Octadecanoic, before the raw total columns
colnames(LCFA)

######################### NORMALIZATIONS ##############################
SCFA_raw<-as.data.frame(SCFA)
MCFA_raw<-as.data.frame(MCFA)
LCFA_raw<- as.data.frame(LCFA)

# #### normalization through relative abundance in each sample, e.g. https://www.nature.com/articles/s41598-021-81420-3
# #### normalization of FFA respect total SCFA, MCFA or LCFA
SCFA_prop<-SCFA
SCFA_prop[ !colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-SCFA_prop[!colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(SCFA_prop[!colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
SCFA_prop<-as.data.frame(SCFA_prop)
# 
MCFA_prop<-MCFA
MCFA_prop[ !colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-MCFA_prop[!colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(MCFA_prop[!colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
MCFA_prop<-as.data.frame(MCFA_prop)
# 
LCFA_prop<-LCFA
LCFA_prop[ !colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-LCFA_prop[!colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(LCFA_prop[!colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
LCFA_prop<-as.data.frame(LCFA_prop)
# rowSums(LCFA[!colnames(LCFA)%in% c("Subject","Condition")])

####################### MANN WITHNEY TESTs ####################
#dir.create("FFA")

##### FFA raw abundances
list<-list("SCFA"=SCFA_raw, "MCFA"=MCFA_raw, "LCFA"=LCFA_raw)
names(list)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("ID","Time","Sample_Name","Group","Type_diet")]) {
    test<-wilcox.test(temp[[x]]~temp$Time)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  write.csv2(results, file=paste0(names(list)[y],"_raw_abund_Mann_Whit.csv"), row.names = F)
  
  table <- melt(data = temp, id.vars = c("ID","Time","Sample_Name","Group","Type_diet"))
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Time)) +
    scale_fill_manual(values =c("T1" = "red3", "T0" = "#FF6666"))  +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  #ggsave(filename = paste0("FFA/",names(list)[y],"_raw_abundance.png"), height=3, width = 6, dpi=300)
  
  if(y==2){ # if MCFA
    table_MCFA<-table }
  
  if(y==3){ # if LCFA
    table_LCFA<-table }
}


# specific plot for SCFA
head(table_SCFA, n=2)
unique(table_SCFA$variable)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("isoValeric", "2MetilButiric"), "", "")
signif_y<-tapply(table_SCFA$value, table_SCFA$variable, max )
table_SCFA$signif_y<-as.numeric(signif_y[table_SCFA$variable]) + 10 # plus 5 to put the sign higher

ggplot(data = table_SCFA, mapping = aes(x = variable, y = value, fill = Time)) +
  scale_fill_manual(values = c("T0" = alpha("#FF6666", 0.5), "T1" = alpha("red3", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label = signif, y = signif_y, group = variable), 
            color = "darkgray", size = 8) +
  theme(axis.text.x = element_text(angle = -20, size = 9.5, vjust = 1, hjust = 0.05),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 11),
        panel.grid.minor.y = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.35),
        title = element_text(size = 10.5)) +
  scale_y_log10(breaks = c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50, seq(100, 500, 100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width = 0.65, size = 0.5, outlier.size = 1) + 
  labs(title = "SCFA concentration", x = "", fill = "", y = "Percentage (%)")
ggsave(filename = "FFA/SCFA_zoom_plot.png", height=4.5, width = 7, dpi=300)


# specific plot for MCFA
head(table_MCFA, n=2)
unique(table_MCFA$variable)
table_MCFA$signif<-ifelse(table_MCFA$variable %in% c("Octanoic","Nonanoic"), "*", "")   
signif_y<-tapply(table_MCFA$value, table_MCFA$variable, max )
table_MCFA$signif_y<-as.numeric(signif_y[table_MCFA$variable]) + 2 # plus 5 to put the sign higher

ggplot(data=table_MCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("#FF6666", 0.5), "T1" = alpha("red3", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 500,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="MCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA/MCFA_zoom_plot.png", height=4.5, width = 7, dpi=300)

# specific plot for LCFA
head(table_LCFA, n=2)
table_LCFA$signif<-ifelse(table_LCFA$variable %in% c("Hexadecanoic","Octadecanoic"), "", "")   
signif_y<-tapply(table_LCFA$value, table_LCFA$variable, max )
table_LCFA$signif_y<-as.numeric(signif_y[table_LCFA$variable]) + 10 # plus 5 to put the sign higher

# View(table_LCFA)
ggplot(data=table_LCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("#FF6666", 0.5), "T1" = alpha("red3", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 700,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="LCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA/LCFA_zoom_plot.png", height=4.5, width = 7, dpi=300)

########### PCoA ########### 
list_PCoA <- list("SCFA" = SCFA_prop, "MCFA" = MCFA_prop, "LCFA" = LCFA_prop)

p_values <- NULL
for (x in 1:length(list_PCoA)) {
  rm(data, coord)
  data <- as.data.frame(list_PCoA[[x]])
  row.names(data) <- data$ID
  
  dist <- vegdist(data[, !colnames(data) %in% c("ID","Time","Sample_Name","Group","Type_diet")], method = "bray") 
  coord <- cmdscale(dist, k = 2, eig = TRUE) # PCoA
  coord_df <- as.data.frame(coord$points)
  colnames(coord_df) <- c("PC1", "PC2")
  coord_df$ID <- row.names(coord_df)
  coord_df$Time <- data$Time
  
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Time)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("T0" = alpha("#FF6666", 0.5), "T1" = alpha("red3", 0.5))) +
    stat_ellipse(size = 0.2) +
    geom_text(aes(label = ID), size = 2, color = "black") +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of Intervention diet subjects\n with Bray-Curtis dissimilarity")
    )
  
  # again but with no names
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Time)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("T0" = alpha("#FF6666", 0.5), "T1" = alpha("red3", 0.5))) +
    stat_ellipse(size = 0.2) +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of Intervention diet\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("PCoA_BRAY_", names(list_PCoA)[x], "Intervention_noname.png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )                        
  
  # checking for significant dispersion
  suppressWarnings(rm(p_value_diver, p_value_disp))
  p_value_diver<-adonis(dist ~ data$Time, data = data, permutations = 9999)
  p_value_disp<-permutest(betadisper(dist, data$Time), permutations = 9999)
  p_value_disp<-as.data.frame(p_value_disp$tab[1,])
  colnames(p_value_disp)<-colnames(p_value_diver$aov.tab[1,])
  
  new_values<-rbind(p_value_diver$aov.tab[1,],p_value_disp)
  row.names(new_values)<-c(paste("Beta_diversity",names(list_PCoA)[x], sep = "_"),
                           paste("Beta_dispersion",names(list_PCoA)[x], sep = "_"))
  p_values<-rbind(p_values,new_values)
}  

p_values$adjusted_BH_p_values<-p.adjust(p_values$`Pr(>F)`, method = "BH")
p_values
write.csv2(p_values, file="FFA/p_values_BRAY_intervention.csv", quote = F,row.names = T)


########### total scfa, mcfa, lcfa #############
FFA <- read_excel("Faecal_FFAs.xlsx", sheet=3) 
colnames(FFA)
FFA$Type_diet
FFA[! colnames(FFA) %in% c("ID","Time","Sample_Name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("ID","Time","Sample_Name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Type_diet<-factor(FFA$Type_diet, levels = c("control","morning","evening"))
FFA<-subset.data.frame(FFA, Group=='intervention')
exclude_cols <- c("ID", "Time", "Sample_Name", "Group", "Type_diet")
FFA[!colnames(FFA) %in% exclude_cols] <- 
  apply(FFA[!colnames(FFA) %in% exclude_cols], 2, function(x) as.numeric(x) * 100)

SCFA<-FFA[,1:6] 
MCFA<-FFA[,c(1:5,7)] 
LCFA<-FFA[,c(1:5,8)] 
colnames(LCFA)

SCFA_raw<-SCFA
MCFA_raw<-MCFA
LCFA_raw<- LCFA

#WILCOXON TEST
list<-list("SCFA"=SCFA_raw, "MCFA"=MCFA_raw, "LCFA"=LCFA_raw)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("ID", "Time", "Sample_Name", "Group", "Type_diet")]) {
    test<-wilcox.test(temp[[x]]~temp$Time)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  results
  #write.csv2(results, file=paste0("FFA/",names(list)[y],"_raw_abund_Mann_Whit_TOT.csv"), row.names = F)
  
  table <- melt(data = temp)
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Time)) +
    scale_fill_manual(values =  c("T0" = alpha("#FF6666", 0.5), "T1" = alpha("red3", 0.5))) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  #ggsave(filename = paste0("FFA_tot/",names(list)[y],"_raw_abundance_TOT.png"), height=3, width = 6, dpi=300)
  
  if(y==2){ # if MCFA
    table_MCFA<-table }
  
  if(y==3){ # if MCFA
    table_LCFA<-table }
}

# specific plot for SCFA
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Acetic", "Propionic","Butyric","Valeric"), "", "")
signif_y<-tapply(table_SCFA$value, table_SCFA$variable, max )
table_SCFA$signif_y<-as.numeric(signif_y[table_SCFA$variable]) + 10 # plus 5 to put the sign higher

ggplot(data = table_SCFA, mapping = aes(x = variable, y = value, fill = Time)) +
  scale_fill_manual(values = c("T0" = alpha("#FF6666", 0.5), "T1" = alpha("red3", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label = signif, y = signif_y, group = variable), 
            color = "darkgray", size = 8) +
  theme(axis.text.x = element_text(angle = -20, size = 9.5, vjust = 1, hjust = 0.05),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 11),
        panel.grid.minor.y = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.35),
        title = element_text(size = 10.5)) +
  scale_y_log10(breaks = c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50, seq(100, 500, 100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width = 0.65, size = 0.5, outlier.size = 1) + 
  labs(title = "SCFA concentration", x = "", fill = "", y = "Percentage (%)")
ggsave(filename = "SCFA_zoom_plot_TOT.png", height=4.5, width = 7, dpi=300)


# specific plot for MCFA
head(table_MCFA, n=2)
unique(table_MCFA$variable)
table_MCFA$signif<-ifelse(table_MCFA$variable %in% c("Octanoic","Decanoic", "Dodecanoic"), "", "")   
signif_y<-tapply(table_MCFA$value, table_MCFA$variable, max )
table_MCFA$signif_y<-as.numeric(signif_y[table_MCFA$variable]) + 2 # plus 5 to put the sign higher

ggplot(data=table_MCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("#FF6666", 0.5), "T1" = alpha("red3", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 500,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="MCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "MCFA_zoom_plot_TOT.png", height=4.5, width = 7, dpi=300)


# specific plot for LCFA
head(table_LCFA, n=2)
table_LCFA$signif<-ifelse(table_LCFA$variable %in% c("Hexadecanoic","Octadecanoic"), "", "")   
signif_y<-tapply(table_LCFA$value, table_LCFA$variable, max )
table_LCFA$signif_y<-as.numeric(signif_y[table_LCFA$variable]) + 10 # plus 5 to put the sign higher

# View(table_LCFA)
ggplot(data=table_LCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("#FF6666", 0.5), "T1" = alpha("red3", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 700,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="LCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "LCFA_zoom_plot_TOT.png", height=4.5, width = 7, dpi=300)


########## CONTROL VS INTERVENTION ##########
options(scipen=100)
#setwd('SARA/Bollous pemphigoid')
FFA <- read_excel("Faecal_FFAs.xlsx", sheet=2) 
colnames(FFA)
FFA$Group
FFA[! colnames(FFA) %in% c("ID","Time","Sample_Name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("ID","Time","Sample_Name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Group<-factor(FFA$Group, levels = c("intervention","control"))
FFA<-subset.data.frame(FFA, Time=='T1')

exclude_cols <- c("ID", "Time", "Sample_Name", "Group", "Type_diet")
FFA[!colnames(FFA) %in% exclude_cols] <- 
  apply(FFA[!colnames(FFA) %in% exclude_cols], 2, function(x) as.numeric(x) * 100)

SCFA<-FFA[,1:12] # from acetic to valeric 
MCFA<-FFA[,c(1:5,13:21)] # to dodecanoic
LCFA<-FFA[,c(1:5,22:24)] # to Octadecanoic, before the raw total columns
colnames(LCFA)

######################### NORMALIZATIONS ##############################
SCFA_raw<-as.data.frame(SCFA)
MCFA_raw<-as.data.frame(MCFA)
LCFA_raw<- as.data.frame(LCFA)

# #### normalization through relative abundance in each sample, e.g. https://www.nature.com/articles/s41598-021-81420-3
# #### normalization of FFA respect total SCFA, MCFA or LCFA
SCFA_prop<-SCFA
SCFA_prop[ !colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-SCFA_prop[!colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(SCFA_prop[!colnames(SCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
SCFA_prop<-as.data.frame(SCFA_prop)
# 
MCFA_prop<-MCFA
MCFA_prop[ !colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-MCFA_prop[!colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(MCFA_prop[!colnames(MCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
MCFA_prop<-as.data.frame(MCFA_prop)
# 
LCFA_prop<-LCFA
LCFA_prop[ !colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]<-LCFA_prop[!colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")]/rowSums(LCFA_prop[!colnames(LCFA_prop)%in% c("ID","Time","Sample_Name","Group","Type_diet")])
LCFA_prop<-as.data.frame(LCFA_prop)
# rowSums(LCFA[!colnames(LCFA)%in% c("Subject","Condition")])

####################### MANN WITHNEY TESTs ####################
#dir.create("FFA")

##### FFA raw abundances
list<-list("SCFA"=SCFA_raw, "MCFA"=MCFA_raw, "LCFA"=LCFA_raw)
names(list)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("ID","Time","Sample_Name","Group","Type_diet")]) {
    test<-wilcox.test(temp[[x]]~temp$Group)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  write.csv2(results, file=paste0(names(list)[y],"_raw_abund_Mann_Whit.csv"), row.names = F)
  
  table <- melt(data = temp, id.vars = c("ID","Time","Sample_Name","Group","Type_diet"))
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Group)) +
    scale_fill_manual(values = c("intervention" = "red3", "control" = '#555555')) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  #ggsave(filename = paste0("FFA/",names(list)[y],"_raw_abundance.png"), height=3, width = 6, dpi=300)
  
  if(y==2){ # if MCFA
    table_MCFA<-table }
  
  if(y==3){ # if LCFA
    table_LCFA<-table }
}


# specific plot for SCFA
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Acetic", "Propionic","Butyric","Valeric"), "", "")
signif_y<-tapply(table_SCFA$value, table_SCFA$variable, max )
table_SCFA$signif_y<-as.numeric(signif_y[table_SCFA$variable]) + 10 # plus 5 to put the sign higher

ggplot(data = table_SCFA, mapping = aes(x = variable, y = value, fill = Group)) +
  scale_fill_manual(values = c("intervention" = alpha("red3", 0.5), "control" = alpha("#555555", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label = signif, y = signif_y, group = variable), 
            color = "darkgray", size = 8) +
  theme(axis.text.x = element_text(angle = -20, size = 9.5, vjust = 1, hjust = 0.05),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 11),
        panel.grid.minor.y = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.35),
        title = element_text(size = 10.5)) +
  scale_y_log10(breaks = c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50, seq(100, 500, 100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width = 0.65, size = 0.5, outlier.size = 1) + 
  labs(title = "SCFA concentration", x = "", fill = "", y = "Percentage (%)")
ggsave(filename = "SCFA_zoom_plot.png", height=4.5, width = 7, dpi=300)


# specific plot for MCFA
head(table_MCFA, n=2)
unique(table_MCFA$variable)
table_MCFA$signif<-ifelse(table_MCFA$variable %in% c("Octanoic","Decanoic", "Dodecanoic"), "", "")   
signif_y<-tapply(table_MCFA$value, table_MCFA$variable, max )
table_MCFA$signif_y<-as.numeric(signif_y[table_MCFA$variable]) + 2 # plus 5 to put the sign higher

ggplot(data=table_MCFA, mapping=aes(x=variable, y=value, fill=Group)) +
  scale_fill_manual(values = c("intervention" = alpha("red3", 0.5), "control" = alpha("#555555", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 500,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="MCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "MCFA_zoom_plot.png", height=4.5, width = 7, dpi=300)


# specific plot for LCFA
head(table_LCFA, n=2)
table_LCFA$signif<-ifelse(table_LCFA$variable %in% c("Hexadecanoic","Octadecanoic"), "", "")   
signif_y<-tapply(table_LCFA$value, table_LCFA$variable, max )
table_LCFA$signif_y<-as.numeric(signif_y[table_LCFA$variable]) + 10 # plus 5 to put the sign higher

# View(table_LCFA)
ggplot(data=table_LCFA, mapping=aes(x=variable, y=value, fill=Group)) +
  scale_fill_manual(values =  c("intervention" = alpha("red3", 0.5), "control" = alpha("#555555", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 700,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="LCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "LCFA_zoom_plot.png", height=4.5, width = 7, dpi=300)

########### PCoA ########### 
list_PCoA <- list("SCFA" = SCFA_prop, "MCFA" = MCFA_prop, "LCFA" = LCFA_prop)

p_values <- NULL
for (x in 1:length(list_PCoA)) {
  rm(data, coord)
  data <- as.data.frame(list_PCoA[[x]])
  row.names(data) <- data$ID
  
  if (x == 3) {
    data <- data[!is.na(data$Octadecanoic),] # rimuovere NA in LCFA
  }
  
  dist <- vegdist(data[, !colnames(data) %in% c("ID","Time","Sample_Name","Group","Type_diet")], method = "bray") 
  coord <- cmdscale(dist, k = 2, eig = TRUE) # PCoA
  
  coord_df <- as.data.frame(coord$points)
  colnames(coord_df) <- c("PC1", "PC2")
  coord_df$Subject <- row.names(coord_df)
  coord_df$Group <- data$Group
  
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Group)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values =  c("intervention" = alpha("red3", 0.5), "control" = alpha("#555555", 0.5))) +
    stat_ellipse(size = 0.2) +
    geom_text(aes(label = Subject), size = 2, color = "black") +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of Morning type of diet subjects\n with Bray-Curtis dissimilarity")
    )
  
  # again but with no names
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Group)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("intervention" = alpha("red3", 0.5), "control" = alpha("#555555", 0.5))) +
    stat_ellipse(size = 0.2) +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of Morning type of diet\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("PCoA_BRAY_", names(list_PCoA)[x], "morning_noname.png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )                        
  
  # checking for significant dispersion
  suppressWarnings(rm(p_value_diver, p_value_disp))
  p_value_diver<-adonis(dist ~ data$Group, data = data, permutations = 9999)
  p_value_disp<-permutest(betadisper(dist, data$Group), permutations = 9999)
  p_value_disp<-as.data.frame(p_value_disp$tab[1,])
  colnames(p_value_disp)<-colnames(p_value_diver$aov.tab[1,])
  
  new_values<-rbind(p_value_diver$aov.tab[1,],p_value_disp)
  row.names(new_values)<-c(paste("Beta_diversity",names(list_PCoA)[x], sep = "_"),
                           paste("Beta_dispersion",names(list_PCoA)[x], sep = "_"))
  p_values<-rbind(p_values,new_values)
}  

p_values$adjusted_BH_p_values<-p.adjust(p_values$`Pr(>F)`, method = "BH")
p_values
write.csv2(p_values, file="FFA/p_values_BRAY.csv", quote = F,row.names = T)


########### total scfa, mcfa, lcfa #############
FFA <- read_excel("Faecal_FFAs.xlsx", sheet=3) 
colnames(FFA)
FFA$Group
FFA[! colnames(FFA) %in% c("ID","Time","Sample_Name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("ID","Time","Sample_Name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Group<-factor(FFA$Group, levels = c("control","intervention"))
FFA<-subset.data.frame(FFA, Time=='T1')

exclude_cols <- c("ID", "Time", "Sample_Name", "Group", "Type_diet")
FFA[!colnames(FFA) %in% exclude_cols] <- 
  apply(FFA[!colnames(FFA) %in% exclude_cols], 2, function(x) as.numeric(x) * 100)

SCFA<-FFA[,1:6] 
MCFA<-FFA[,c(1:5,7)] 
LCFA<-FFA[,c(1:5,8)] 
colnames(LCFA)

SCFA_raw<-SCFA
MCFA_raw<-MCFA
LCFA_raw<- LCFA

#WILCOXON TEST
list<-list("SCFA"=SCFA_raw, "MCFA"=MCFA_raw, "LCFA"=LCFA_raw)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("ID", "Time", "Sample_Name", "Group", "Type_diet")]) {
    test<-wilcox.test(temp[[x]]~temp$Group)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  results
  write.csv2(results, file=paste0(names(list)[y],"_raw_abund_Mann_Whit_TOT.csv"), row.names = F)
  
  table <- melt(data = temp)
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Group)) +
    scale_fill_manual(values =   c("intervention" = alpha("red3", 0.5), "control" = alpha("#555555", 0.5))) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  #ggsave(filename = paste0("FFA_tot/",names(list)[y],"_raw_abundance_TOT.png"), height=3, width = 6, dpi=300)
  
  if(y==2){ # if MCFA
    table_MCFA<-table }
  
  if(y==3){ # if MCFA
    table_LCFA<-table }
}

# specific plot for SCFA
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Total_SCFA"), "*", "")
signif_y<-tapply(table_SCFA$value, table_SCFA$variable, max )
table_SCFA$signif_y<-as.numeric(signif_y[table_SCFA$variable]) + 10 # plus 5 to put the sign higher

ggplot(data = table_SCFA, mapping = aes(x = variable, y = value, fill = Group)) +
  scale_fill_manual(values =  c("intervention" = alpha("red3", 0.5), "control" = alpha("#555555", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label = signif, y = signif_y, group = variable), 
            color = "darkgray", size = 8) +
  theme(axis.text.x = element_text(angle = -20, size = 9.5, vjust = 1, hjust = 0.05),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 11),
        panel.grid.minor.y = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.35),
        title = element_text(size = 10.5)) +
  scale_y_log10(breaks = c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50, seq(100, 500, 100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width = 0.65, size = 0.5, outlier.size = 1) + 
  labs(title = "SCFA concentration", x = "", fill = "", y = "Percentage (%)")
ggsave(filename = "SCFA_zoom_plot_TOT.png", height=4.5, width = 7, dpi=300)


# specific plot for MCFA
head(table_MCFA, n=2)
unique(table_MCFA$variable)
table_MCFA$signif<-ifelse(table_MCFA$variable %in% c("Total_MCFA"), "", "")   
signif_y<-tapply(table_MCFA$value, table_MCFA$variable, max )
table_MCFA$signif_y<-as.numeric(signif_y[table_MCFA$variable]) + 2 # plus 5 to put the sign higher

ggplot(data=table_MCFA, mapping=aes(x=variable, y=value, fill=Group)) +
  scale_fill_manual(values =  c("intervention" = alpha("red3", 0.5), "control" = alpha("#555555", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 500,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="MCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "MCFA_zoom_plot_TOT.png", height=4.5, width = 7, dpi=300)


# specific plot for LCFA
head(table_LCFA, n=2)
table_LCFA$signif<-ifelse(table_LCFA$variable %in% c("Total_LCFA"), "*", "")   
signif_y<-tapply(table_LCFA$value, table_LCFA$variable, max )
table_LCFA$signif_y<-as.numeric(signif_y[table_LCFA$variable]) + 10 # plus 5 to put the sign higher

# View(table_LCFA)
ggplot(data=table_LCFA, mapping=aes(x=variable, y=value, fill=Group)) +
  scale_fill_manual(values =  c("intervention" = alpha("red3", 0.5), "control" = alpha("#555555", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 700,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="LCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "LCFA_zoom_plot_TOT.png", height=4.5, width = 7, dpi=300)

