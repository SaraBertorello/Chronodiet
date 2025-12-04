################## ***CRONOTIPO COLIC ACIDS*** ################## 
# mmol/L

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
FFA <- read_excel("Colici_prove/Colic_acids.xlsx", sheet=1) 
colnames(FFA)
FFA$Type_diet
FFA[! colnames(FFA) %in% c("Time","Sample_name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("Time","Sample_name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Type_diet<-factor(FFA$Type_diet, levels = c("Control","Morning","Evening"))

##### ***T0*** ####
FFA<-subset.data.frame(FFA, Time=='T0')

primary <- FFA[,c("Time","Sample_name","Group","Type_diet","Cholic", "Chenodeoxycholic", "Glycocholic", "Glycochenodeoxycholic", "Taurocholic", "Taurochenodeoxycholic")]
secondary <-FFA[,c("Time","Sample_name","Group","Type_diet","Deoxycholic", "Glycodeoxycholic", "Taurodeoxycholic", "Lithocholic")]  
colnames(primary)
colnames(secondary)

######################### NORMALIZATIONS ##############################
primary_raw<-as.data.frame(primary)
secondary_raw<-as.data.frame(secondary)

# #### normalization through relative abundance in each sample, e.g. https://www.nature.com/articles/s41598-021-81420-3
# #### normalization of colic acids
primary_prop<-primary
primary_prop[ !colnames(primary_prop)%in% c("Time","Sample_name","Group","Type_diet")]<-primary_prop[!colnames(primary_prop)%in% c("Time","Sample_name","Group","Type_diet")]/rowSums(primary_prop[!colnames(primary_prop)%in% c("Time","Sample_name","Group","Type_diet")])
primary_prop<-as.data.frame(primary_prop)
# 
secondary_prop<-secondary
secondary_prop[ !colnames(secondary_prop)%in% c("Time","Sample_name","Group","Type_diet")]<-secondary_prop[!colnames(secondary_prop)%in% c("Time","Sample_name","Group","Type_diet")]/rowSums(secondary_prop[!colnames(secondary_prop)%in% c("Time","Sample_name","Group","Type_diet")])
secondary_prop<-as.data.frame(secondary_prop)
#tot
FFA_prop<-FFA
FFA_prop[ !colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")]<-FFA_prop[!colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")]/rowSums(FFA_prop[!colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")])
FFA_prop<-as.data.frame(FFA_prop)

####################### KRUSKAL WALLIS TESTs T0 ####################
dir.create("Colici_prove/Results")
library(FSA)
##### FFA raw abundances
list <- list("primary" = primary_raw, "secondary" = secondary_raw)
names(list)

for (y in 1:length(list)) {
  temp <- list[[y]]
  results <- NULL
  posthoc_results <- NULL
  
  for (x in colnames(temp)[!colnames(temp) %in% c("Time","Sample_name","Group","Type_diet")]) {
    test <- kruskal.test(temp[[x]] ~ temp$Type_diet)
    kruskal <- cbind(paste(x), test$statistic, test$p.value)
    results <- rbind.data.frame(results, kruskal)
    
    dunn <- dunnTest(temp[[x]] ~ temp$Type_diet, method = "bh")
    dunn_res <- cbind(variable = paste(x), dunn$res)
    posthoc_results <- rbind.data.frame(posthoc_results, dunn_res)
  }
  
  colnames(results) <- c(paste(names(list)[y]), "Chi-squared", "p.value")
  results$'p.adj(BH)' <- p.adjust(results$p.value, method = "BH")
  row.names(results) <- NULL
  write.csv2(results, file = paste0("Colici_prove/Results/", names(list)[y], "_raw_abund_KW.csv"), row.names = F)
  write.csv2(posthoc_results, file = paste0("Colici_prove/Results/", names(list)[y], "_raw_abund_Dunn.csv"), row.names = F)
  table <- melt(data = temp, id.vars = c("Time","Sample_name","Group","Type_diet"))
  if (y == 1) {
    table_SCFA <- table 
  }
  
  ggplot(data = table, mapping = aes(x = variable, y = value, fill = Type_diet)) +
    scale_fill_manual(values = c("Control" = "#CCCCCC", "Morning" = "#FFD700", "Evening" = '#4682B4')) +
    theme(axis.text.x = element_text(angle = -25, size = 8, vjust = 1, hjust = 0.05)) +
    geom_boxplot() +
    labs(title = paste("Raw", names(list)[y], "quantities", sep = " "), x = "", y = "Quantity")
  
  #ggsave(filename = paste0("FFA/", names(list)[y], "_raw_abundance.png"), height = 3, width = 6, dpi = 300)
  
  if (y == 2) {
    table_MCFA <- table
  }
}

####################### Kruskal-wallis TESTs ####################
###### loop for each colic acid
library(stringr)
### Stool
Stool<-FFA
Stool$Type_diet<-trimws(Stool$Type_diet)
head(Stool)
Stool<-as.data.frame(Stool)
Stool<-reshape::melt(Stool)
Stool$Type_diet
Stool$variable<-factor(Stool$variable, levels=unique(Stool$variable))
Stool$Type_diet<-factor(Stool$Type_diet, levels = c("Control","Morning","Evening"))

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
      scale_fill_manual(values=c("Control" = "#CCCCCC", "Morning" = "#FFD700", "Evening" = '#4682B4')) +
      scale_color_manual(values=c("Control" = "#CCCCCC", "Morning" = "#FFD700", "Evening" = '#4682B4')) +
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
      labs(y="Concentration (mmol/L)", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("Colici_prove/Results/Acids_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }
  
  if(length(which( grepl("*",test$p_adj_bh, fixed=T) )) == 0 ) { 
    
    kruskal_res<-kruskal.test(target$value~target$Type_diet)
    kruskal_res<-round(kruskal_res$p.value, 3)
    kruskal_res[kruskal_res==0]<-0.001 # if...
    
    ggboxplot(data=target, x="Type_diet", y="value", fill="Type_diet", 
              width=0.6, size= 0.15, alpha= 0.2, outlier.shape = NA) +
      scale_fill_manual(values=c("Control" = "#CCCCCC", "Morning" = "#FFD700", "Evening" = '#4682B4')) +
      scale_color_manual(values=c("Control" = "#CCCCCC", "Morning" = "#FFD700", "Evening" = '#4682B4')) +
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
      labs(y="Concentration (mmol/L)", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("Colici_prove/Results/Acids_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }
  
}

########### PCoA ########### 
list_PCoA <- list("primary" = FFA_prop, "secondary" = secondary_prop, "all"=FFA_prop)
p_values <- NULL
for (x in 1:length(list_PCoA)) {
  rm(data, coord)
  data <- as.data.frame(list_PCoA[[x]])
  row.names(data) <- data$Sample_name
  
  dist <- vegdist(data[, !colnames(data) %in% c("Time","Sample_name","Group","Type_diet")], method = "bray") 
  coord <- cmdscale(dist, k = 2, eig = TRUE) 

  coord_df <- as.data.frame(coord$points)
  colnames(coord_df) <- c("PC1", "PC2")
  coord_df$ID <- row.names(coord_df)
  coord_df$Type_diet <- data$Type_diet
  
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Type_diet)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("Control" = alpha("#CCCCCC", 0.5), "Morning" = alpha("#FFD700", 0.5), "Evening"=alpha('#4682B4',0.5))) +
    stat_ellipse(size = 0.2) +
    geom_text(aes(label = ID), size = 2, color = "black") +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of different type of diet\n with Bray-Curtis dissimilarity")
    )
  
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Type_diet)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("Control" = alpha("#CCCCCC", 0.5), "Morning" = alpha("#FFD700", 0.5), "Evening"=alpha('#4682B4',0.5))) +
    stat_ellipse(size = 0.2) +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of different type of diet\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("Colici_prove/Results/PCoA_BRAY_", names(list_PCoA)[x], "T0_noname.png", sep = "_"), 
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
write.csv2(p_values, file="Colici_prove/Results/p_values_BRAY_T0.csv", quote = F,row.names = T)


##### ***T1*** ####
FFA <- read_excel("Colici_prove/Colic_acids.xlsx", sheet=1) 
colnames(FFA)
FFA$Type_diet
FFA[! colnames(FFA) %in% c("Time","Sample_name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("Time","Sample_name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Type_diet<-factor(FFA$Type_diet, levels = c("Control","Morning","Evening"))

FFA<-subset.data.frame(FFA, Time=='T1')

primary <- FFA[,c("Time","Sample_name","Group","Type_diet","Cholic", "Chenodeoxycholic", "Glycocholic", "Glycochenodeoxycholic", "Taurocholic", "Taurochenodeoxycholic")]
secondary <-FFA[,c("Time","Sample_name","Group","Type_diet","Deoxycholic", "Glycodeoxycholic", "Taurodeoxycholic", "Lithocholic")]  
colnames(primary)
colnames(secondary)

######################### NORMALIZATIONS ##############################
primary_raw<-as.data.frame(primary)
secondary_raw<-as.data.frame(secondary)

# #### normalization through relative abundance in each sample, e.g. https://www.nature.com/articles/s41598-021-81420-3
# #### normalization of colic acids
primary_prop<-primary
primary_prop[ !colnames(primary_prop)%in% c("ID","Time","Sample_name","Group","Type_diet")]<-primary_prop[!colnames(primary_prop)%in% c("ID","Time","Sample_name","Group","Type_diet")]/rowSums(primary_prop[!colnames(primary_prop)%in% c("Time","Sample_name","Group","Type_diet")])
primary_prop<-as.data.frame(primary_prop)
# 
secondary_prop<-secondary
secondary_prop[ !colnames(secondary_prop)%in% c("ID","Time","Sample_name","Group","Type_diet")]<-secondary_prop[!colnames(secondary_prop)%in% c("Time","Sample_name","Group","Type_diet")]/rowSums(secondary_prop[!colnames(secondary_prop)%in% c("Time","Sample_name","Group","Type_diet")])
secondary_prop<-as.data.frame(secondary_prop)
#all
FFA_prop<-FFA
FFA_prop[ !colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")]<-FFA_prop[!colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")]/rowSums(FFA_prop[!colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")])
FFA_prop<-as.data.frame(FFA_prop)

####################### KRUSKAL WALLIS TESTs T1 ####################
dir.create("Colici_prove/Results")
library(FSA)
##### FFA raw abundances
list <- list("primary" = primary_raw, "secondary" = secondary_raw)
names(list)

for (y in 1:length(list)) {
  temp <- list[[y]]
  results <- NULL
  posthoc_results <- NULL
  
  for (x in colnames(temp)[!colnames(temp) %in% c("Time","Sample_name","Group","Type_diet")]) {
    # Kruskal-Wallis Test
    test <- kruskal.test(temp[[x]] ~ temp$Type_diet)
    kruskal <- cbind(paste(x), test$statistic, test$p.value)
    results <- rbind.data.frame(results, kruskal)
    
    dunn <- dunnTest(temp[[x]] ~ temp$Type_diet, method = "bh")
    dunn_res <- cbind(variable = paste(x), dunn$res)
    posthoc_results <- rbind.data.frame(posthoc_results, dunn_res)
  }
  
  colnames(results) <- c(paste(names(list)[y]), "Chi-squared", "p.value")
  results$'p.adj(BH)' <- p.adjust(results$p.value, method = "BH")
  row.names(results) <- NULL
  write.csv2(results, file = paste0("Colici_prove/Results/", names(list)[y], "_raw_abund_KW.csv"), row.names = F)
  write.csv2(posthoc_results, file = paste0("Colici_prove/Results/", names(list)[y], "_raw_abund_Dunn.csv"), row.names = F)
  
  table <- melt(data = temp, id.vars = c("Time","Sample_name","Group","Type_diet"))
  if (y == 1) {
    table_SCFA <- table 
  }
  
  ggplot(data = table, mapping = aes(x = variable, y = value, fill = Type_diet)) +
    scale_fill_manual(values = c("Control" = "#CCCCCC", "Morning" = "#FFD700", "Evening" = '#4682B4')) +
    theme(axis.text.x = element_text(angle = -25, size = 8, vjust = 1, hjust = 0.05)) +
    geom_boxplot() +
    labs(title = paste("Raw", names(list)[y], "quantities", sep = " "), x = "", y = "Quantity")
  
  #ggsave(filename = paste0("FFA/", names(list)[y], "_raw_abundance.png"), height = 3, width = 6, dpi = 300)

  if (y == 2) {
    table_MCFA <- table
  }
}

####################### Kruskal-wallis TESTs ####################
###### loop for each colic acid
library(stringr)
### Stool
Stool<-FFA
Stool$Type_diet<-trimws(Stool$Type_diet)
head(Stool)
Stool<-as.data.frame(Stool)
Stool<-reshape::melt(Stool)
Stool$Type_diet
Stool$variable<-factor(Stool$variable, levels=unique(Stool$variable))
Stool$Type_diet<-factor(Stool$Type_diet, levels = c("Control","Morning","Evening"))

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
      scale_fill_manual(values=c("Control" = "#CCCCCC", "Morning" = "#FFD700", "Evening" = '#4682B4')) +
      scale_color_manual(values=c("Control" = "#CCCCCC", "Morning" = "#FFD700", "Evening" = '#4682B4')) +
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
      labs(y="Concentration (mmol/L)", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("Colici_prove/Results/Acids_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }
  
  if(length(which( grepl("*",test$p_adj_bh, fixed=T) )) == 0 ) { # otherwise, if there are not pairwise results...
    
    kruskal_res<-kruskal.test(target$value~target$Type_diet)
    kruskal_res<-round(kruskal_res$p.value, 3)
    kruskal_res[kruskal_res==0]<-0.001 # if...
    
    ggboxplot(data=target, x="Type_diet", y="value", fill="Type_diet", 
              width=0.6, size= 0.15, alpha= 0.2, outlier.shape = NA) +
      scale_fill_manual(values=c("Control" = "#CCCCCC", "Morning" = "#FFD700", "Evening" = '#4682B4')) +
      scale_color_manual(values=c("Control" = "#CCCCCC", "Morning" = "#FFD700", "Evening" = '#4682B4')) +
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
      labs(y="Concentration (mmol/L)", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("Colici_prove/Results/Acids_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }
  
}

########### PCoA ########### 
list_PCoA <- list("primary" = primary_prop, "secondary" = secondary_prop, "all"=FFA_prop)
p_values <- NULL
for (x in 1:length(list_PCoA)) {
  rm(data, coord)
  data <- as.data.frame(list_PCoA[[x]])
  row.names(data) <- data$Sample_name
  
  dist <- vegdist(data[, !colnames(data) %in% c("Time","Sample_name","Group","Type_diet")], method = "bray") 
  coord <- cmdscale(dist, k = 2, eig = TRUE) # PCoA
  
  coord_df <- as.data.frame(coord$points)
  colnames(coord_df) <- c("PC1", "PC2")
  coord_df$ID <- row.names(coord_df)
  coord_df$Type_diet <- data$Type_diet
  
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Type_diet)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("Control" = alpha("#CCCCCC", 0.5), "Morning" = alpha("#FFD700", 0.5), "Evening"=alpha('#4682B4',0.5))) +
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
    scale_color_manual(values = c("Control" = alpha("#CCCCCC", 0.5), "Morning" = alpha("#FFD700", 0.5), "Evening"=alpha('#4682B4',0.5))) +
    stat_ellipse(size = 0.2) +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of different type of diet\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("Colici_prove/Results/PCoA_BRAY_", names(list_PCoA)[x], "T0_noname.png", sep = "_"), 
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
write.csv2(p_values, file="Colici_prove/Results/p_values_BRAY_T1.csv", quote = F,row.names = T)

########## ***TYPE DIET*** ##########
##### ***MORNING*** ####
FFA <- read_excel("Colici_prove/Colic_acids.xlsx", sheet=1) 
colnames(FFA)
FFA$Type_diet
FFA[! colnames(FFA) %in% c("Time","Sample_name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("Time","Sample_name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Type_diet<-factor(FFA$Type_diet, levels = c("Control","Morning","Evening"))

FFA<-subset.data.frame(FFA, Type_diet=='Morning')

primary <- FFA[,c("Time","Sample_name","Group","Type_diet","Cholic", "Chenodeoxycholic", "Glycocholic", "Glycochenodeoxycholic", "Taurocholic", "Taurochenodeoxycholic")]
secondary <-FFA[,c("Time","Sample_name","Group","Type_diet","Deoxycholic", "Glycodeoxycholic", "Taurodeoxycholic", "Lithocholic")]  
colnames(primary)
colnames(secondary)

######################### NORMALIZATIONS ##############################
primary_raw<-as.data.frame(primary)
secondary_raw<-as.data.frame(secondary)

# #### normalization through relative abundance in each sample, e.g. https://www.nature.com/articles/s41598-021-81420-3
# #### normalization of colic acids
primary_prop<-primary
primary_prop[ !colnames(primary_prop)%in% c("ID","Time","Sample_name","Group","Type_diet")]<-primary_prop[!colnames(primary_prop)%in% c("ID","Time","Sample_name","Group","Type_diet")]/rowSums(primary_prop[!colnames(primary_prop)%in% c("Time","Sample_name","Group","Type_diet")])
primary_prop<-as.data.frame(primary_prop)
# 
secondary_prop<-secondary
secondary_prop[ !colnames(secondary_prop)%in% c("ID","Time","Sample_name","Group","Type_diet")]<-secondary_prop[!colnames(secondary_prop)%in% c("Time","Sample_name","Group","Type_diet")]/rowSums(secondary_prop[!colnames(secondary_prop)%in% c("Time","Sample_name","Group","Type_diet")])
secondary_prop<-as.data.frame(secondary_prop)
#all
FFA_prop<-FFA
FFA_prop[ !colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")]<-FFA_prop[!colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")]/rowSums(FFA_prop[!colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")])
FFA_prop<-as.data.frame(FFA_prop)

####################### MANN WITHNEY TESTs ####################
##### FFA raw abundances
list<-list("all"=FFA)
FFA
names(list)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("Time","Sample_name","Group","Type_diet")]) {
    test<-wilcox.test(temp[[x]]~temp$Time)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  write.csv2(results, file=paste0("Colici_prove/Results/",names(list)[y],"_raw_abund_Mann_Whit.csv"), row.names = F)
  
  table <- melt(data = temp, id.vars = c("Time","Sample_name","Group","Type_diet"))
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Time)) +
    scale_fill_manual(values =c("T0" = "#FFD700", "T1" = "#FFB347")) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  #ggsave(filename = paste0("FFA/",names(list)[y],"_raw_abundance.png"), height=3, width = 6, dpi=300)
  
}


# specific plot for ALL colic acids
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Taurochenodeoxycholic", "Taurocholic","Taurodeoxycholic","Glycocholic"), "*", "")
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
  scale_y_log10(breaks = c(0.02, 0.5, 1, 2.5, 5, 10, 25, 50, seq(100, 500, 100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width = 0.65, size = 0.5, outlier.size = 1) + 
  labs(title = "Cholic acids concentration", x = "", fill = "", y = "Concentration (mmol/L)")
ggsave(filename = "Colici_prove/Results/Morning_zoom_plot.png", height=4.5, width = 7, dpi=300)


########### PCoA ########### 
list_PCoA <- list("primary" = primary_prop, "secondary" = secondary_prop, "all"=FFA_prop)
p_values <- NULL
for (x in 1:length(list_PCoA)) {
  rm(data, coord)
  data <- as.data.frame(list_PCoA[[x]])
  row.names(data) <- data$Sample_name
  
  dist <- vegdist(data[, !colnames(data) %in% c("Time","Sample_name","Group","Type_diet")], method = "bray") # calcolare dissimilarità Bray-Curtis
  coord <- cmdscale(dist, k = 2, eig = TRUE) # PCoA
  
  coord_df <- as.data.frame(coord$points)
  colnames(coord_df) <- c("PC1", "PC2")
  coord_df$ID <- row.names(coord_df)
  coord_df$Time <- data$Time
  
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Time)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("T0" = alpha("#FFD700", 0.5), "T1" = alpha("#FFB347", 0.5))) +
    stat_ellipse(size = 0.2) +
    geom_text(aes(label = ID), size = 2, color = "black") +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of different type of diet\n with Bray-Curtis dissimilarity")
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
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of different type of diet\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("Colici_prove/Results/PCoA_BRAY_", names(list_PCoA)[x], "T0_noname.png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )                         
  

  p_value_disp <- permutest(betadisper(dist, group = data$Time), permutations = 9999)
  p_value_diver <- adonis(dist ~ Time, data = data, permutations = 9999)
  p_value_disp <- as.data.frame(p_value_disp$tab[1, ])
  colnames(p_value_disp) <- colnames(p_value_diver$aov.tab[1, ])
  new_values <- rbind(p_value_diver$aov.tab[1, ], p_value_disp)
  row.names(new_values) <- c(paste("Beta_diversity", names(list_PCoA)[x], sep = "_"),
                             paste("Beta_dispersion", names(list_PCoA)[x], sep = "_"))
  p_values <- rbind(p_values, new_values)
  
}

p_values$adjusted_BH_p_values<-p.adjust(p_values$`Pr(>F)`, method = "BH")
p_values
write.csv2(p_values, file="Colici_prove/Results/p_values_BRAY_Morning.csv", quote = F,row.names = T)


##### ***EVENING*** ####
FFA <- read_excel("Colici_prove/Colic_acids.xlsx", sheet=1) 
colnames(FFA)
FFA$Type_diet
FFA[! colnames(FFA) %in% c("Time","Sample_name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("Time","Sample_name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Type_diet<-factor(FFA$Type_diet, levels = c("Control","Morning","Evening"))

FFA<-subset.data.frame(FFA, Type_diet=='Evening')

primary <- FFA[,c("Time","Sample_name","Group","Type_diet","Cholic", "Chenodeoxycholic", "Glycocholic", "Glycochenodeoxycholic", "Taurocholic", "Taurochenodeoxycholic")]
secondary <-FFA[,c("Time","Sample_name","Group","Type_diet","Deoxycholic", "Glycodeoxycholic", "Taurodeoxycholic", "Lithocholic")]  
colnames(primary)
colnames(secondary)

######################### NORMALIZATIONS ##############################
primary_raw<-as.data.frame(primary)
secondary_raw<-as.data.frame(secondary)

# #### normalization through relative abundance in each sample, e.g. https://www.nature.com/articles/s41598-021-81420-3
# #### normalization of colic acids
primary_prop<-primary
primary_prop[ !colnames(primary_prop)%in% c("ID","Time","Sample_name","Group","Type_diet")]<-primary_prop[!colnames(primary_prop)%in% c("ID","Time","Sample_name","Group","Type_diet")]/rowSums(primary_prop[!colnames(primary_prop)%in% c("Time","Sample_name","Group","Type_diet")])
primary_prop<-as.data.frame(primary_prop)
# 
secondary_prop<-secondary
secondary_prop[ !colnames(secondary_prop)%in% c("ID","Time","Sample_name","Group","Type_diet")]<-secondary_prop[!colnames(secondary_prop)%in% c("Time","Sample_name","Group","Type_diet")]/rowSums(secondary_prop[!colnames(secondary_prop)%in% c("Time","Sample_name","Group","Type_diet")])
secondary_prop<-as.data.frame(secondary_prop)
#all
FFA_prop<-FFA
FFA_prop[ !colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")]<-FFA_prop[!colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")]/rowSums(FFA_prop[!colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")])
FFA_prop<-as.data.frame(FFA_prop)

####################### MANN WITHNEY TESTs ####################
##### FFA raw abundances
list<-list("all"=FFA)
FFA
names(list)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("Time","Sample_name","Group","Type_diet")]) {
    test<-wilcox.test(temp[[x]]~temp$Time)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  write.csv2(results, file=paste0("Colici_prove/Results/",names(list)[y],"_raw_abund_Mann_Whit.csv"), row.names = F)
  
  table <- melt(data = temp, id.vars = c("Time","Sample_name","Group","Type_diet"))
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Time)) +
    scale_fill_manual(values =c("T0" = "#4682B4", "T1" = "#0047AB")) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  #ggsave(filename = paste0("FFA/",names(list)[y],"_raw_abundance.png"), height=3, width = 6, dpi=300)
  
}


# specific plot for ALL colic acids
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Taurochenodeoxycholic", "Taurocholic","Taurodeoxycholic"), "*", "")
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
  scale_y_log10(breaks = c(0.02, 0.5, 1, 2.5, 5, 10, 25, 50, seq(100, 500, 100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width = 0.65, size = 0.5, outlier.size = 1) + 
  labs(title = "Cholic acids concentration", x = "", fill = "", y = "Concentration (mmol/L)")
ggsave(filename = "Colici_prove/Results/Evening_zoom_plot.png", height=4.5, width = 7, dpi=300)


########### PCoA ########### 
list_PCoA <- list("primary" = primary_prop, "secondary" = secondary_prop, "all"=FFA_prop)
p_values <- NULL
for (x in 1:length(list_PCoA)) {
  rm(data, coord)
  data <- as.data.frame(list_PCoA[[x]])
  row.names(data) <- data$Sample_name
  
  dist <- vegdist(data[, !colnames(data) %in% c("Time","Sample_name","Group","Type_diet")], method = "bray") # calcolare dissimilarità Bray-Curtis
  coord <- cmdscale(dist, k = 2, eig = TRUE) # PCoA
  
  # Convertire le coordinate in un dataframe e aggiungere ID e Time
  coord_df <- as.data.frame(coord$points)
  colnames(coord_df) <- c("PC1", "PC2")
  coord_df$ID <- row.names(coord_df)
  coord_df$Time <- data$Time
  
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Time)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("T0" = alpha("#4682B4", 0.5), "T1" = alpha("#0047AB", 0.5))) +
    stat_ellipse(size = 0.2) +
    geom_text(aes(label = ID), size = 2, color = "black") +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of different type of diet\n with Bray-Curtis dissimilarity")
    )
  
  # again but with no names
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Time)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("T0" = alpha("#4682B4", 0.5), "T1" = alpha("#0047AB", 0.5))) +
    stat_ellipse(size = 0.2) +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of different type of diet\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("Colici_prove/Results/PCoA_BRAY_", names(list_PCoA)[x], "T0_noname.png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )                         
  
  # Calcolare la dispersione beta per 'Time'
  p_value_disp <- permutest(betadisper(dist, group = data$Time), permutations = 9999)
  p_value_diver <- adonis(dist ~ Time, data = data, permutations = 9999)
  p_value_disp <- as.data.frame(p_value_disp$tab[1, ])
  colnames(p_value_disp) <- colnames(p_value_diver$aov.tab[1, ])
  new_values <- rbind(p_value_diver$aov.tab[1, ], p_value_disp)
  row.names(new_values) <- c(paste("Beta_diversity", names(list_PCoA)[x], sep = "_"),
                             paste("Beta_dispersion", names(list_PCoA)[x], sep = "_"))
  p_values <- rbind(p_values, new_values)
  
}

p_values$adjusted_BH_p_values<-p.adjust(p_values$`Pr(>F)`, method = "BH")
p_values
write.csv2(p_values, file="Colici_prove/Results/p_values_BRAY_Evening.csv", quote = F,row.names = T)

##### ***CONTROL*** ####
FFA <- read_excel("Colici_prove/Colic_acids.xlsx", sheet=1) 
colnames(FFA)
FFA$Type_diet
FFA[! colnames(FFA) %in% c("Time","Sample_name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("Time","Sample_name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Type_diet<-factor(FFA$Type_diet, levels = c("Control","Morning","Evening"))

FFA<-subset.data.frame(FFA, Group=='Control')

primary <- FFA[,c("Time","Sample_name","Group","Type_diet","Cholic", "Chenodeoxycholic", "Glycocholic", "Glycochenodeoxycholic", "Taurocholic", "Taurochenodeoxycholic")]
secondary <-FFA[,c("Time","Sample_name","Group","Type_diet","Deoxycholic", "Glycodeoxycholic", "Taurodeoxycholic", "Lithocholic")]  
colnames(primary)
colnames(secondary)

######################### NORMALIZATIONS ##############################
primary_raw<-as.data.frame(primary)
secondary_raw<-as.data.frame(secondary)

# #### normalization through relative abundance in each sample, e.g. https://www.nature.com/articles/s41598-021-81420-3
# #### normalization of colic acids
primary_prop<-primary
primary_prop[ !colnames(primary_prop)%in% c("ID","Time","Sample_name","Group","Type_diet")]<-primary_prop[!colnames(primary_prop)%in% c("ID","Time","Sample_name","Group","Type_diet")]/rowSums(primary_prop[!colnames(primary_prop)%in% c("Time","Sample_name","Group","Type_diet")])
primary_prop<-as.data.frame(primary_prop)
# 
secondary_prop<-secondary
secondary_prop[ !colnames(secondary_prop)%in% c("ID","Time","Sample_name","Group","Type_diet")]<-secondary_prop[!colnames(secondary_prop)%in% c("Time","Sample_name","Group","Type_diet")]/rowSums(secondary_prop[!colnames(secondary_prop)%in% c("Time","Sample_name","Group","Type_diet")])
secondary_prop<-as.data.frame(secondary_prop)
#all
FFA_prop<-FFA
FFA_prop[ !colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")]<-FFA_prop[!colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")]/rowSums(FFA_prop[!colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")])
FFA_prop<-as.data.frame(FFA_prop)

####################### MANN WITHNEY TESTs ####################
##### FFA raw abundances
list<-list("all"=FFA)
FFA
names(list)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("Time","Sample_name","Group","Type_diet")]) {
    test<-wilcox.test(temp[[x]]~temp$Time)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  write.csv2(results, file=paste0("Colici_prove/Results/",names(list)[y],"_raw_abund_Mann_Whit.csv"), row.names = F)
  
  table <- melt(data = temp, id.vars = c("Time","Sample_name","Group","Type_diet"))
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Time)) +
    scale_fill_manual(values =c("T0" = "#CCCCCC", "T1" = "#555555")) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  #ggsave(filename = paste0("FFA/",names(list)[y],"_raw_abundance.png"), height=3, width = 6, dpi=300)
  
}


# specific plot for ALL colic acids
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Taurochenodeoxycholic", "Taurocholic","Taurodeoxycholic"), "*", "")
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
  scale_y_log10(breaks = c(0.02, 0.5, 1, 2.5, 5, 10, 25, 50, seq(100, 500, 100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width = 0.65, size = 0.5, outlier.size = 1) + 
  labs(title = "Cholic acids concentration", x = "", fill = "", y = "Concentration (mmol/L)")
ggsave(filename = "Colici_prove/Results/Control_zoom_plot.png", height=4.5, width = 7, dpi=300)


########### PCoA ########### 
list_PCoA <- list("primary" = primary_prop, "secondary" = secondary_prop, "all"=FFA_prop)
p_values <- NULL
for (x in 1:length(list_PCoA)) {
  rm(data, coord)
  data <- as.data.frame(list_PCoA[[x]])
  row.names(data) <- data$Sample_name
  
  dist <- vegdist(data[, !colnames(data) %in% c("Time","Sample_name","Group","Type_diet")], method = "bray") # calcolare dissimilarità Bray-Curtis
  coord <- cmdscale(dist, k = 2, eig = TRUE) # PCoA
  
  # Convertire le coordinate in un dataframe e aggiungere ID e Time
  coord_df <- as.data.frame(coord$points)
  colnames(coord_df) <- c("PC1", "PC2")
  coord_df$ID <- row.names(coord_df)
  coord_df$Time <- data$Time
  
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Time)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("T0" = alpha("#CCCCCC", 0.5), "T1" = alpha("#555555", 0.5))) +
    stat_ellipse(size = 0.2) +
    geom_text(aes(label = ID), size = 2, color = "black") +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of different type of diet\n with Bray-Curtis dissimilarity")
    )
  
  # again but with no names
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Time)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("T0" = alpha("#CCCCCC", 0.5), "T1" = alpha("#555555", 0.5))) +
    stat_ellipse(size = 0.2) +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of different type of diet\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("Colici_prove/Results/PCoA_BRAY_", names(list_PCoA)[x], "T0_noname.png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )                         
  
  p_value_disp <- permutest(betadisper(dist, group = data$Time), permutations = 9999)
  p_value_diver <- adonis(dist ~ Time, data = data, permutations = 9999)
  p_value_disp <- as.data.frame(p_value_disp$tab[1, ])
  colnames(p_value_disp) <- colnames(p_value_diver$aov.tab[1, ])
  new_values <- rbind(p_value_diver$aov.tab[1, ], p_value_disp)
  row.names(new_values) <- c(paste("Beta_diversity", names(list_PCoA)[x], sep = "_"),
                             paste("Beta_dispersion", names(list_PCoA)[x], sep = "_"))
  p_values <- rbind(p_values, new_values)
  
}

p_values$adjusted_BH_p_values<-p.adjust(p_values$`Pr(>F)`, method = "BH")
p_values
write.csv2(p_values, file="Colici_prove/Results/p_values_BRAY_Control.csv", quote = F,row.names = T)

##### ***INTERVENTION*** ####
FFA <- read_excel("Colici_prove/Colic_acids.xlsx", sheet=1) 
colnames(FFA)
FFA$Type_diet
FFA[! colnames(FFA) %in% c("Time","Sample_name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("Time","Sample_name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Type_diet<-factor(FFA$Type_diet, levels = c("Control","Morning","Evening"))

FFA<-subset.data.frame(FFA, Group=='Intervention')

primary <- FFA[,c("Time","Sample_name","Group","Type_diet","Cholic", "Chenodeoxycholic", "Glycocholic", "Glycochenodeoxycholic", "Taurocholic", "Taurochenodeoxycholic")]
secondary <-FFA[,c("Time","Sample_name","Group","Type_diet","Deoxycholic", "Glycodeoxycholic", "Taurodeoxycholic", "Lithocholic")]  
colnames(primary)
colnames(secondary)

######################### NORMALIZATIONS ##############################
primary_raw<-as.data.frame(primary)
secondary_raw<-as.data.frame(secondary)

# #### normalization through relative abundance in each sample, e.g. https://www.nature.com/articles/s41598-021-81420-3
# #### normalization of colic acids
primary_prop<-primary
primary_prop[ !colnames(primary_prop)%in% c("ID","Time","Sample_name","Group","Type_diet")]<-primary_prop[!colnames(primary_prop)%in% c("ID","Time","Sample_name","Group","Type_diet")]/rowSums(primary_prop[!colnames(primary_prop)%in% c("Time","Sample_name","Group","Type_diet")])
primary_prop<-as.data.frame(primary_prop)
# 
secondary_prop<-secondary
secondary_prop[ !colnames(secondary_prop)%in% c("ID","Time","Sample_name","Group","Type_diet")]<-secondary_prop[!colnames(secondary_prop)%in% c("Time","Sample_name","Group","Type_diet")]/rowSums(secondary_prop[!colnames(secondary_prop)%in% c("Time","Sample_name","Group","Type_diet")])
secondary_prop<-as.data.frame(secondary_prop)
#all
FFA_prop<-FFA
FFA_prop[ !colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")]<-FFA_prop[!colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")]/rowSums(FFA_prop[!colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")])
FFA_prop<-as.data.frame(FFA_prop)

####################### MANN WITHNEY TESTs ####################
##### FFA raw abundances
list<-list("all"=FFA)
FFA
names(list)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("Time","Sample_name","Group","Type_diet")]) {
    test<-wilcox.test(temp[[x]]~temp$Time)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  write.csv2(results, file=paste0("Colici_prove/Results/",names(list)[y],"_raw_abund_Mann_Whit.csv"), row.names = F)
  
  table <- melt(data = temp, id.vars = c("Time","Sample_name","Group","Type_diet"))
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Time)) +
    scale_fill_manual(values =c("T1" = "red3", "T0" = "#FF6666"))  +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  #ggsave(filename = paste0("FFA/",names(list)[y],"_raw_abundance.png"), height=3, width = 6, dpi=300)
  
}


# specific plot for ALL colic acids
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Taurochenodeoxycholic", "Taurocholic","Taurodeoxycholic","Glycocholic"), "*", "")
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
  scale_y_log10(breaks = c(0.02, 0.5, 1, 2.5, 5, 10, 25, 50, seq(100, 500, 100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width = 0.65, size = 0.5, outlier.size = 1) + 
  labs(title = "Cholic acids concentration", x = "", fill = "", y = "Concentration (mmol/L)")
ggsave(filename = "Colici_prove/Results/intervention_zoom_plot.png", height=4.5, width = 7, dpi=300)


########### PCoA ########### 
list_PCoA <- list("primary" = primary_prop, "secondary" = secondary_prop, "all"=FFA_prop)
p_values <- NULL
for (x in 1:length(list_PCoA)) {
  rm(data, coord)
  data <- as.data.frame(list_PCoA[[x]])
  row.names(data) <- data$Sample_name
  
  dist <- vegdist(data[, !colnames(data) %in% c("Time","Sample_name","Group","Type_diet")], method = "bray") 
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
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of different type of diet\n with Bray-Curtis dissimilarity")
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
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of different type of diet\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("Colici_prove/Results/PCoA_BRAY_", names(list_PCoA)[x], "T0_noname.png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )                         
  
  p_value_disp <- permutest(betadisper(dist, group = data$Time), permutations = 9999)
  p_value_diver <- adonis(dist ~ Time, data = data, permutations = 9999)
  p_value_disp <- as.data.frame(p_value_disp$tab[1, ])
  colnames(p_value_disp) <- colnames(p_value_diver$aov.tab[1, ])
  new_values <- rbind(p_value_diver$aov.tab[1, ], p_value_disp)
  row.names(new_values) <- c(paste("Beta_diversity", names(list_PCoA)[x], sep = "_"),
                             paste("Beta_dispersion", names(list_PCoA)[x], sep = "_"))
  p_values <- rbind(p_values, new_values)
  
}

p_values$adjusted_BH_p_values<-p.adjust(p_values$`Pr(>F)`, method = "BH")
p_values
write.csv2(p_values, file="Colici_prove/Results/p_values_BRAY_Intervention.csv", quote = F,row.names = T)

##### ***CONTROL VS INTERVENTION*** ####
FFA <- read_excel("Colici_prove/Colic_acids.xlsx", sheet=1) 
colnames(FFA)
FFA$Type_diet
FFA[! colnames(FFA) %in% c("Time","Sample_name","Group","Type_diet")]<-apply(FFA[!colnames(FFA)%in% c("Time","Sample_name","Group","Type_diet")], MARGIN = 2, as.numeric)
FFA$Type_diet<-factor(FFA$Type_diet, levels = c("Control","Morning","Evening"))

FFA<-subset.data.frame(FFA, Time=='T1')

primary <- FFA[,c("Time","Sample_name","Group","Type_diet","Cholic", "Chenodeoxycholic", "Glycocholic", "Glycochenodeoxycholic", "Taurocholic", "Taurochenodeoxycholic")]
secondary <-FFA[,c("Time","Sample_name","Group","Type_diet","Deoxycholic", "Glycodeoxycholic", "Taurodeoxycholic", "Lithocholic")]  
colnames(primary)
colnames(secondary)

######################### NORMALIZATIONS ##############################
primary_raw<-as.data.frame(primary)
secondary_raw<-as.data.frame(secondary)

# #### normalization through relative abundance in each sample, e.g. https://www.nature.com/articles/s41598-021-81420-3
# #### normalization of colic acids
primary_prop<-primary
primary_prop[ !colnames(primary_prop)%in% c("ID","Time","Sample_name","Group","Type_diet")]<-primary_prop[!colnames(primary_prop)%in% c("ID","Time","Sample_name","Group","Type_diet")]/rowSums(primary_prop[!colnames(primary_prop)%in% c("Time","Sample_name","Group","Type_diet")])
primary_prop<-as.data.frame(primary_prop)
# 
secondary_prop<-secondary
secondary_prop[ !colnames(secondary_prop)%in% c("ID","Time","Sample_name","Group","Type_diet")]<-secondary_prop[!colnames(secondary_prop)%in% c("Time","Sample_name","Group","Type_diet")]/rowSums(secondary_prop[!colnames(secondary_prop)%in% c("Time","Sample_name","Group","Type_diet")])
secondary_prop<-as.data.frame(secondary_prop)
#all
FFA_prop<-FFA
FFA_prop[ !colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")]<-FFA_prop[!colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")]/rowSums(FFA_prop[!colnames(FFA_prop)%in% c("Time","Sample_name","Group","Type_diet")])
FFA_prop<-as.data.frame(FFA_prop)

####################### MANN WITHNEY TESTs ####################
##### FFA raw abundances
list<-list("all"=FFA)
FFA
names(list)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("Time","Sample_name","Group","Type_diet")]) {
    test<-wilcox.test(temp[[x]]~temp$Group)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  write.csv2(results, file=paste0("Colici_prove/Results/",names(list)[y],"_raw_abund_Mann_Whit.csv"), row.names = F)
  
  table <- melt(data = temp, id.vars = c("Time","Sample_name","Group","Type_diet"))
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Group)) +
    scale_fill_manual(values =c("Intervention" = "red3", "Control" = '#555555')) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  #ggsave(filename = paste0("FFA/",names(list)[y],"_raw_abundance.png"), height=3, width = 6, dpi=300)
  
}


# specific plot for ALL colic acids
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Taurochenodeoxycholic", "Taurocholic","Taurodeoxycholic","Glycocholic"), "", "")
signif_y<-tapply(table_SCFA$value, table_SCFA$variable, max )
table_SCFA$signif_y<-as.numeric(signif_y[table_SCFA$variable]) + 10 # plus 5 to put the sign higher

ggplot(data = table_SCFA, mapping = aes(x = variable, y = value, fill = Group)) +
  scale_fill_manual(values = c("Intervention" = alpha("red3", 0.5), "Control" = alpha("#555555", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label = signif, y = signif_y, group = variable), 
            color = "darkgray", size = 8) +
  theme(axis.text.x = element_text(angle = -20, size = 9.5, vjust = 1, hjust = 0.05),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 11),
        panel.grid.minor.y = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.35),
        title = element_text(size = 10.5)) +
  scale_y_log10(breaks = c(0.02, 0.5, 1, 2.5, 5, 10, 25, 50, seq(100, 500, 100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width = 0.65, size = 0.5, outlier.size = 1) + 
  labs(title = "Cholic acids concentration", x = "", fill = "", y = "Concentration (mmol/L)")
ggsave(filename = "Colici_prove/Results/intervention_cnt_zoom_plot.png", height=4.5, width = 7, dpi=300)


########### PCoA ########### 
list_PCoA <- list("primary" = primary_prop, "secondary" = secondary_prop, "all"=FFA_prop)
p_values <- NULL
for (x in 1:length(list_PCoA)) {
  rm(data, coord)
  data <- as.data.frame(list_PCoA[[x]])
  row.names(data) <- data$Sample_name
  
  dist <- vegdist(data[, !colnames(data) %in% c("Time","Sample_name","Group","Type_diet")], method = "bray") # calcolare dissimilarità Bray-Curtis
  coord <- cmdscale(dist, k = 2, eig = TRUE) # PCoA
  
  # Convertire le coordinate in un dataframe e aggiungere ID e Group
  coord_df <- as.data.frame(coord$points)
  colnames(coord_df) <- c("PC1", "PC2")
  coord_df$ID <- row.names(coord_df)
  coord_df$Group <- data$Group
  
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Group)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("Control" = alpha("#555555", 0.5), "Intervention" = alpha("red3", 0.5))) +
    stat_ellipse(size = 0.2) +
    geom_text(aes(label = ID), size = 2, color = "black") +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of different type of diet\n with Bray-Curtis dissimilarity")
    )
  
  # again but with no names
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Group)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("Control" = alpha("#555555", 0.5), "Intervention" = alpha("red3", 0.5))) +
    stat_ellipse(size = 0.2) +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of different type of diet\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("Colici_prove/Results/PCoA_BRAY_", names(list_PCoA)[x], "T0_noname.png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )                         
  
  p_value_disp <- permutest(betadisper(dist, group = data$Group), permutations = 9999)
  p_value_diver <- adonis(dist ~ Group, data = data, permutations = 9999)
  p_value_disp <- as.data.frame(p_value_disp$tab[1, ])
  colnames(p_value_disp) <- colnames(p_value_diver$aov.tab[1, ])
  new_values <- rbind(p_value_diver$aov.tab[1, ], p_value_disp)
  row.names(new_values) <- c(paste("Beta_diversity", names(list_PCoA)[x], sep = "_"),
                             paste("Beta_dispersion", names(list_PCoA)[x], sep = "_"))
  p_values <- rbind(p_values, new_values)
  
}

p_values$adjusted_BH_p_values<-p.adjust(p_values$`Pr(>F)`, method = "BH")
p_values
write.csv2(p_values, file="Colici_prove/Results/p_values_BRAY_Intervention_cnt.csv", quote = F,row.names = T)

