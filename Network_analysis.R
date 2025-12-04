###################### CHRONODIET NETWORK ANALYSIS ###################### 
{ 
  library("phyloseq")
  library("ggplot2")
  library("ggvenn")
  library("egg")
  library("ggh4x")
  library("ggbreak") 
  library("vegan")
  library("stringr")
  library("xlsx")  
  library("qiime2R")
  library(phyloseq)
  library(ggplot2)
  library(vegan)
  library(readxl)
  library(dplyr)
  library(DESeq2)
  library(ggpubr)
  library(devtools)
  library(tidyr)
  library(qiime2R)
  library(dunn.test)
  library(pairwiseAdonis)
  library(cluster)
  library(ggh4x) 
  library(NetCoMi)
  library(phyloseq)
  library(SpiecEasi)
}
options(scipen=100)
#setwd('../../Cronotipo/Cronotipo_2025')
load('data_Cronotipo_filtered.RData')

###### ON TOP 20 GENERA MORNING ###### 
sample_data(data)
morning_evening<-subset_samples(data, Time=='T1')
morning_evening<-subset_samples(morning_evening, Type_diet=='morning')
sample_data(morning_evening)

data.phy <- tax_glom(morning_evening, taxrank = "Genus", NArm = F)
data.phy.prop <- transform_sample_counts(data.phy, function(x) x/sum(x)*100)
top5_Genus <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:20] 
prune.dat_top5_Genus <- prune_taxa(top5_Genus, data.phy.prop)
otu_table(prune.dat_top5_Genus)

sample_data(data.phy.prop) #morning e top 15 genera
data.phy.prop<-prune.dat_top5_Genus
sample_data(data.phy.prop)
#data.phy.prop <- subset_taxa(data.phy.prop, !grepl("uncultured", tax_table(data.phy.prop)[, "Genus"])) 

###### NETWORK ANALYSIS ###### 
amgut_genus <- data.phy.prop
tax_table(data.phy.prop)

taxtab <- as(tax_table(amgut_genus), "matrix")

amgut_genus_renamed <- renameTaxa(amgut_genus, 
                                  pat = "<name>", 
                                  substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Genus")
tax_table(amgut_genus_renamed)
# Network construction and analysis
net_genus <- netConstruct(amgut_genus_renamed,
                          taxRank = "Genus",
                          measure = "spearman",
                          zeroMethod = "multRepl",
                          normMethod = "clr",
                          sparsMethod = "threshold",
                          thresh = 0.3,
                          verbose = 3)
props_genus <- netAnalyze(net_genus, clustMethod = "cluster_fast_greedy")

#CENTRALITY
centro<-props_genus$centralities$eigenv1 
df<-as.data.frame(centro)

library(dplyr)
df_ordinato <- df %>%
  arrange(desc(centro))
                                         
print(df_ordinato) 
props_genus$compSize1

#AVERAGE PATH LENGHT
#install.packages('igraph')
library(igraph)

names(net_genus)

# --- Metric ---
edge_list <- net_genus$edgelist1  
edge_list_simple <- edge_list[, 1:2]
network_igraph <- graph_from_edgelist(as.matrix(edge_list_simple), directed = FALSE)

# --- Network metric ---
avg_path_length <- mean_distance(network_igraph, directed = FALSE)
edge_density_value <- edge_density(network_igraph, loops = FALSE)
edge_connectivity_value <- edge_connectivity(network_igraph)
components_count <- components(network_igraph)$no

# --- Centrality ---
degree_centrality <- degree(network_igraph, mode = "all", normalized = TRUE)
betweenness_centrality <- betweenness(network_igraph, normalized = TRUE)

# --- Correlations ---
cor <- net_genus$edgelist1
significant_correlations <- cor[abs(cor$asso) > 0.3, ]

# --- Stampa in console ---
cat("Average Path Length:", avg_path_length, "\n")
cat("Edge Density:", edge_density_value, "\n")
cat("Edge Connectivity:", edge_connectivity_value, "\n")
cat("Number of connected components:", components_count, "\n")

# --- Save in Excel ---
output_list <- list(
  "Metrics" = data.frame(
    Metric = c("Average Path Length", "Edge Density", "Edge Connectivity", "Connected Components"),
    Value = c(avg_path_length, edge_density_value, edge_connectivity_value, components_count)
  ),
  "Degree Centrality" = data.frame(Node = names(degree_centrality), DegreeCentrality = degree_centrality),
  "Betweenness Centrality" = data.frame(Node = names(betweenness_centrality), BetweennessCentrality = betweenness_centrality),
  "Significant Correlations" = significant_correlations
)

writexl::write_xlsx (output_list, path = "network_analysis_results_morning_T1.xlsx")

library(ggplot2)
library(reshape2)

correlation_matrix <- dcast(significant_correlations, v1 ~ v2, value.var = "asso", fun.aggregate = mean, fill = 0)
melt_correlation <- melt(correlation_matrix, id.vars = "v1")
ggplot(melt_correlation, aes(v1, variable, fill = value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), 
                       name = "Correlation") +
  theme_minimal() + 
  labs(x = "Genera", y = "Genera") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


############# GRAFICI ############# 
# Compute layout
graph3 <- igraph::graph_from_adjacency_matrix(net_genus$adjaMat1, 
                                              weighted = TRUE)
set.seed(123456)
lay_fr <- igraph::layout_with_fr(graph3)

# Row names of the layout matrix must match the node names
rownames(lay_fr) <- rownames(net_genus$adjaMat1)

plot(props_genus,
     layout = "layout_with_fr",
     shortenLabels = "none",
     labelLength = 50,
     labelPattern = c(5, "'", 3, "'", 3),
     labelScale = FALSE,
     rmSingles = TRUE,
     nodeSize = "clr",
     nodeColor = "cluster",
     hubBorderCol = "darkgray",
     cexNodes = 2,
     cexLabels = 0.5,
     cexHubLabels = 1.2,
     title1 = "Network on genus level with Pearson correlations", 
     showTitle = TRUE,
     cexTitle = 1.5)

legend(0.7, 0.9, cex = 1, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)
sort(colSums(net_genus$normCounts1), decreasing = TRUE)[1:10]

### ADD PHYLUM COLORS
taxtab <- as(tax_table(amgut_genus_renamed), "matrix")
phyla <- as.factor(gsub("p__", "", taxtab[, "Phylum"]))
names(phyla) <- taxtab[, "Genus"]
#table(phyla)

phylcol <- c("cyan", "blue3", "red", "lawngreen", "yellow", "deeppink", 'green3', 'brown', 'lightblue', 'red3')

png("Network_morning_T1.png", width = 35, height = 25, units = "cm", res = 300)

plot(props_genus,
     layout = "spring",
     repulsion = 0.84,
     shortenLabels = "none",
     charToRm = "g__",
     labelScale = FALSE,
     rmSingles = TRUE,
     nodeSize = "clr",
     nodeSizeSpread = 4,
     nodeColor = "feature", 
     featVecCol = phyla, 
     colorVec =  phylcol,
     posCol = "darkturquoise", 
     negCol = "orange",
     edgeTranspLow = 0,
     edgeTranspHigh = 40,
     cexNodes = 2,
     cexLabels = 1,
     cexHubLabels = 1.5,
     title1 = "Network analysis on Genus level with Speraman correlation\nMorning Diet (T1)", 
     showTitle = TRUE,
     cexTitle = 1.5)

# Colors used in the legend should be equally transparent as in the plot
phylcol_transp <- colToTransp(phylcol, 60)

legend(-1, 1.1, cex = 0.7, pt.cex = 2.5, title = "Phylum:", 
       legend = levels(phyla), col = phylcol_transp, bty = "n", pch = 16)

legend(0.8, 1.1, cex = 0.7, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("darkturquoise","orange"), 
       bty = "n", horiz = TRUE)

dev.off()

########## EVENING ########## 
sample_data(data)
morning_evening<-subset_samples(data, Time=='T1')
morning_evening<-subset_samples(morning_evening, Type_diet=='evening')
sample_data(morning_evening)

data.phy <- tax_glom(morning_evening, taxrank = "Genus", NArm = F)
data.phy.prop <- transform_sample_counts(data.phy, function(x) x/sum(x)*100)
top5_Genus <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:20] 
prune.dat_top5_Genus <- prune_taxa(top5_Genus, data.phy.prop)
otu_table(prune.dat_top5_Genus)

sample_data(data.phy.prop) #morning e top 15 genera
data.phy.prop<-prune.dat_top5_Genus
sample_data(data.phy.prop)
#data.phy.prop <- subset_taxa(data.phy.prop, !grepl("uncultured", tax_table(data.phy.prop)[, "Genus"])) 

###### NETWORK ANALISI ###### 
amgut_genus <- data.phy.prop
tax_table(data.phy.prop)

taxtab <- as(tax_table(amgut_genus), "matrix")

amgut_genus_renamed <- renameTaxa(amgut_genus, 
                                  pat = "<name>", 
                                  substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Genus")
tax_table(amgut_genus_renamed)

net_genus <- netConstruct(amgut_genus_renamed,
                          taxRank = "Genus",
                          measure = "spearman",
                          zeroMethod = "multRepl",
                          normMethod = "clr",
                          sparsMethod = "threshold",
                          thresh = 0.3,
                          verbose = 3)
props_genus <- netAnalyze(net_genus, clustMethod = "cluster_fast_greedy")

#CENTRALITY
centro<-props_genus$centralities$eigenv1 
centro

df<-as.data.frame(centro)

library(dplyr)
df_ordinato <- df %>%
  arrange(desc(centro))

print(df_ordinato)

props_genus$compSize1

#AVERAGE PATH LENGHT
#install.packages('igraph')
library(igraph)
names(net_genus)
                                         
# --- Network ---
edge_list <- net_genus$edgelist1  # Adatta al nome reale se diverso
edge_list_simple <- edge_list[, 1:2]
network_igraph <- graph_from_edgelist(as.matrix(edge_list_simple), directed = FALSE)

# --- Network metric ---
avg_path_length <- mean_distance(network_igraph, directed = FALSE)
edge_density_value <- edge_density(network_igraph, loops = FALSE)
edge_connectivity_value <- edge_connectivity(network_igraph)
components_count <- components(network_igraph)$no

# --- Centrality ---
degree_centrality <- degree(network_igraph, mode = "all", normalized = TRUE)
betweenness_centrality <- betweenness(network_igraph, normalized = TRUE)

# --- Correlations ---
cor <- net_genus$edgelist1
significant_correlations <- cor[abs(cor$asso) > 0.3, ]

# --- Print in console ---
cat("Average Path Length:", avg_path_length, "\n")
cat("Edge Density:", edge_density_value, "\n")
cat("Edge Connectivity:", edge_connectivity_value, "\n")
cat("Number of connected components:", components_count, "\n")

# --- Save in Excel ---
output_list <- list(
  "Metrics" = data.frame(
    Metric = c("Average Path Length", "Edge Density", "Edge Connectivity", "Connected Components"),
    Value = c(avg_path_length, edge_density_value, edge_connectivity_value, components_count)
  ),
  "Degree Centrality" = data.frame(Node = names(degree_centrality), DegreeCentrality = degree_centrality),
  "Betweenness Centrality" = data.frame(Node = names(betweenness_centrality), BetweennessCentrality = betweenness_centrality),
  "Significant Correlations" = significant_correlations
)

writexl::write_xlsx (output_list, path = "network_analysis_results_evening_T1.xlsx")


library(ggplot2)
library(reshape2)


correlation_matrix <- dcast(significant_correlations, v1 ~ v2, value.var = "asso", fun.aggregate = mean, fill = 0)
melt_correlation <- melt(correlation_matrix, id.vars = "v1")


ggplot(melt_correlation, aes(v1, variable, fill = value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), 
                       name = "Correlation") +
  theme_minimal() + 
  labs(x = "Genera", y = "Genera") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


############# GRAFICI ############# 
# Compute layout
graph3 <- igraph::graph_from_adjacency_matrix(net_genus$adjaMat1, 
                                              weighted = TRUE)
set.seed(123456)
lay_fr <- igraph::layout_with_fr(graph3)

# Row names of the layout matrix must match the node names
rownames(lay_fr) <- rownames(net_genus$adjaMat1)

plot(props_genus,
     layout = "layout_with_fr",
     shortenLabels = "none",
     labelLength = 50,
     labelPattern = c(5, "'", 3, "'", 3),
     labelScale = FALSE,
     rmSingles = TRUE,
     nodeSize = "clr",
     nodeColor = "cluster",
     hubBorderCol = "darkgray",
     cexNodes = 2,
     cexLabels = 0.5,
     cexHubLabels = 1.2,
     title1 = "Network on genus level with Pearson correlations", 
     showTitle = TRUE,
     cexTitle = 1.5)

legend(0.7, 0.9, cex = 1, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)
sort(colSums(net_genus$normCounts1), decreasing = TRUE)[1:10]

### ADD PHYLUM COLORS
taxtab <- as(tax_table(amgut_genus_renamed), "matrix")
phyla <- as.factor(gsub("p__", "", taxtab[, "Phylum"]))
names(phyla) <- taxtab[, "Genus"]
#table(phyla)

png("Network_evening_T1.png", width = 35, height = 25, units = "cm", res = 300)
phylcol <- c("cyan", "blue3", "red", "yellow", "deeppink", 'green3', 'brown', 'lightblue', 'red3')

plot(props_genus,
     layout = "spring",
     repulsion = 0.84,
     shortenLabels = "none",
     charToRm = "g__",
     labelScale = FALSE,
     rmSingles = TRUE,
     nodeSize = "clr",
     nodeSizeSpread = 4,
     nodeColor = "feature", 
     featVecCol = phyla, 
     colorVec =  phylcol,
     posCol = "darkturquoise", 
     negCol = "orange",
     edgeTranspLow = 0,
     edgeTranspHigh = 40,
     cexNodes = 2,
     cexLabels = 1,
     cexHubLabels = 1.5,
     title1 = "Network analysis on Genus level with Speraman correlation\nEvening Diet (T1)", 
     showTitle = TRUE,
     cexTitle = 1.5)

# Colors used in the legend should be equally transparent as in the plot
phylcol_transp <- colToTransp(phylcol, 60)

legend(-1, 1.1, cex = 0.7, pt.cex = 2.5, title = "Phylum:", 
       legend = levels(phyla), col = phylcol_transp, bty = "n", pch = 16)

legend(0.6, 1.1, cex = 0.7, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("darkturquoise","orange"), 
       bty = "n", horiz = TRUE)
dev.off()

########## CONTROL ########## 
sample_data(data)
morning_evening<-subset_samples(data, Time=='T1')
morning_evening<-subset_samples(morning_evening, Type_diet=='control')
sample_data(morning_evening)

data.phy <- tax_glom(morning_evening, taxrank = "Genus", NArm = F)
data.phy.prop <- transform_sample_counts(data.phy, function(x) x/sum(x)*100)
top5_Genus <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:20] 
prune.dat_top5_Genus <- prune_taxa(top5_Genus, data.phy.prop)
otu_table(prune.dat_top5_Genus)

sample_data(data.phy.prop) #morning e top 15 genera
data.phy.prop<-prune.dat_top5_Genus
sample_data(data.phy.prop)
#data.phy.prop <- subset_taxa(data.phy.prop, !grepl("uncultured", tax_table(data.phy.prop)[, "Genus"])) 

###### NETWORK ANALISI ###### 
amgut_genus <- data.phy.prop
tax_table(data.phy.prop)

taxtab <- as(tax_table(amgut_genus), "matrix")

# Rename taxonomic table and make Genus unique (rinomina duplicati es. se ci fossero due uncultured li chiamerebbe uncultured_01 e 02)
amgut_genus_renamed <- renameTaxa(amgut_genus, 
                                  pat = "<name>", 
                                  substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Genus")
tax_table(amgut_genus_renamed)
# Network construction and analysis
net_genus <- netConstruct(amgut_genus_renamed,
                          taxRank = "Genus",
                          measure = "spearman",
                          zeroMethod = "multRepl",
                          normMethod = "clr",
                          sparsMethod = "threshold",
                          thresh = 0.3,
                          verbose = 3)
props_genus <- netAnalyze(net_genus, clustMethod = "cluster_fast_greedy")

#CENTRALITY
centro<-props_genus$centralities$eigenv1 
centro

df<-as.data.frame(centro)

library(dplyr)
df_ordinato <- df %>%
  arrange(desc(centro))

print(df_ordinato)

props_genus$compSize1

#AVERAGE PATH LENGHT
#install.packages('igraph')
library(igraph)

names(net_genus)
# --- Network ---
edge_list <- net_genus$edgelist1  # Adatta al nome reale se diverso
edge_list_simple <- edge_list[, 1:2]
network_igraph <- graph_from_edgelist(as.matrix(edge_list_simple), directed = FALSE)

# --- Network metric ---
avg_path_length <- mean_distance(network_igraph, directed = FALSE)
edge_density_value <- edge_density(network_igraph, loops = FALSE)
edge_connectivity_value <- edge_connectivity(network_igraph)
components_count <- components(network_igraph)$no

# --- Centrality ---
degree_centrality <- degree(network_igraph, mode = "all", normalized = TRUE)
betweenness_centrality <- betweenness(network_igraph, normalized = TRUE)

# --- Correlations ---
cor <- net_genus$edgelist1
significant_correlations <- cor[abs(cor$asso) > 0.3, ]

# --- Print in console ---
cat("Average Path Length:", avg_path_length, "\n")
cat("Edge Density:", edge_density_value, "\n")
cat("Edge Connectivity:", edge_connectivity_value, "\n")
cat("Number of connected components:", components_count, "\n")

# --- Save in Excel ---
output_list <- list(
  "Metrics" = data.frame(
    Metric = c("Average Path Length", "Edge Density", "Edge Connectivity", "Connected Components"),
    Value = c(avg_path_length, edge_density_value, edge_connectivity_value, components_count)
  ),
  "Degree Centrality" = data.frame(Node = names(degree_centrality), DegreeCentrality = degree_centrality),
  "Betweenness Centrality" = data.frame(Node = names(betweenness_centrality), BetweennessCentrality = betweenness_centrality),
  "Significant Correlations" = significant_correlations
)

writexl::write_xlsx (output_list, path = "network_analysis_results_control_T1.xlsx")

library(ggplot2)
library(reshape2)


correlation_matrix <- dcast(significant_correlations, v1 ~ v2, value.var = "asso", fun.aggregate = mean, fill = 0)
melt_correlation <- melt(correlation_matrix, id.vars = "v1")
ggplot(melt_correlation, aes(v1, variable, fill = value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), 
                       name = "Correlation") +
  theme_minimal() + 
  labs(x = "Genera", y = "Genera") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


############# Graphs ############# 
# Compute layout
graph3 <- igraph::graph_from_adjacency_matrix(net_genus$adjaMat1, 
                                              weighted = TRUE)
set.seed(123456)
lay_fr <- igraph::layout_with_fr(graph3)

# Row names of the layout matrix must match the node names
rownames(lay_fr) <- rownames(net_genus$adjaMat1)

plot(props_genus,
     layout = "layout_with_fr",
     shortenLabels = "none",
     labelLength = 50,
     labelPattern = c(5, "'", 3, "'", 3),
     labelScale = FALSE,
     rmSingles = TRUE,
     nodeSize = "clr",
     nodeColor = "cluster",
     hubBorderCol = "darkgray",
     cexNodes = 2,
     cexLabels = 0.5,
     cexHubLabels = 1.2,
     title1 = "Network analysis on Genus level with Speraman correlation\nControl Diet (T1)", 
     showTitle = TRUE,
     cexTitle = 1.5)

legend(0.7, 0.9, cex = 1, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)
sort(colSums(net_genus$normCounts1), decreasing = TRUE)[1:10]

### ADD PHYLUM COLORS
taxtab <- as(tax_table(amgut_genus_renamed), "matrix")
phyla <- as.factor(gsub("p__", "", taxtab[, "Phylum"]))
names(phyla) <- taxtab[, "Genus"]
#table(phyla)

png("Network_control_T1.png", width = 35, height = 25, units = "cm", res = 300)
phylcol <- c("cyan", "blue3", "red", "lawngreen", "yellow", "deeppink", 'green3', 'brown', 'lightblue', 'red3')

plot(props_genus,
     layout = "spring",
     repulsion = 0.84,
     shortenLabels = "none",
     charToRm = "g__",
     labelScale = FALSE,
     rmSingles = TRUE,
     nodeSize = "clr",
     nodeSizeSpread = 4,
     nodeColor = "feature", 
     featVecCol = phyla, 
     colorVec =  phylcol,
     posCol = "darkturquoise", 
     negCol = "orange",
     edgeTranspLow = 0,
     edgeTranspHigh = 40,
     cexNodes = 2,
     cexLabels = 1,
     cexHubLabels = 1.5,
     title1 = "Network analysis on Genus level with Speraman correlation\nControl Diet (T1)", 
     showTitle = TRUE,
     cexTitle = 1.5)

# Colors used in the legend should be equally transparent as in the plot
phylcol_transp <- colToTransp(phylcol, 60)

legend(-1.1, 1.15, cex = 0.7, pt.cex = 2.5, title = "Phylum:", 
       legend = levels(phyla), col = phylcol_transp, bty = "n", pch = 16)

legend(0.7, 1.1, cex = 0.7, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("darkturquoise","orange"), 
       bty = "n", horiz = TRUE)

dev.off()

###### ON TOP 20 GENERA MORNING ###### 
sample_data(data)
morning_evening<-subset_samples(data, Time=='T0')
morning_evening<-subset_samples(morning_evening, Type_diet=='morning')
sample_data(morning_evening)

data.phy <- tax_glom(morning_evening, taxrank = "Genus", NArm = F)
data.phy.prop <- transform_sample_counts(data.phy, function(x) x/sum(x)*100)
top5_Genus <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:20] 
prune.dat_top5_Genus <- prune_taxa(top5_Genus, data.phy.prop)
otu_table(prune.dat_top5_Genus)

sample_data(data.phy.prop) #morning e top 15 genera
data.phy.prop<-prune.dat_top5_Genus
sample_data(data.phy.prop)
#data.phy.prop <- subset_taxa(data.phy.prop, !grepl("uncultured", tax_table(data.phy.prop)[, "Genus"])) 

###### NETWORK ANALYSIS ###### 
amgut_genus <- data.phy.prop
tax_table(data.phy.prop)

taxtab <- as(tax_table(amgut_genus), "matrix")

amgut_genus_renamed <- renameTaxa(amgut_genus, 
                                  pat = "<name>", 
                                  substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Genus")
tax_table(amgut_genus_renamed)
                                         
# Network construction and analysis
net_genus <- netConstruct(amgut_genus_renamed,
                          taxRank = "Genus",
                          measure = "spearman",
                          zeroMethod = "multRepl",
                          normMethod = "clr",
                          sparsMethod = "threshold",
                          thresh = 0.3,
                          verbose = 3)
props_genus <- netAnalyze(net_genus, clustMethod = "cluster_fast_greedy")

#CENTRALITY
centro<-props_genus$centralities$eigenv1 
centro
df<-as.data.frame(centro)

library(dplyr)
df_ordinato <- df %>%
  arrange(desc(centro))

print(df_ordinato) 

props_genus$compSize1

#AVERAGE PATH LENGHT
#install.packages('igraph')
library(igraph)

names(net_genus)

# --- Network ---
edge_list <- net_genus$edgelist1  # Adatta al nome reale se diverso
edge_list_simple <- edge_list[, 1:2]
network_igraph <- graph_from_edgelist(as.matrix(edge_list_simple), directed = FALSE)

# --- Network metric ---
avg_path_length <- mean_distance(network_igraph, directed = FALSE)
edge_density_value <- edge_density(network_igraph, loops = FALSE)
edge_connectivity_value <- edge_connectivity(network_igraph)
components_count <- components(network_igraph)$no

# --- Centrality ---
degree_centrality <- degree(network_igraph, mode = "all", normalized = TRUE)
betweenness_centrality <- betweenness(network_igraph, normalized = TRUE)

# --- Correlations ---
cor <- net_genus$edgelist1
significant_correlations <- cor[abs(cor$asso) > 0.3, ]

# --- Print in console ---
cat("Average Path Length:", avg_path_length, "\n")
cat("Edge Density:", edge_density_value, "\n")
cat("Edge Connectivity:", edge_connectivity_value, "\n")
cat("Number of connected components:", components_count, "\n")

# --- Save in Excel ---
output_list <- list(
  "Metrics" = data.frame(
    Metric = c("Average Path Length", "Edge Density", "Edge Connectivity", "Connected Components"),
    Value = c(avg_path_length, edge_density_value, edge_connectivity_value, components_count)
  ),
  "Degree Centrality" = data.frame(Node = names(degree_centrality), DegreeCentrality = degree_centrality),
  "Betweenness Centrality" = data.frame(Node = names(betweenness_centrality), BetweennessCentrality = betweenness_centrality),
  "Significant Correlations" = significant_correlations
)

writexl::write_xlsx (output_list, path = "network_analysis_results_morning_T0.xlsx")

library(ggplot2)
library(reshape2)


correlation_matrix <- dcast(significant_correlations, v1 ~ v2, value.var = "asso", fun.aggregate = mean, fill = 0)
melt_correlation <- melt(correlation_matrix, id.vars = "v1")
ggplot(melt_correlation, aes(v1, variable, fill = value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), 
                       name = "Correlation") +
  theme_minimal() + 
  labs(x = "Genera", y = "Genera") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


############# GRAPHS ############# 
# Compute layout
graph3 <- igraph::graph_from_adjacency_matrix(net_genus$adjaMat1, 
                                              weighted = TRUE)
set.seed(123456)
lay_fr <- igraph::layout_with_fr(graph3)

# Row names of the layout matrix must match the node names
rownames(lay_fr) <- rownames(net_genus$adjaMat1)

plot(props_genus,
     layout = "layout_with_fr",
     shortenLabels = "none",
     labelLength = 50,
     labelPattern = c(5, "'", 3, "'", 3),
     labelScale = FALSE,
     rmSingles = TRUE,
     nodeSize = "clr",
     nodeColor = "cluster",
     hubBorderCol = "darkgray",
     cexNodes = 2,
     cexLabels = 0.5,
     cexHubLabels = 1.2,
     title1 = "Network on genus level with Pearson correlations", 
     showTitle = TRUE,
     cexTitle = 1.5)

legend(0.7, 0.9, cex = 1, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)
sort(colSums(net_genus$normCounts1), decreasing = TRUE)[1:10]

### ADD PHYLUM COLORS
taxtab <- as(tax_table(amgut_genus_renamed), "matrix")
phyla <- as.factor(gsub("p__", "", taxtab[, "Phylum"]))
names(phyla) <- taxtab[, "Genus"]
#table(phyla)

png("Network_morning_T0.png", width = 35, height = 25, units = "cm", res = 300)
phylcol <- c("cyan", "blue3", "red", "lawngreen", "yellow", "deeppink", 'green3', 'brown', 'lightblue', 'red3')

plot(props_genus,
     layout = "spring",
     repulsion = 0.84,
     shortenLabels = "none",
     charToRm = "g__",
     labelScale = FALSE,
     rmSingles = TRUE,
     nodeSize = "clr",
     nodeSizeSpread = 4,
     nodeColor = "feature", 
     featVecCol = phyla, 
     colorVec =  phylcol,
     posCol = "darkturquoise", 
     negCol = "orange",
     edgeTranspLow = 0,
     edgeTranspHigh = 40,
     cexNodes = 2,
     cexLabels = 1,
     cexHubLabels = 1.5,
     title1 = "Network analysis on Genus level with Speraman correlation\nMorning Diet (T0)", 
     showTitle = TRUE,
     cexTitle = 1.5)

# Colors used in the legend should be equally transparent as in the plot
phylcol_transp <- colToTransp(phylcol, 60)

legend(-1, 1.1, cex = 0.7, pt.cex = 2.5, title = "Phylum:", 
       legend = levels(phyla), col = phylcol_transp, bty = "n", pch = 16)

legend(0.7, 1.1, cex = 0.7, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("darkturquoise","orange"), 
       bty = "n", horiz = TRUE)
dev.off()

########## EVENING ########## 
sample_data(data)
morning_evening<-subset_samples(data, Time=='T0')
morning_evening<-subset_samples(morning_evening, Type_diet=='evening')
sample_data(morning_evening)

data.phy <- tax_glom(morning_evening, taxrank = "Genus", NArm = F)
data.phy.prop <- transform_sample_counts(data.phy, function(x) x/sum(x)*100)
top5_Genus <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:20] 
prune.dat_top5_Genus <- prune_taxa(top5_Genus, data.phy.prop)
otu_table(prune.dat_top5_Genus)

sample_data(data.phy.prop) #morning e top 15 genera
data.phy.prop<-prune.dat_top5_Genus
sample_data(data.phy.prop)
#data.phy.prop <- subset_taxa(data.phy.prop, !grepl("uncultured", tax_table(data.phy.prop)[, "Genus"])) 

###### NETWORK ANALYSIS ###### 
amgut_genus <- data.phy.prop
tax_table(data.phy.prop)

taxtab <- as(tax_table(amgut_genus), "matrix")

# Rename taxonomic table and make Genus unique (rinomina duplicati es. se ci fossero due uncultured li chiamerebbe uncultured_01 e 02)
amgut_genus_renamed <- renameTaxa(amgut_genus, 
                                  pat = "<name>", 
                                  substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Genus")
tax_table(amgut_genus_renamed)
# Network construction and analysis
net_genus <- netConstruct(amgut_genus_renamed,
                          taxRank = "Genus",
                          measure = "spearman",
                          zeroMethod = "multRepl",
                          normMethod = "clr",
                          sparsMethod = "threshold",
                          thresh = 0.3,
                          verbose = 3)
props_genus <- netAnalyze(net_genus, clustMethod = "cluster_fast_greedy")

#CENTRALITY
centro<-props_genus$centralities$eigenv1 
centro
df<-as.data.frame(centro)

library(dplyr)
df_ordinato <- df %>%
  arrange(desc(centro))

print(df_ordinato) 

props_genus$compSize1

#AVERAGE PATH LENGHT
#install.packages('igraph')
library(igraph)

names(net_genus)
# --- Network ---
edge_list <- net_genus$edgelist1  # Adatta al nome reale se diverso
edge_list_simple <- edge_list[, 1:2]
network_igraph <- graph_from_edgelist(as.matrix(edge_list_simple), directed = FALSE)

# --- Network metric ---
avg_path_length <- mean_distance(network_igraph, directed = FALSE)
edge_density_value <- edge_density(network_igraph, loops = FALSE)
edge_connectivity_value <- edge_connectivity(network_igraph)
components_count <- components(network_igraph)$no

# --- Centrality ---
degree_centrality <- degree(network_igraph, mode = "all", normalized = TRUE)
betweenness_centrality <- betweenness(network_igraph, normalized = TRUE)

# --- Correlations ---
cor <- net_genus$edgelist1
significant_correlations <- cor[abs(cor$asso) > 0.3, ]

# --- Print in console ---
cat("Average Path Length:", avg_path_length, "\n")
cat("Edge Density:", edge_density_value, "\n")
cat("Edge Connectivity:", edge_connectivity_value, "\n")
cat("Number of connected components:", components_count, "\n")

# --- Save in Excel ---
output_list <- list(
  "Metrics" = data.frame(
    Metric = c("Average Path Length", "Edge Density", "Edge Connectivity", "Connected Components"),
    Value = c(avg_path_length, edge_density_value, edge_connectivity_value, components_count)
  ),
  "Degree Centrality" = data.frame(Node = names(degree_centrality), DegreeCentrality = degree_centrality),
  "Betweenness Centrality" = data.frame(Node = names(betweenness_centrality), BetweennessCentrality = betweenness_centrality),
  "Significant Correlations" = significant_correlations
)

writexl::write_xlsx (output_list, path = "network_analysis_results_evening_T0.xlsx")


library(ggplot2)
library(reshape2)


correlation_matrix <- dcast(significant_correlations, v1 ~ v2, value.var = "asso", fun.aggregate = mean, fill = 0)
melt_correlation <- melt(correlation_matrix, id.vars = "v1")

ggplot(melt_correlation, aes(v1, variable, fill = value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), 
                       name = "Correlation") +
  theme_minimal() + 
  labs(x = "Genera", y = "Genera") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


############# GRAPHS ############# 
# Compute layout
graph3 <- igraph::graph_from_adjacency_matrix(net_genus$adjaMat1, 
                                              weighted = TRUE)
set.seed(123456)
lay_fr <- igraph::layout_with_fr(graph3)

# Row names of the layout matrix must match the node names
rownames(lay_fr) <- rownames(net_genus$adjaMat1)


plot(props_genus,
     layout = "layout_with_fr",
     shortenLabels = "none",
     labelLength = 50,
     labelPattern = c(5, "'", 3, "'", 3),
     labelScale = FALSE,
     rmSingles = TRUE,
     nodeSize = "clr",
     nodeColor = "cluster",
     hubBorderCol = "darkgray",
     cexNodes = 2,
     cexLabels = 0.5,
     cexHubLabels = 1.2,
     title1 = "Network on genus level with Pearson correlations", 
     showTitle = TRUE,
     cexTitle = 1.5)

legend(0.7, 0.9, cex = 1, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)
sort(colSums(net_genus$normCounts1), decreasing = TRUE)[1:10]

### ADD PHYLUM COLORS
taxtab <- as(tax_table(amgut_genus_renamed), "matrix")
phyla <- as.factor(gsub("p__", "", taxtab[, "Phylum"]))
names(phyla) <- taxtab[, "Genus"]
#table(phyla)

png("Network_evening_T0.png", width = 35, height = 25, units = "cm", res = 300)
phylcol <- c("cyan", "blue3", "red", "lawngreen", "deeppink", 'green3', 'brown', 'lightblue', 'red3')

plot(props_genus,
     layout = "spring",
     repulsion = 0.84,
     shortenLabels = "none",
     charToRm = "g__",
     labelScale = FALSE,
     rmSingles = TRUE,
     nodeSize = "clr",
     nodeSizeSpread = 4,
     nodeColor = "feature", 
     featVecCol = phyla, 
     colorVec =  phylcol,
     posCol = "darkturquoise", 
     negCol = "orange",
     edgeTranspLow = 0,
     edgeTranspHigh = 40,
     cexNodes = 2,
     cexLabels = 1,
     cexHubLabels = 1.5,
     title1 = "Network analysis on Genus level with Speraman correlation\nEvening Diet (T0)", 
     showTitle = TRUE,
     cexTitle = 1.5)

# Colors used in the legend should be equally transparent as in the plot
phylcol_transp <- colToTransp(phylcol, 60)

legend(-1, 1.1, cex = 0.7, pt.cex = 2.5, title = "Phylum:", 
       legend = levels(phyla), col = phylcol_transp, bty = "n", pch = 16)

legend(0.6, 1.1, cex = 0.7, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("darkturquoise","orange"), 
       bty = "n", horiz = TRUE)
dev.off()

########## CONTROL ########## 
sample_data(data)
morning_evening<-subset_samples(data, Time=='T0')
morning_evening<-subset_samples(morning_evening, Type_diet=='control')
sample_data(morning_evening)

data.phy <- tax_glom(morning_evening, taxrank = "Genus", NArm = F)
data.phy.prop <- transform_sample_counts(data.phy, function(x) x/sum(x)*100)
top5_Genus <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:20] 
prune.dat_top5_Genus <- prune_taxa(top5_Genus, data.phy.prop)
otu_table(prune.dat_top5_Genus)

sample_data(data.phy.prop) #morning e top 15 genera
data.phy.prop<-prune.dat_top5_Genus
sample_data(data.phy.prop)
#data.phy.prop <- subset_taxa(data.phy.prop, !grepl("uncultured", tax_table(data.phy.prop)[, "Genus"])) 

###### NETWORK ANALISI ###### 
amgut_genus <- data.phy.prop
tax_table(data.phy.prop)

taxtab <- as(tax_table(amgut_genus), "matrix")

# Rename taxonomic table and make Genus unique (rinomina duplicati es. se ci fossero due uncultured li chiamerebbe uncultured_01 e 02)
amgut_genus_renamed <- renameTaxa(amgut_genus, 
                                  pat = "<name>", 
                                  substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Genus")
tax_table(amgut_genus_renamed)
# Network construction and analysis
net_genus <- netConstruct(amgut_genus_renamed,
                          taxRank = "Genus",
                          measure = "spearman",
                          zeroMethod = "multRepl",
                          normMethod = "clr",
                          sparsMethod = "threshold",
                          thresh = 0.3,
                          verbose = 3)
props_genus <- netAnalyze(net_genus, clustMethod = "cluster_fast_greedy")

#CENTRALITY
centro<-props_genus$centralities$eigenv1 
centro

df<-as.data.frame(centro)

library(dplyr)
df_ordinato <- df %>%
  arrange(desc(centro))

print(df_ordinato) 

props_genus$compSize1

#AVERAGE PATH LENGHT
#install.packages('igraph')
library(igraph)

names(net_genus)
# --- Network ---
edge_list <- net_genus$edgelist1  # Adatta al nome reale se diverso
edge_list_simple <- edge_list[, 1:2]
network_igraph <- graph_from_edgelist(as.matrix(edge_list_simple), directed = FALSE)

# --- Network metric ---
avg_path_length <- mean_distance(network_igraph, directed = FALSE)
edge_density_value <- edge_density(network_igraph, loops = FALSE)
edge_connectivity_value <- edge_connectivity(network_igraph)
components_count <- components(network_igraph)$no

# --- Centrality ---
degree_centrality <- degree(network_igraph, mode = "all", normalized = TRUE)
betweenness_centrality <- betweenness(network_igraph, normalized = TRUE)

# --- Correlations ---
cor <- net_genus$edgelist1
significant_correlations <- cor[abs(cor$asso) > 0.3, ]

# --- Print in console ---
cat("Average Path Length:", avg_path_length, "\n")
cat("Edge Density:", edge_density_value, "\n")
cat("Edge Connectivity:", edge_connectivity_value, "\n")
cat("Number of connected components:", components_count, "\n")

# --- Save in Excel ---
output_list <- list(
  "Metrics" = data.frame(
    Metric = c("Average Path Length", "Edge Density", "Edge Connectivity", "Connected Components"),
    Value = c(avg_path_length, edge_density_value, edge_connectivity_value, components_count)
  ),
  "Degree Centrality" = data.frame(Node = names(degree_centrality), DegreeCentrality = degree_centrality),
  "Betweenness Centrality" = data.frame(Node = names(betweenness_centrality), BetweennessCentrality = betweenness_centrality),
  "Significant Correlations" = significant_correlations
)

writexl::write_xlsx (output_list, path = "network_analysis_results_control_T0.xlsx")

library(ggplot2)
library(reshape2)

correlation_matrix <- dcast(significant_correlations, v1 ~ v2, value.var = "asso", fun.aggregate = mean, fill = 0)
melt_correlation <- melt(correlation_matrix, id.vars = "v1")

ggplot(melt_correlation, aes(v1, variable, fill = value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), 
                       name = "Correlation") +
  theme_minimal() + 
  labs(x = "Genera", y = "Genera") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


############# GRAPHS ############# 
# Compute layout
graph3 <- igraph::graph_from_adjacency_matrix(net_genus$adjaMat1, 
                                              weighted = TRUE)
set.seed(123456)
lay_fr <- igraph::layout_with_fr(graph3)

# Row names of the layout matrix must match the node names
rownames(lay_fr) <- rownames(net_genus$adjaMat1)

plot(props_genus,
     layout = "layout_with_fr",
     shortenLabels = "none",
     labelLength = 50,
     labelPattern = c(5, "'", 3, "'", 3),
     labelScale = FALSE,
     rmSingles = TRUE,
     nodeSize = "clr",
     nodeColor = "cluster",
     hubBorderCol = "darkgray",
     cexNodes = 2,
     cexLabels = 0.5,
     cexHubLabels = 1.2,
     title1 = "Network analysis on Genus level with Speraman correlation\nControl Diet (T1)", 
     showTitle = TRUE,
     cexTitle = 1.5)

legend(0.7, 0.9, cex = 1, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)
sort(colSums(net_genus$normCounts1), decreasing = TRUE)[1:10]

### ADD PHYLUM COLORS
taxtab <- as(tax_table(amgut_genus_renamed), "matrix")
phyla <- as.factor(gsub("p__", "", taxtab[, "Phylum"]))
names(phyla) <- taxtab[, "Genus"]
#table(phyla)

png("Network_control_T0.png", width = 35, height = 25, units = "cm", res = 300)
phylcol <- c("cyan", "blue3", "red", "lawngreen", "yellow", "deeppink", 'green3', 'brown', 'lightblue', 'red3')

plot(props_genus,
     layout = "spring",
     repulsion = 0.84,
     shortenLabels = "none",
     charToRm = "g__",
     labelScale = FALSE,
     rmSingles = TRUE,
     nodeSize = "clr",
     nodeSizeSpread = 4,
     nodeColor = "feature", 
     featVecCol = phyla, 
     colorVec =  phylcol,
     posCol = "darkturquoise", 
     negCol = "orange",
     edgeTranspLow = 0,
     edgeTranspHigh = 40,
     cexNodes = 2,
     cexLabels = 1,
     cexHubLabels = 1.5,
     title1 = "Network analysis on Genus level with Speraman correlation\nControl Diet (T0)", 
     showTitle = TRUE,
     cexTitle = 1.5)

# Colors used in the legend should be equally transparent as in the plot
phylcol_transp <- colToTransp(phylcol, 60)

legend(-1, 1.1, cex = 0.7, pt.cex = 2.5, title = "Phylum:", 
       legend = levels(phyla), col = phylcol_transp, bty = "n", pch = 16)

legend(0.7, 1.1, cex = 0.7, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("darkturquoise","orange"), 
       bty = "n", horiz = TRUE)
dev.off()

