########## LONGITUDINAL MORNING vs EVENING: ALPHA, BETA, DESeq2 ########## 

library(phyloseq)
library(vegan)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(permute)

set.seed(1994)
outdir <- "Results_Cronotipo/Longitudinal_Morning_Evening"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "Alpha"), showWarnings = FALSE)
dir.create(file.path(outdir, "Beta"), showWarnings = FALSE)
dir.create(file.path(outdir, "DESeq2"), showWarnings = FALSE)

COLORS <- c(evening = "#0047AB", morning = "#FFB347")
GROUP_LEVELS <- c("evening", "morning")
TIME_LEVELS <- c("T0", "T1")

#### 1. Select only morning/evening with T0 and T1 ####

meta0 <- data.frame(sample_data(data))
required_cols <- c("Sample_Name", "Time", "Type_diet")
if (!all(required_cols %in% colnames(meta0))) stop("No metadata: Sample_Name, Time o Type_diet.")

meta0$Sample_Name <- as.character(meta0$Sample_Name)
meta0$Time <- as.character(meta0$Time)
meta0$Type_diet <- tolower(as.character(meta0$Type_diet))

meta0 <- meta0[meta0$Type_diet %in% GROUP_LEVELS & meta0$Time %in% TIME_LEVELS, , drop = FALSE]
if (nrow(meta0) == 0) stop("No samples morning/evening with Time T0/T1. Control metadata.")

groups_per_subject <- tapply(meta0$Type_diet, meta0$Sample_Name, function(x) length(unique(x)))
if (any(groups_per_subject > 1)) stop("Problem")

tab <- table(meta0$Sample_Name, factor(meta0$Time, levels = TIME_LEVELS))
paired_ids <- rownames(tab)[tab[, "T0"] == 1 & tab[, "T1"] == 1]
meta0 <- meta0[meta0$Sample_Name %in% paired_ids, , drop = FALSE]
if (length(paired_ids) < 2) stop("Few samples")

ps <- prune_samples(rownames(meta0), data)
sample_data(ps)$Sample_Name <- factor(as.character(sample_data(ps)$Sample_Name))
sample_data(ps)$Time <- factor(as.character(sample_data(ps)$Time), levels = TIME_LEVELS)
sample_data(ps)$Type_diet <- factor(tolower(as.character(sample_data(ps)$Type_diet)), levels = GROUP_LEVELS)

meta <- data.frame(sample_data(ps))
write.csv(as.data.frame(table(meta$Type_diet, meta$Time)), file.path(outdir, "Sample_counts_by_group_time.csv"), row.names = FALSE)

cat("\nSamples for group and time:\n")
print(table(sample_data(ps)$Type_diet, sample_data(ps)$Time))
cat("\nNumber of samples:\n")
print(table(tapply(as.character(sample_data(ps)$Type_diet), sample_data(ps)$Sample_Name, unique)))

##### 2. ALPHA DIVERSITY #####
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
alpha <- estimate_richness(ps_genus, measures = c("Observed", "Shannon"))
alpha$SampleID <- rownames(alpha)
md <- data.frame(sample_data(ps_genus))
md$SampleID <- rownames(md)
alpha <- left_join(alpha, md, by = "SampleID")
alpha$Observed <- as.numeric(alpha$Observed)
alpha$Shannon <- as.numeric(alpha$Shannon)
alpha$Evenness <- ifelse(alpha$Observed > 1, alpha$Shannon / log(alpha$Observed), NA_real_)
alpha$Time <- factor(as.character(alpha$Time), levels = TIME_LEVELS)
alpha$Type_diet <- factor(as.character(alpha$Type_diet), levels = GROUP_LEVELS)
alpha$Sample_Name <- factor(alpha$Sample_Name)

alpha_metrics <- c("Observed", "Shannon", "Evenness")
alpha_models <- list()
alpha_interaction_results <- list()
alpha_change_results <- list()

for (metric in alpha_metrics) {
  dat <- alpha[!is.na(alpha[[metric]]) & is.finite(alpha[[metric]]), , drop = FALSE]
  form_full <- as.formula(paste0(metric, " ~ Type_diet * Time + (1 | Sample_Name)"))
  form_reduced <- as.formula(paste0(metric, " ~ Type_diet + Time + (1 | Sample_Name)"))
  model_full <- lmer(form_full, data = dat, REML = FALSE)
  model_reduced <- lmer(form_reduced, data = dat, REML = FALSE)
  lrt <- anova(model_reduced, model_full)
  p_int <- lrt$`Pr(>Chisq)`[2]
  alpha_models[[metric]] <- model_full
  alpha_interaction_results[[metric]] <- data.frame(Metric = metric, Chi_square = lrt$Chisq[2], df = lrt$Df[2], p_value = lrt$`Pr(>Chisq)`[2])

  emm_time <- emmeans(model_full, ~ Time | Type_diet)
  changes <- as.data.frame(contrast(emm_time, method = "revpairwise", adjust = "none"))
  changes$Metric <- metric
  alpha_change_results[[metric]] <- changes

  p <- ggplot(dat, aes(x = Time, y = .data[[metric]])) +
    geom_line(aes(group = Sample_Name), color = "grey65", alpha = 0.55, linewidth = 0.45) +
    geom_boxplot(aes(group = interaction(Type_diet, Time), fill = Type_diet), width = 0.52, alpha = 0.22, outlier.shape = NA, position = position_dodge(width = 0.16)) +
    geom_point(aes(color = Type_diet), position = position_jitter(width = 0.035, height = 0, seed = 1994), size = 2.1, alpha = 0.85) +
    stat_summary(aes(group = Type_diet, color = Type_diet), fun = mean, geom = "line", linewidth = 1.05, position = position_dodge(width = 0.05)) +
    stat_summary(aes(group = Type_diet, color = Type_diet), fun = mean, geom = "point", size = 3.1, shape = 18) +
    scale_color_manual(values = COLORS, breaks = GROUP_LEVELS) +
    scale_fill_manual(values = COLORS, breaks = GROUP_LEVELS) +
    labs(title = paste0(metric, " diversity at genus level"), subtitle = paste0("Group × time interaction: p = ", signif(p_int, 3)), x = NULL, y = metric, color = NULL, fill = NULL) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "plain", size = 14), plot.subtitle = element_text(size = 10), legend.position = "top", axis.title.y = element_text(margin = margin(r = 8)))
  print(p)
  ggsave(file.path(outdir, "Alpha", paste0("Alpha_", metric, ".png")), p, width = 7.2, height = 5.2, dpi = 320)
}
p

alpha_interaction_df <- bind_rows(alpha_interaction_results)
alpha_interaction_df$p_adj_BH <- p.adjust(alpha_interaction_df$p_value, method = "BH")
write.csv(alpha_interaction_df, file.path(outdir, "Alpha", "Alpha_interaction_tests.csv"), row.names = FALSE)
write.csv(bind_rows(alpha_change_results), file.path(outdir, "Alpha", "Alpha_T1_vs_T0_within_group.csv"), row.names = FALSE)

#### 3. BETA DIVERSITY ####

ps_beta <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
ps_rel <- transform_sample_counts(ps_beta, function(x) if (sum(x) > 0) x / sum(x) else x)
ps_sqrt <- transform_sample_counts(ps_rel, sqrt)
dist_bc <- phyloseq::distance(ps_sqrt, method = "bray")
meta_beta <- data.frame(sample_data(ps_sqrt))
meta_beta <- meta_beta[sample_names(ps_sqrt), , drop = FALSE]
meta_beta$Sample_Name <- factor(meta_beta$Sample_Name)
meta_beta$Time <- factor(as.character(meta_beta$Time), levels = TIME_LEVELS)
meta_beta$Type_diet <- factor(as.character(meta_beta$Type_diet), levels = GROUP_LEVELS)

perm_control <- permute::how(nperm = 9999, blocks = meta_beta$Sample_Name)
perm_beta <- adonis2(dist_bc ~ Type_diet * Time, data = meta_beta, permutations = perm_control, by = "margin")
print(perm_beta)
write.csv(as.data.frame(perm_beta), file.path(outdir, "Beta", "Beta_PERMANOVA.csv"), row.names = TRUE)

interaction_rows <- grep(":", rownames(perm_beta))
if (length(interaction_rows) == 0) interaction_rows <- grep("Time:Type_diet|Type_diet:Time", rownames(perm_beta))
p_beta_int <- if (length(interaction_rows) > 0) perm_beta$`Pr(>F)`[interaction_rows[1]] else NA_real_

meta_beta$Group_Time <- interaction(meta_beta$Type_diet, meta_beta$Time, drop = TRUE)
bd <- betadisper(dist_bc, group = meta_beta$Group_Time)
disp_test <- permutest(bd, permutations = perm_control)
print(disp_test)
write.csv(as.data.frame(disp_test$tab), file.path(outdir, "Beta", "Beta_dispersion_test.csv"), row.names = TRUE)

ord <- ordinate(ps_sqrt, method = "PCoA", distance = dist_bc)
ord_df <- as.data.frame(ord$vectors[, 1:2, drop = FALSE])
ord_df$SampleID <- rownames(ord_df)
ord_df <- left_join(ord_df, tibble::rownames_to_column(meta_beta, "SampleID"), by = "SampleID")
eig <- ord$values$Eigenvalues
if (is.null(eig)) eig <- rep(NA_real_, 2)
positive_eig <- eig[eig > 0]
var1 <- if (length(positive_eig) >= 1 && sum(positive_eig) > 0) round(100 * positive_eig[1] / sum(positive_eig), 1) else NA
var2 <- if (length(positive_eig) >= 2 && sum(positive_eig) > 0) round(100 * positive_eig[2] / sum(positive_eig), 1) else NA

p_beta <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2)) +
  geom_line(aes(group = Sample_Name), color = "grey65", alpha = 0.5, linewidth = 0.45) +
  geom_point(aes(color = Type_diet, shape = Time), size = 3, alpha = 0.9) +
  scale_color_manual(values = COLORS, breaks = GROUP_LEVELS) +
  scale_shape_manual(values = c(T0 = 16, T1 = 17)) +
  labs(title = "Genus-level beta diversity", subtitle = paste0("Bray–Curtis PCoA | Group × time PERMANOVA p = ", signif(p_beta_int, 3)), x = paste0("PCoA1", if (!is.na(var1)) paste0(" (", var1, "%)") else ""), y = paste0("PCoA2", if (!is.na(var2)) paste0(" (", var2, "%)") else ""), color = NULL, shape = "Time") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "plain", size = 14), plot.subtitle = element_text(size = 10), legend.position = "top")
print(p_beta)
ggsave(file.path(outdir, "Beta", "Beta_PCoA_longitudinal.png"), p_beta, width = 7.5, height = 6, dpi = 320)


#### 4. DESeq2 #####

ps_dds <- ps
sample_data(ps_dds)$Time <- factor(as.character(sample_data(ps_dds)$Time), levels = TIME_LEVELS)
sample_data(ps_dds)$Type_diet <- factor(as.character(sample_data(ps_dds)$Type_diet), levels = GROUP_LEVELS)
sample_data(ps_dds)$Sample_Name <- factor(as.character(sample_data(ps_dds)$Sample_Name))
sample_data(ps_dds)$Interaction <- factor(ifelse(sample_data(ps_dds)$Time == "T1" & sample_data(ps_dds)$Type_diet == "morning", "morning_T1", "other"), levels = c("other", "morning_T1"))

meta_dds <- data.frame(sample_data(ps_dds))
X <- model.matrix(~ Sample_Name + Time + Interaction, data = meta_dds)
if (qr(X)$rank != ncol(X)) stop("DESeq2 is not è full rank: control")

tax_ranks <- intersect(c("Phylum", "Class", "Order", "Family", "Genus"), colnames(tax_table(ps_dds)))
if (length(tax_ranks) == 0) stop("Problem")

deseq_results <- list()
for (rank in tax_ranks) {
  message("DESeq2 - livello tassonomico: ", rank)
  ps_rank <- tax_glom(ps_dds, taxrank = rank, NArm = FALSE)
  dds <- phyloseq_to_deseq2(ps_rank, ~ Sample_Name + Time + Interaction)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  if (nrow(dds) == 0) next
  dds <- DESeq(dds, test = "Wald", quiet = TRUE)
  coef_interaction <- grep("^Interaction_", resultsNames(dds), value = TRUE)
  if (length(coef_interaction) != 1) {
    warning("Coefficient Interaction ", rank, ". resultsNames: ", paste(resultsNames(dds), collapse = ", "))
    next
  }
  res <- as.data.frame(results(dds, name = coef_interaction, alpha = 0.05))
  res$Taxon_ID <- rownames(res)
  res$Taxonomic_rank <- rank
  tax_df <- as.data.frame(tax_table(ps_rank))
  tax_df$Taxon_ID <- rownames(tax_df)
  res$Taxon_name <- as.character(tax_df[[rank]][match(res$Taxon_ID, tax_df$Taxon_ID)])
  missing_name <- is.na(res$Taxon_name) | res$Taxon_name == ""
  res$Taxon_name[missing_name] <- paste0("Unclassified_", rank, "_", res$Taxon_ID[missing_name])
  deseq_results[[rank]] <- res
}

if (length(deseq_results) > 0) {
  deseq_all <- bind_rows(deseq_results)
  deseq_all$padj_global_BH <- p.adjust(deseq_all$pvalue, method = "BH")
  deseq_all <- deseq_all %>% arrange(padj_global_BH, pvalue)
  write.csv(deseq_all, file.path(outdir, "DESeq2", "DESeq2_all_interaction_results.csv"), row.names = FALSE)
  sig <- deseq_all %>% filter(!is.na(padj_global_BH), padj_global_BH < 0.05)
  write.csv(sig, file.path(outdir, "DESeq2", "DESeq2_significant_interaction_results.csv"), row.names = FALSE)
  } 

################### DESEQ2: GRAPHS ###################
COLORS <- c("evening" = "#0047AB", "morning" = "#FFB347")
GROUP_LEVELS <- c("evening", "morning")
TIME_LEVELS <- c("T0", "T1")

plot_dir <- file.path(outdir, "DESeq2", "Individual_Genus_Plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
sig_genus <- sig %>%
  filter(Taxonomic_rank == "Genus",
         !is.na(padj_global_BH),
         padj_global_BH < 0.05)

if (nrow(sig_genus) == 0) {
  message("Nessun genere significativo dopo correzione globale BH.")
} else {
  ps_genus <- tax_glom(ps_dds, taxrank = "Genus", NArm = FALSE)
  ps_genus_rel <- transform_sample_counts(
    ps_genus,
    function(x) if (sum(x) > 0) x / sum(x) * 100 else x
  )
  taxa_keep <- intersect(taxa_names(ps_genus_rel), sig_genus$Taxon_ID)
  ps_genus_rel <- prune_taxa(taxa_keep, ps_genus_rel)
  plot_data <- psmelt(ps_genus_rel) %>%
    mutate(
      Time = factor(as.character(Time), levels = TIME_LEVELS),
      Type_diet = factor(as.character(Type_diet), levels = GROUP_LEVELS),
      Sample_Name = as.character(Sample_Name),
      Time_num = as.numeric(Time),
      x_plot = Time_num + ifelse(Type_diet == "evening", -0.09, 0.09)
    ) %>%
    filter(!is.na(Time), !is.na(Type_diet))
  for (taxon_id in taxa_keep) {
    
    result_taxon <- sig_genus %>%
      filter(Taxon_ID == taxon_id)
    
    taxon_name <- result_taxon$Taxon_name[1]
    fdr <- result_taxon$padj_global_BH[1]
    log2fc <- result_taxon$log2FoldChange[1]
    
    dat_taxon <- plot_data %>%
      filter(OTU == taxon_id)
    
    if (nrow(dat_taxon) == 0) next
    
    p <- ggplot(dat_taxon, aes(x = x_plot, y = Abundance)) +
      geom_line(
        aes(group = Sample_Name),
        color = "grey65",
        alpha = 0.55,
        linewidth = 0.45
      ) +
      geom_boxplot(
        aes(group = interaction(Type_diet, Time), fill = Type_diet),
        width = 0.13,
        alpha = 0.22,
        outlier.shape = NA,
        linewidth = 0.45
      ) +
      geom_point(
        aes(color = Type_diet),
        position = position_jitter(width = 0.025, height = 0, seed = 1994),
        size = 2.1,
        alpha = 0.85
      ) +
      stat_summary(
        aes(group = Type_diet, color = Type_diet),
        fun = mean,
        geom = "line",
        linewidth = 1.05
      ) +
      
      stat_summary(
        aes(group = Type_diet, color = Type_diet),
        fun = mean,
        geom = "point",
        size = 3.1,
        shape = 18
      ) +
      
      scale_color_manual(values = COLORS, breaks = GROUP_LEVELS) +
      scale_fill_manual(values = COLORS, breaks = GROUP_LEVELS) +
      
      scale_x_continuous(
        breaks = c(1, 2),
        labels = TIME_LEVELS,
        expand = expansion(mult = c(0.12, 0.12))
      ) +
      
      scale_y_continuous(
        expand = expansion(mult = c(0.04, 0.12))
      ) +
      
      labs(
        title = taxon_name,
        subtitle = paste0(
          "Group × time interaction: FDR = ", signif(fdr, 3),
          " | log2FC interaction = ", signif(log2fc, 3)
        ),
        x = NULL,
        y = "Relative abundance (%)",
        color = NULL,
        fill = NULL
      ) +
      
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(face = "italic", size = 14),
        plot.subtitle = element_text(size = 9.5),
        legend.position = "top",
        axis.title.y = element_text(margin = margin(r = 8))
      )
    
    print(p)
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", taxon_name)
    safe_name <- substr(safe_name, 1, 80)
    
    ggsave(
      filename = file.path(
        plot_dir,
        paste0("DESeq2_Genus_", safe_name, ".png")
      ),
      plot = p,
      width = 7.2,
      height = 5.2,
      dpi = 320
    )
  }
}


########## LONGITUDINAL MORNING vs CNTRL: ALPHA, BETA, DESeq2 ########## 

library(phyloseq)
library(vegan)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(permute)

set.seed(1994)
outdir <- "Results_Cronotipo/Longitudinal_Morning_CNT"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "Alpha"), showWarnings = FALSE)
dir.create(file.path(outdir, "Beta"), showWarnings = FALSE)
dir.create(file.path(outdir, "DESeq2"), showWarnings = FALSE)
sample_data(data)
COLORS <- c(control = "gray", morning = "#FFB347")
GROUP_LEVELS <- c("control", "morning")
TIME_LEVELS <- c("T0", "T1")


#### 1. Select partecipants ####

meta0 <- data.frame(sample_data(data))
required_cols <- c("Sample_Name", "Time", "Type_diet")
if (!all(required_cols %in% colnames(meta0))) stop("Mancano metadati richiesti: Sample_Name, Time o Type_diet.")

meta0$Sample_Name <- as.character(meta0$Sample_Name)
meta0$Time <- as.character(meta0$Time)
meta0$Type_diet <- tolower(as.character(meta0$Type_diet))

meta0 <- meta0[meta0$Type_diet %in% GROUP_LEVELS & meta0$Time %in% TIME_LEVELS, , drop = FALSE]
if (nrow(meta0) == 0) stop("Problem")

# Assunzione di questo script: ogni partecipante appartiene a un solo gruppo
groups_per_subject <- tapply(meta0$Type_diet, meta0$Sample_Name, function(x) length(unique(x)))
if (any(groups_per_subject > 1)) stop("Problem")

# Mantieni solo partecipanti con esattamente un campione a T0 e uno a T1
tab <- table(meta0$Sample_Name, factor(meta0$Time, levels = TIME_LEVELS))
paired_ids <- rownames(tab)[tab[, "T0"] == 1 & tab[, "T1"] == 1]
meta0 <- meta0[meta0$Sample_Name %in% paired_ids, , drop = FALSE]
if (length(paired_ids) < 2) stop("Few paired partecipants")

ps <- prune_samples(rownames(meta0), data)
sample_data(ps)$Sample_Name <- factor(as.character(sample_data(ps)$Sample_Name))
sample_data(ps)$Time <- factor(as.character(sample_data(ps)$Time), levels = TIME_LEVELS)
sample_data(ps)$Type_diet <- factor(tolower(as.character(sample_data(ps)$Type_diet)), levels = GROUP_LEVELS)

meta <- data.frame(sample_data(ps))
write.csv(as.data.frame(table(meta$Type_diet, meta$Time)), file.path(outdir, "Sample_counts_by_group_time.csv"), row.names = FALSE)

cat("\nSamples for group and time:\n")
print(table(sample_data(ps)$Type_diet, sample_data(ps)$Time))
cat("\nSamples for group and time\n")
print(table(tapply(as.character(sample_data(ps)$Type_diet), sample_data(ps)$Sample_Name, unique)))


#### 2. ALPHA DIVERSITY #### 

ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
alpha <- estimate_richness(ps_genus, measures = c("Observed", "Shannon"))
alpha$SampleID <- rownames(alpha)
md <- data.frame(sample_data(ps_genus))
md$SampleID <- rownames(md)
alpha <- left_join(alpha, md, by = "SampleID")
alpha$Observed <- as.numeric(alpha$Observed)
alpha$Shannon <- as.numeric(alpha$Shannon)
alpha$Evenness <- ifelse(alpha$Observed > 1, alpha$Shannon / log(alpha$Observed), NA_real_)
alpha$Time <- factor(as.character(alpha$Time), levels = TIME_LEVELS)
alpha$Type_diet <- factor(as.character(alpha$Type_diet), levels = GROUP_LEVELS)
alpha$Sample_Name <- factor(alpha$Sample_Name)

alpha_metrics <- c("Observed", "Shannon", "Evenness")
alpha_models <- list()
alpha_interaction_results <- list()
alpha_change_results <- list()

for (metric in alpha_metrics) {
  dat <- alpha[!is.na(alpha[[metric]]) & is.finite(alpha[[metric]]), , drop = FALSE]
  form_full <- as.formula(paste0(metric, " ~ Type_diet * Time + (1 | Sample_Name)"))
  form_reduced <- as.formula(paste0(metric, " ~ Type_diet + Time + (1 | Sample_Name)"))
  model_full <- lmer(form_full, data = dat, REML = FALSE)
  model_reduced <- lmer(form_reduced, data = dat, REML = FALSE)
  lrt <- anova(model_reduced, model_full)
  p_int <- lrt$`Pr(>Chisq)`[2]
  alpha_models[[metric]] <- model_full
  alpha_interaction_results[[metric]] <- data.frame(Metric = metric, Chi_square = lrt$Chisq[2], df = lrt$Df[2], p_value = lrt$`Pr(>Chisq)`[2])
  
  emm_time <- emmeans(model_full, ~ Time | Type_diet)
  changes <- as.data.frame(contrast(emm_time, method = "revpairwise", adjust = "none"))
  changes$Metric <- metric
  alpha_change_results[[metric]] <- changes
  
  p <- ggplot(dat, aes(x = Time, y = .data[[metric]])) +
    geom_line(aes(group = Sample_Name), color = "grey65", alpha = 0.55, linewidth = 0.45) +
    geom_boxplot(aes(group = interaction(Type_diet, Time), fill = Type_diet), width = 0.52, alpha = 0.22, outlier.shape = NA, position = position_dodge(width = 0.16)) +
    geom_point(aes(color = Type_diet), position = position_jitter(width = 0.035, height = 0, seed = 1994), size = 2.1, alpha = 0.85) +
    stat_summary(aes(group = Type_diet, color = Type_diet), fun = mean, geom = "line", linewidth = 1.05, position = position_dodge(width = 0.05)) +
    stat_summary(aes(group = Type_diet, color = Type_diet), fun = mean, geom = "point", size = 3.1, shape = 18) +
    scale_color_manual(values = COLORS, breaks = GROUP_LEVELS) +
    scale_fill_manual(values = COLORS, breaks = GROUP_LEVELS) +
    labs(title = paste0(metric, " diversity at genus level"), subtitle = paste0("Group × time interaction: p = ", signif(p_int, 3)), x = NULL, y = metric, color = NULL, fill = NULL) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "plain", size = 14), plot.subtitle = element_text(size = 10), legend.position = "top", axis.title.y = element_text(margin = margin(r = 8)))
  print(p)
  ggsave(file.path(outdir, "Alpha", paste0("Alpha_", metric, ".png")), p, width = 7.2, height = 5.2, dpi = 320)
}
p

alpha_interaction_df <- bind_rows(alpha_interaction_results)
alpha_interaction_df$p_adj_BH <- p.adjust(alpha_interaction_df$p_value, method = "BH")
write.csv(alpha_interaction_df, file.path(outdir, "Alpha", "Alpha_interaction_tests.csv"), row.names = FALSE)
write.csv(bind_rows(alpha_change_results), file.path(outdir, "Alpha", "Alpha_T1_vs_T0_within_group.csv"), row.names = FALSE)


#### 3. BETA DIVERSITY ####

ps_beta <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
ps_rel <- transform_sample_counts(ps_beta, function(x) if (sum(x) > 0) x / sum(x) else x)
ps_sqrt <- transform_sample_counts(ps_rel, sqrt)
dist_bc <- phyloseq::distance(ps_sqrt, method = "bray")
meta_beta <- data.frame(sample_data(ps_sqrt))
meta_beta <- meta_beta[sample_names(ps_sqrt), , drop = FALSE]
meta_beta$Sample_Name <- factor(meta_beta$Sample_Name)
meta_beta$Time <- factor(as.character(meta_beta$Time), levels = TIME_LEVELS)
meta_beta$Type_diet <- factor(as.character(meta_beta$Type_diet), levels = GROUP_LEVELS)

perm_control <- permute::how(nperm = 9999, blocks = meta_beta$Sample_Name)
perm_beta <- adonis2(dist_bc ~ Type_diet * Time, data = meta_beta, permutations = perm_control, by = "margin")
print(perm_beta)
write.csv(as.data.frame(perm_beta), file.path(outdir, "Beta", "Beta_PERMANOVA.csv"), row.names = TRUE)

interaction_rows <- grep(":", rownames(perm_beta))
if (length(interaction_rows) == 0) interaction_rows <- grep("Time:Type_diet|Type_diet:Time", rownames(perm_beta))
p_beta_int <- if (length(interaction_rows) > 0) perm_beta$`Pr(>F)`[interaction_rows[1]] else NA_real_

meta_beta$Group_Time <- interaction(meta_beta$Type_diet, meta_beta$Time, drop = TRUE)
bd <- betadisper(dist_bc, group = meta_beta$Group_Time)
disp_test <- permutest(bd, permutations = perm_control)
print(disp_test)
write.csv(as.data.frame(disp_test$tab), file.path(outdir, "Beta", "Beta_dispersion_test.csv"), row.names = TRUE)

ord <- ordinate(ps_sqrt, method = "PCoA", distance = dist_bc)
ord_df <- as.data.frame(ord$vectors[, 1:2, drop = FALSE])
ord_df$SampleID <- rownames(ord_df)
ord_df <- left_join(ord_df, tibble::rownames_to_column(meta_beta, "SampleID"), by = "SampleID")
eig <- ord$values$Eigenvalues
if (is.null(eig)) eig <- rep(NA_real_, 2)
positive_eig <- eig[eig > 0]
var1 <- if (length(positive_eig) >= 1 && sum(positive_eig) > 0) round(100 * positive_eig[1] / sum(positive_eig), 1) else NA
var2 <- if (length(positive_eig) >= 2 && sum(positive_eig) > 0) round(100 * positive_eig[2] / sum(positive_eig), 1) else NA

p_beta <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2)) +
  geom_line(aes(group = Sample_Name), color = "grey65", alpha = 0.5, linewidth = 0.45) +
  geom_point(aes(color = Type_diet, shape = Time), size = 3, alpha = 0.9) +
  scale_color_manual(values = COLORS, breaks = GROUP_LEVELS) +
  scale_shape_manual(values = c(T0 = 16, T1 = 17)) +
  labs(title = "Genus-level beta diversity", subtitle = paste0("Bray–Curtis PCoA | Group × time PERMANOVA p = ", signif(p_beta_int, 3)), x = paste0("PCoA1", if (!is.na(var1)) paste0(" (", var1, "%)") else ""), y = paste0("PCoA2", if (!is.na(var2)) paste0(" (", var2, "%)") else ""), color = NULL, shape = "Time") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "plain", size = 14), plot.subtitle = element_text(size = 10), legend.position = "top")
print(p_beta)
ggsave(file.path(outdir, "Beta", "Beta_PCoA_longitudinal.png"), p_beta, width = 7.5, height = 6, dpi = 320)

### 4. DESeq2 #####

ps_dds <- ps
sample_data(ps_dds)$Time <- factor(as.character(sample_data(ps_dds)$Time), levels = TIME_LEVELS)
sample_data(ps_dds)$Type_diet <- factor(as.character(sample_data(ps_dds)$Type_diet), levels = GROUP_LEVELS)
sample_data(ps_dds)$Sample_Name <- factor(as.character(sample_data(ps_dds)$Sample_Name))
sample_data(ps_dds)$Interaction <- factor(ifelse(sample_data(ps_dds)$Time == "T1" & sample_data(ps_dds)$Type_diet == "morning", "morning_T1", "other"), levels = c("other", "morning_T1"))

meta_dds <- data.frame(sample_data(ps_dds))
X <- model.matrix(~ Sample_Name + Time + Interaction, data = meta_dds)
if (qr(X)$rank != ncol(X)) stop("Deseq2 model is not full rank: control")

tax_ranks <- intersect(c("Phylum", "Class", "Order", "Family", "Genus"), colnames(tax_table(ps_dds)))
if (length(tax_ranks) == 0) stop("Problem")

deseq_results <- list()
for (rank in tax_ranks) {
  message("DESeq2 - livello tassonomico: ", rank)
  ps_rank <- tax_glom(ps_dds, taxrank = rank, NArm = FALSE)
  dds <- phyloseq_to_deseq2(ps_rank, ~ Sample_Name + Time + Interaction)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  if (nrow(dds) == 0) next
  dds <- DESeq(dds, test = "Wald", quiet = TRUE)
  coef_interaction <- grep("^Interaction_", resultsNames(dds), value = TRUE)
  if (length(coef_interaction) != 1) {
    warning("No Interaction Coefficient", rank, ". resultsNames: ", paste(resultsNames(dds), collapse = ", "))
    next
  }
  res <- as.data.frame(results(dds, name = coef_interaction, alpha = 0.05))
  res$Taxon_ID <- rownames(res)
  res$Taxonomic_rank <- rank
  tax_df <- as.data.frame(tax_table(ps_rank))
  tax_df$Taxon_ID <- rownames(tax_df)
  res$Taxon_name <- as.character(tax_df[[rank]][match(res$Taxon_ID, tax_df$Taxon_ID)])
  missing_name <- is.na(res$Taxon_name) | res$Taxon_name == ""
  res$Taxon_name[missing_name] <- paste0("Unclassified_", rank, "_", res$Taxon_ID[missing_name])
  deseq_results[[rank]] <- res
}

if (length(deseq_results) > 0) {
  deseq_all <- bind_rows(deseq_results)
  deseq_all$padj_global_BH <- p.adjust(deseq_all$pvalue, method = "BH")
  deseq_all <- deseq_all %>% arrange(padj_global_BH, pvalue)
  write.csv(deseq_all, file.path(outdir, "DESeq2", "DESeq2_all_interaction_results.csv"), row.names = FALSE)
  sig <- deseq_all %>% filter(!is.na(padj_global_BH), padj_global_BH < 0.05)
  write.csv(sig, file.path(outdir, "DESeq2", "DESeq2_significant_interaction_results.csv"), row.names = FALSE)
} 


################### DESEQ2: GRAPHS ##################
COLORS <- c("control" = "gray", "morning" = "#FFB347")
GROUP_LEVELS <- c("control", "morning")
TIME_LEVELS <- c("T0", "T1")

plot_dir <- file.path(outdir, "DESeq2", "Individual_Taxa_Plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

sig_taxa <- sig %>%
  filter(
    !is.na(padj_global_BH),
    padj_global_BH < 0.05,
    Taxonomic_rank %in% c("Phylum", "Class", "Order", "Family", "Genus")
  )

if (nrow(sig_taxa) == 0) {
  message("No significance after BH")
} else {
  
  for (rank in unique(sig_taxa$Taxonomic_rank)) {
    
    message("Graphs for: ", rank)
    ps_rank <- tax_glom(ps_dds, taxrank = rank, NArm = FALSE)
    ps_rank_rel <- transform_sample_counts(
      ps_rank,
      function(x) if (sum(x) > 0) x / sum(x) * 100 else x
    )
    sig_rank <- sig_taxa %>%
      filter(Taxonomic_rank == rank)
    taxa_keep <- intersect(taxa_names(ps_rank_rel), sig_rank$Taxon_ID)
    
    if (length(taxa_keep) == 0) next
    
    ps_rank_rel <- prune_taxa(taxa_keep, ps_rank_rel)
    plot_data <- psmelt(ps_rank_rel) %>%
      mutate(
        Time = factor(as.character(Time), levels = TIME_LEVELS),
        Type_diet = factor(as.character(Type_diet), levels = GROUP_LEVELS),
        Sample_Name = as.character(Sample_Name),
        Time_num = as.numeric(Time),
        x_plot = Time_num + ifelse(Type_diet == "control", -0.09, 0.09)
      ) %>%
      filter(!is.na(Time), !is.na(Type_diet))
    for (taxon_id in taxa_keep) {
      
      result_taxon <- sig_rank %>%
        filter(Taxon_ID == taxon_id)
      
      taxon_name <- result_taxon$Taxon_name[1]
      fdr <- result_taxon$padj_global_BH[1]
      log2fc <- result_taxon$log2FoldChange[1]
      
      dat_taxon <- plot_data %>%
        filter(OTU == taxon_id)
      
      if (nrow(dat_taxon) == 0) next
      
      p <- ggplot(dat_taxon, aes(x = x_plot, y = Abundance)) +
        geom_line(
          aes(group = Sample_Name),
          color = "grey65",
          alpha = 0.55,
          linewidth = 0.45
        ) +
        geom_boxplot(
          aes(group = interaction(Type_diet, Time), fill = Type_diet),
          width = 0.13,
          alpha = 0.22,
          outlier.shape = NA,
          linewidth = 0.45
        ) +
        geom_point(
          aes(color = Type_diet),
          position = position_jitter(width = 0.025, height = 0, seed = 1994),
          size = 2.1,
          alpha = 0.85
        ) +
        stat_summary(
          aes(group = Type_diet, color = Type_diet),
          fun = mean,
          geom = "line",
          linewidth = 1.05
        ) +
        
        stat_summary(
          aes(group = Type_diet, color = Type_diet),
          fun = mean,
          geom = "point",
          size = 3.1,
          shape = 18
        ) +
        
        scale_color_manual(values = COLORS, breaks = GROUP_LEVELS) +
        scale_fill_manual(values = COLORS, breaks = GROUP_LEVELS) +
        
        scale_x_continuous(
          breaks = c(1, 2),
          labels = TIME_LEVELS,
          expand = expansion(mult = c(0.12, 0.12))
        ) +
        
        scale_y_continuous(
          expand = expansion(mult = c(0.04, 0.12))
        ) +
        
        labs(
          title = taxon_name,
          subtitle = paste0(
            rank,
            " | Group × time interaction: FDR = ", signif(fdr, 3),
            " | log2FC interaction = ", signif(log2fc, 3)
          ),
          x = NULL,
          y = "Relative abundance (%)",
          color = NULL,
          fill = NULL
        ) +
        
        theme_classic(base_size = 12) +
        theme(
          plot.title = element_text(face = "italic", size = 14),
          plot.subtitle = element_text(size = 9.5),
          legend.position = "top",
          axis.title.y = element_text(margin = margin(r = 8))
        )
      
      print(p)
      safe_name <- gsub("[^A-Za-z0-9_-]", "_", taxon_name)
      safe_name <- substr(safe_name, 1, 80)
      
      ggsave(
        filename = file.path(
          plot_dir,
          paste0("DESeq2_", rank, "_", safe_name, ".png")
        ),
        plot = p,
        width = 7.2,
        height = 5.2,
        dpi = 320
      )
    }
  }
}


########## LONGITUDINAL EVENING vs CINTRL: ALPHA, BETA, DESeq2 ########## 

library(phyloseq)
library(vegan)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(permute)

set.seed(1994)
outdir <- "Results_Cronotipo/Longitudinal_evening_CNT"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "Alpha"), showWarnings = FALSE)
dir.create(file.path(outdir, "Beta"), showWarnings = FALSE)
dir.create(file.path(outdir, "DESeq2"), showWarnings = FALSE)
sample_data(data)
COLORS <- c(control = "gray", evening = "blue")
GROUP_LEVELS <- c("control", "evening")
TIME_LEVELS <- c("T0", "T1")

#### 1. Select partecipants ####
meta0 <- data.frame(sample_data(data))
required_cols <- c("Sample_Name", "Time", "Type_diet")
if (!all(required_cols %in% colnames(meta0))) stop("No metadata for: Sample_Name, Time o Type_diet.")

meta0$Sample_Name <- as.character(meta0$Sample_Name)
meta0$Time <- as.character(meta0$Time)
meta0$Type_diet <- tolower(as.character(meta0$Type_diet))

meta0 <- meta0[meta0$Type_diet %in% GROUP_LEVELS & meta0$Time %in% TIME_LEVELS, , drop = FALSE]
if (nrow(meta0) == 0) stop("Problem")

groups_per_subject <- tapply(meta0$Type_diet, meta0$Sample_Name, function(x) length(unique(x)))
if (any(groups_per_subject > 1)) stop("Problem")
tab <- table(meta0$Sample_Name, factor(meta0$Time, levels = TIME_LEVELS))
paired_ids <- rownames(tab)[tab[, "T0"] == 1 & tab[, "T1"] == 1]
meta0 <- meta0[meta0$Sample_Name %in% paired_ids, , drop = FALSE]
if (length(paired_ids) < 2) stop("Few paired partecipants")

ps <- prune_samples(rownames(meta0), data)
sample_data(ps)$Sample_Name <- factor(as.character(sample_data(ps)$Sample_Name))
sample_data(ps)$Time <- factor(as.character(sample_data(ps)$Time), levels = TIME_LEVELS)
sample_data(ps)$Type_diet <- factor(tolower(as.character(sample_data(ps)$Type_diet)), levels = GROUP_LEVELS)

meta <- data.frame(sample_data(ps))
write.csv(as.data.frame(table(meta$Type_diet, meta$Time)), file.path(outdir, "Sample_counts_by_group_time.csv"), row.names = FALSE)

cat("\nSamples for group and time:\n")
print(table(sample_data(ps)$Type_diet, sample_data(ps)$Time))
cat("\nNumber of partecipants:\n")
print(table(tapply(as.character(sample_data(ps)$Type_diet), sample_data(ps)$Sample_Name, unique)))

##### 2. ALPHA DIVERSITY ######
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
alpha <- estimate_richness(ps_genus, measures = c("Observed", "Shannon"))
alpha$SampleID <- rownames(alpha)
md <- data.frame(sample_data(ps_genus))
md$SampleID <- rownames(md)
alpha <- left_join(alpha, md, by = "SampleID")
alpha$Observed <- as.numeric(alpha$Observed)
alpha$Shannon <- as.numeric(alpha$Shannon)
alpha$Evenness <- ifelse(alpha$Observed > 1, alpha$Shannon / log(alpha$Observed), NA_real_)
alpha$Time <- factor(as.character(alpha$Time), levels = TIME_LEVELS)
alpha$Type_diet <- factor(as.character(alpha$Type_diet), levels = GROUP_LEVELS)
alpha$Sample_Name <- factor(alpha$Sample_Name)

alpha_metrics <- c("Observed", "Shannon", "Evenness")
alpha_models <- list()
alpha_interaction_results <- list()
alpha_change_results <- list()

for (metric in alpha_metrics) {
  dat <- alpha[!is.na(alpha[[metric]]) & is.finite(alpha[[metric]]), , drop = FALSE]
  form_full <- as.formula(paste0(metric, " ~ Type_diet * Time + (1 | Sample_Name)"))
  form_reduced <- as.formula(paste0(metric, " ~ Type_diet + Time + (1 | Sample_Name)"))
  model_full <- lmer(form_full, data = dat, REML = FALSE)
  model_reduced <- lmer(form_reduced, data = dat, REML = FALSE)
  lrt <- anova(model_reduced, model_full)
  p_int <- lrt$`Pr(>Chisq)`[2]
  alpha_models[[metric]] <- model_full
  alpha_interaction_results[[metric]] <- data.frame(Metric = metric, Chi_square = lrt$Chisq[2], df = lrt$Df[2], p_value = lrt$`Pr(>Chisq)`[2])
  
  emm_time <- emmeans(model_full, ~ Time | Type_diet)
  changes <- as.data.frame(contrast(emm_time, method = "revpairwise", adjust = "none"))
  changes$Metric <- metric
  alpha_change_results[[metric]] <- changes
  
  p <- ggplot(dat, aes(x = Time, y = .data[[metric]])) +
    geom_line(aes(group = Sample_Name), color = "grey65", alpha = 0.55, linewidth = 0.45) +
    geom_boxplot(aes(group = interaction(Type_diet, Time), fill = Type_diet), width = 0.52, alpha = 0.22, outlier.shape = NA, position = position_dodge(width = 0.16)) +
    geom_point(aes(color = Type_diet), position = position_jitter(width = 0.035, height = 0, seed = 1994), size = 2.1, alpha = 0.85) +
    stat_summary(aes(group = Type_diet, color = Type_diet), fun = mean, geom = "line", linewidth = 1.05, position = position_dodge(width = 0.05)) +
    stat_summary(aes(group = Type_diet, color = Type_diet), fun = mean, geom = "point", size = 3.1, shape = 18) +
    scale_color_manual(values = COLORS, breaks = GROUP_LEVELS) +
    scale_fill_manual(values = COLORS, breaks = GROUP_LEVELS) +
    labs(title = paste0(metric, " diversity at genus level"), subtitle = paste0("Group × time interaction: p = ", signif(p_int, 3)), x = NULL, y = metric, color = NULL, fill = NULL) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "plain", size = 14), plot.subtitle = element_text(size = 10), legend.position = "top", axis.title.y = element_text(margin = margin(r = 8)))
  print(p)
  ggsave(file.path(outdir, "Alpha", paste0("Alpha_", metric, ".png")), p, width = 7.2, height = 5.2, dpi = 320)
}
p

alpha_interaction_df <- bind_rows(alpha_interaction_results)
alpha_interaction_df$p_adj_BH <- p.adjust(alpha_interaction_df$p_value, method = "BH")
write.csv(alpha_interaction_df, file.path(outdir, "Alpha", "Alpha_interaction_tests.csv"), row.names = FALSE)
write.csv(bind_rows(alpha_change_results), file.path(outdir, "Alpha", "Alpha_T1_vs_T0_within_group.csv"), row.names = FALSE)


##### 3. BETA DIVERSITY ##### 
ps_beta <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
ps_rel <- transform_sample_counts(ps_beta, function(x) if (sum(x) > 0) x / sum(x) else x)
ps_sqrt <- transform_sample_counts(ps_rel, sqrt)
dist_bc <- phyloseq::distance(ps_sqrt, method = "bray")
meta_beta <- data.frame(sample_data(ps_sqrt))
meta_beta <- meta_beta[sample_names(ps_sqrt), , drop = FALSE]
meta_beta$Sample_Name <- factor(meta_beta$Sample_Name)
meta_beta$Time <- factor(as.character(meta_beta$Time), levels = TIME_LEVELS)
meta_beta$Type_diet <- factor(as.character(meta_beta$Type_diet), levels = GROUP_LEVELS)

perm_control <- permute::how(nperm = 9999, blocks = meta_beta$Sample_Name)
perm_beta <- adonis2(dist_bc ~ Type_diet * Time, data = meta_beta, permutations = perm_control, by = "margin")
print(perm_beta)
write.csv(as.data.frame(perm_beta), file.path(outdir, "Beta", "Beta_PERMANOVA.csv"), row.names = TRUE)

interaction_rows <- grep(":", rownames(perm_beta))
if (length(interaction_rows) == 0) interaction_rows <- grep("Time:Type_diet|Type_diet:Time", rownames(perm_beta))
p_beta_int <- if (length(interaction_rows) > 0) perm_beta$`Pr(>F)`[interaction_rows[1]] else NA_real_

meta_beta$Group_Time <- interaction(meta_beta$Type_diet, meta_beta$Time, drop = TRUE)
bd <- betadisper(dist_bc, group = meta_beta$Group_Time)
disp_test <- permutest(bd, permutations = perm_control)
print(disp_test)
write.csv(as.data.frame(disp_test$tab), file.path(outdir, "Beta", "Beta_dispersion_test.csv"), row.names = TRUE)

ord <- ordinate(ps_sqrt, method = "PCoA", distance = dist_bc)
ord_df <- as.data.frame(ord$vectors[, 1:2, drop = FALSE])
ord_df$SampleID <- rownames(ord_df)
ord_df <- left_join(ord_df, tibble::rownames_to_column(meta_beta, "SampleID"), by = "SampleID")
eig <- ord$values$Eigenvalues
if (is.null(eig)) eig <- rep(NA_real_, 2)
positive_eig <- eig[eig > 0]
var1 <- if (length(positive_eig) >= 1 && sum(positive_eig) > 0) round(100 * positive_eig[1] / sum(positive_eig), 1) else NA
var2 <- if (length(positive_eig) >= 2 && sum(positive_eig) > 0) round(100 * positive_eig[2] / sum(positive_eig), 1) else NA

p_beta <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2)) +
  geom_line(aes(group = Sample_Name), color = "grey65", alpha = 0.5, linewidth = 0.45) +
  geom_point(aes(color = Type_diet, shape = Time), size = 3, alpha = 0.9) +
  scale_color_manual(values = COLORS, breaks = GROUP_LEVELS) +
  scale_shape_manual(values = c(T0 = 16, T1 = 17)) +
  labs(title = "Genus-level beta diversity", subtitle = paste0("Bray–Curtis PCoA | Group × time PERMANOVA p = ", signif(p_beta_int, 3)), x = paste0("PCoA1", if (!is.na(var1)) paste0(" (", var1, "%)") else ""), y = paste0("PCoA2", if (!is.na(var2)) paste0(" (", var2, "%)") else ""), color = NULL, shape = "Time") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "plain", size = 14), plot.subtitle = element_text(size = 10), legend.position = "top")
print(p_beta)
ggsave(file.path(outdir, "Beta", "Beta_PCoA_longitudinal.png"), p_beta, width = 7.5, height = 6, dpi = 320)


######## 4. DESeq2 #########
ps_dds <- ps
sample_data(ps_dds)$Time <- factor(as.character(sample_data(ps_dds)$Time), levels = TIME_LEVELS)
sample_data(ps_dds)$Type_diet <- factor(as.character(sample_data(ps_dds)$Type_diet), levels = GROUP_LEVELS)
sample_data(ps_dds)$Sample_Name <- factor(as.character(sample_data(ps_dds)$Sample_Name))
sample_data(ps_dds)$Interaction <- factor(ifelse(sample_data(ps_dds)$Time == "T1" & sample_data(ps_dds)$Type_diet == "evening", "evening_T1", "other"), levels = c("other", "evening_T1"))

meta_dds <- data.frame(sample_data(ps_dds))
X <- model.matrix(~ Sample_Name + Time + Interaction, data = meta_dds)
if (qr(X)$rank != ncol(X)) stop("DESeq2 model is not full rank")

tax_ranks <- intersect(c("Phylum", "Class", "Order", "Family", "Genus"), colnames(tax_table(ps_dds)))
if (length(tax_ranks) == 0) stop("Problem")

deseq_results <- list()
for (rank in tax_ranks) {
  message("DESeq2 - Taxa level: ", rank)
  ps_rank <- tax_glom(ps_dds, taxrank = rank, NArm = FALSE)
  dds <- phyloseq_to_deseq2(ps_rank, ~ Sample_Name + Time + Interaction)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  if (nrow(dds) == 0) next
  dds <- DESeq(dds, test = "Wald", quiet = TRUE)
  coef_interaction <- grep("^Interaction_", resultsNames(dds), value = TRUE)
  if (length(coef_interaction) != 1) {
    warning("No Interaction coefficent for", rank, ". resultsNames: ", paste(resultsNames(dds), collapse = ", "))
    next
  }
  res <- as.data.frame(results(dds, name = coef_interaction, alpha = 0.05))
  res$Taxon_ID <- rownames(res)
  res$Taxonomic_rank <- rank
  tax_df <- as.data.frame(tax_table(ps_rank))
  tax_df$Taxon_ID <- rownames(tax_df)
  res$Taxon_name <- as.character(tax_df[[rank]][match(res$Taxon_ID, tax_df$Taxon_ID)])
  missing_name <- is.na(res$Taxon_name) | res$Taxon_name == ""
  res$Taxon_name[missing_name] <- paste0("Unclassified_", rank, "_", res$Taxon_ID[missing_name])
  deseq_results[[rank]] <- res
}

if (length(deseq_results) > 0) {
  deseq_all <- bind_rows(deseq_results)
  deseq_all$padj_global_BH <- p.adjust(deseq_all$pvalue, method = "BH")
  deseq_all <- deseq_all %>% arrange(padj_global_BH, pvalue)
  write.csv(deseq_all, file.path(outdir, "DESeq2", "DESeq2_all_interaction_results.csv"), row.names = FALSE)
  sig <- deseq_all %>% filter(!is.na(padj_global_BH), padj_global_BH < 0.05)
  write.csv(sig, file.path(outdir, "DESeq2", "DESeq2_significant_interaction_results.csv"), row.names = FALSE)
} 


########## LONGITUDINAL intervention vs CONTRL: ALPHA, BETA, DESeq2 ########## 

library(phyloseq)
library(vegan)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(permute)

set.seed(1994)
outdir <- "Results_Cronotipo/Longitudinal_intervention_CNT"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "Alpha"), showWarnings = FALSE)
dir.create(file.path(outdir, "Beta"), showWarnings = FALSE)
dir.create(file.path(outdir, "DESeq2"), showWarnings = FALSE)
sample_data(data)
COLORS <- c(control = "gray", intervention = "red")
GROUP_LEVELS <- c("control", "intervention")
TIME_LEVELS <- c("T0", "T1")

##### 1. Select partecipants ########

meta0 <- data.frame(sample_data(data))
required_cols <- c("Sample_Name", "Time", "Group")
if (!all(required_cols %in% colnames(meta0))) stop("No metadata: Sample_Name, Time o Group.")

meta0$Sample_Name <- as.character(meta0$Sample_Name)
meta0$Time <- as.character(meta0$Time)
meta0$Group <- tolower(as.character(meta0$Group))

meta0 <- meta0[meta0$Group %in% GROUP_LEVELS & meta0$Time %in% TIME_LEVELS, , drop = FALSE]
if (nrow(meta0) == 0) stop("Problem")

groups_per_subject <- tapply(meta0$Group, meta0$Sample_Name, function(x) length(unique(x)))
if (any(groups_per_subject > 1)) stop("Problem")
tab <- table(meta0$Sample_Name, factor(meta0$Time, levels = TIME_LEVELS))
paired_ids <- rownames(tab)[tab[, "T0"] == 1 & tab[, "T1"] == 1]
meta0 <- meta0[meta0$Sample_Name %in% paired_ids, , drop = FALSE]
if (length(paired_ids) < 2) stop("Few paired partecipants")

ps <- prune_samples(rownames(meta0), data)
sample_data(ps)$Sample_Name <- factor(as.character(sample_data(ps)$Sample_Name))
sample_data(ps)$Time <- factor(as.character(sample_data(ps)$Time), levels = TIME_LEVELS)
sample_data(ps)$Group <- factor(tolower(as.character(sample_data(ps)$Group)), levels = GROUP_LEVELS)

meta <- data.frame(sample_data(ps))
write.csv(as.data.frame(table(meta$Group, meta$Time)), file.path(outdir, "Sample_counts_by_group_time.csv"), row.names = FALSE)

cat("\nSamples for group and time:\n")
print(table(sample_data(ps)$Group, sample_data(ps)$Time))
cat("\nNumber of partecipants:\n")
print(table(tapply(as.character(sample_data(ps)$Group), sample_data(ps)$Sample_Name, unique)))

######## 2. ALPHA DIVERSITY #######

ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
alpha <- estimate_richness(ps_genus, measures = c("Observed", "Shannon"))
alpha$SampleID <- rownames(alpha)
md <- data.frame(sample_data(ps_genus))
md$SampleID <- rownames(md)
alpha <- left_join(alpha, md, by = "SampleID")
alpha$Observed <- as.numeric(alpha$Observed)
alpha$Shannon <- as.numeric(alpha$Shannon)
alpha$Evenness <- ifelse(alpha$Observed > 1, alpha$Shannon / log(alpha$Observed), NA_real_)
alpha$Time <- factor(as.character(alpha$Time), levels = TIME_LEVELS)
alpha$Group <- factor(as.character(alpha$Group), levels = GROUP_LEVELS)
alpha$Sample_Name <- factor(alpha$Sample_Name)

alpha_metrics <- c("Observed", "Shannon", "Evenness")
alpha_models <- list()
alpha_interaction_results <- list()
alpha_change_results <- list()

for (metric in alpha_metrics) {
  dat <- alpha[!is.na(alpha[[metric]]) & is.finite(alpha[[metric]]), , drop = FALSE]
  form_full <- as.formula(paste0(metric, " ~ Group * Time + (1 | Sample_Name)"))
  form_reduced <- as.formula(paste0(metric, " ~ Group + Time + (1 | Sample_Name)"))
  model_full <- lmer(form_full, data = dat, REML = FALSE)
  model_reduced <- lmer(form_reduced, data = dat, REML = FALSE)
  lrt <- anova(model_reduced, model_full)
  p_int <- lrt$`Pr(>Chisq)`[2]
  alpha_models[[metric]] <- model_full
  alpha_interaction_results[[metric]] <- data.frame(Metric = metric, Chi_square = lrt$Chisq[2], df = lrt$Df[2], p_value = lrt$`Pr(>Chisq)`[2])
  
  emm_time <- emmeans(model_full, ~ Time | Group)
  changes <- as.data.frame(contrast(emm_time, method = "revpairwise", adjust = "none"))
  changes$Metric <- metric
  alpha_change_results[[metric]] <- changes
  
  p <- ggplot(dat, aes(x = Time, y = .data[[metric]])) +
    geom_line(aes(group = Sample_Name), color = "grey65", alpha = 0.55, linewidth = 0.45) +
    geom_boxplot(aes(group = interaction(Group, Time), fill = Group), width = 0.52, alpha = 0.22, outlier.shape = NA, position = position_dodge(width = 0.16)) +
    geom_point(aes(color = Group), position = position_jitter(width = 0.035, height = 0, seed = 1994), size = 2.1, alpha = 0.85) +
    stat_summary(aes(group = Group, color = Group), fun = mean, geom = "line", linewidth = 1.05, position = position_dodge(width = 0.05)) +
    stat_summary(aes(group = Group, color = Group), fun = mean, geom = "point", size = 3.1, shape = 18) +
    scale_color_manual(values = COLORS, breaks = GROUP_LEVELS) +
    scale_fill_manual(values = COLORS, breaks = GROUP_LEVELS) +
    labs(title = paste0(metric, " diversity at genus level"), subtitle = paste0("Group × time interaction: p = ", signif(p_int, 3)), x = NULL, y = metric, color = NULL, fill = NULL) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "plain", size = 14), plot.subtitle = element_text(size = 10), legend.position = "top", axis.title.y = element_text(margin = margin(r = 8)))
  print(p)
  ggsave(file.path(outdir, "Alpha", paste0("Alpha_", metric, ".png")), p, width = 7.2, height = 5.2, dpi = 320)
}
p

alpha_interaction_df <- bind_rows(alpha_interaction_results)
alpha_interaction_df$p_adj_BH <- p.adjust(alpha_interaction_df$p_value, method = "BH")
write.csv(alpha_interaction_df, file.path(outdir, "Alpha", "Alpha_interaction_tests.csv"), row.names = FALSE)
write.csv(bind_rows(alpha_change_results), file.path(outdir, "Alpha", "Alpha_T1_vs_T0_within_group.csv"), row.names = FALSE)

###### 3. BETA DIVERSITY ########

ps_beta <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
ps_rel <- transform_sample_counts(ps_beta, function(x) if (sum(x) > 0) x / sum(x) else x)
ps_sqrt <- transform_sample_counts(ps_rel, sqrt)
dist_bc <- phyloseq::distance(ps_sqrt, method = "bray")
meta_beta <- data.frame(sample_data(ps_sqrt))
meta_beta <- meta_beta[sample_names(ps_sqrt), , drop = FALSE]
meta_beta$Sample_Name <- factor(meta_beta$Sample_Name)
meta_beta$Time <- factor(as.character(meta_beta$Time), levels = TIME_LEVELS)
meta_beta$Group <- factor(as.character(meta_beta$Group), levels = GROUP_LEVELS)

perm_control <- permute::how(nperm = 9999, blocks = meta_beta$Sample_Name)
perm_beta <- adonis2(dist_bc ~ Group * Time, data = meta_beta, permutations = perm_control, by = "margin")
print(perm_beta)
write.csv(as.data.frame(perm_beta), file.path(outdir, "Beta", "Beta_PERMANOVA.csv"), row.names = TRUE)

interaction_rows <- grep(":", rownames(perm_beta))
if (length(interaction_rows) == 0) interaction_rows <- grep("Time:Group|Group:Time", rownames(perm_beta))
p_beta_int <- if (length(interaction_rows) > 0) perm_beta$`Pr(>F)`[interaction_rows[1]] else NA_real_

meta_beta$Group_Time <- interaction(meta_beta$Group, meta_beta$Time, drop = TRUE)
bd <- betadisper(dist_bc, group = meta_beta$Group_Time)
disp_test <- permutest(bd, permutations = perm_control)
print(disp_test)
write.csv(as.data.frame(disp_test$tab), file.path(outdir, "Beta", "Beta_dispersion_test.csv"), row.names = TRUE)

ord <- ordinate(ps_sqrt, method = "PCoA", distance = dist_bc)
ord_df <- as.data.frame(ord$vectors[, 1:2, drop = FALSE])
ord_df$SampleID <- rownames(ord_df)
ord_df <- left_join(ord_df, tibble::rownames_to_column(meta_beta, "SampleID"), by = "SampleID")
eig <- ord$values$Eigenvalues
if (is.null(eig)) eig <- rep(NA_real_, 2)
positive_eig <- eig[eig > 0]
var1 <- if (length(positive_eig) >= 1 && sum(positive_eig) > 0) round(100 * positive_eig[1] / sum(positive_eig), 1) else NA
var2 <- if (length(positive_eig) >= 2 && sum(positive_eig) > 0) round(100 * positive_eig[2] / sum(positive_eig), 1) else NA

p_beta <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2)) +
  geom_line(aes(group = Sample_Name), color = "grey65", alpha = 0.5, linewidth = 0.45) +
  geom_point(aes(color = Group, shape = Time), size = 3, alpha = 0.9) +
  scale_color_manual(values = COLORS, breaks = GROUP_LEVELS) +
  scale_shape_manual(values = c(T0 = 16, T1 = 17)) +
  labs(title = "Genus-level beta diversity", subtitle = paste0("Bray–Curtis PCoA | Group × time PERMANOVA p = ", signif(p_beta_int, 3)), x = paste0("PCoA1", if (!is.na(var1)) paste0(" (", var1, "%)") else ""), y = paste0("PCoA2", if (!is.na(var2)) paste0(" (", var2, "%)") else ""), color = NULL, shape = "Time") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "plain", size = 14), plot.subtitle = element_text(size = 10), legend.position = "top")
print(p_beta)
ggsave(file.path(outdir, "Beta", "Beta_PCoA_longitudinal.png"), p_beta, width = 7.5, height = 6, dpi = 320)

######## 4. DESeq2 #########
ps_dds <- ps
sample_data(ps_dds)$Time <- factor(as.character(sample_data(ps_dds)$Time), levels = TIME_LEVELS)
sample_data(ps_dds)$Group <- factor(as.character(sample_data(ps_dds)$Group), levels = GROUP_LEVELS)
sample_data(ps_dds)$Sample_Name <- factor(as.character(sample_data(ps_dds)$Sample_Name))
sample_data(ps_dds)$Interaction <- factor(ifelse(sample_data(ps_dds)$Time == "T1" & sample_data(ps_dds)$Group == "intervention", "intervention_T1", "other"), levels = c("other", "intervention_T1"))

meta_dds <- data.frame(sample_data(ps_dds))
X <- model.matrix(~ Sample_Name + Time + Interaction, data = meta_dds)
if (qr(X)$rank != ncol(X)) stop("DESeq2 model is not full rank")

tax_ranks <- intersect(c("Phylum", "Class", "Order", "Family", "Genus"), colnames(tax_table(ps_dds)))
if (length(tax_ranks) == 0) stop("Problem")

deseq_results <- list()
for (rank in tax_ranks) {
  message("DESeq2 - taxonomic level: ", rank)
  ps_rank <- tax_glom(ps_dds, taxrank = rank, NArm = FALSE)
  dds <- phyloseq_to_deseq2(ps_rank, ~ Sample_Name + Time + Interaction)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  if (nrow(dds) == 0) next
  dds <- DESeq(dds, test = "Wald", quiet = TRUE)
  coef_interaction <- grep("^Interaction_", resultsNames(dds), value = TRUE)
  if (length(coef_interaction) != 1) {
    warning("No Interaction Coefficent for", rank, ". resultsNames: ", paste(resultsNames(dds), collapse = ", "))
    next
  }
  res <- as.data.frame(results(dds, name = coef_interaction, alpha = 0.05))
  res$Taxon_ID <- rownames(res)
  res$Taxonomic_rank <- rank
  tax_df <- as.data.frame(tax_table(ps_rank))
  tax_df$Taxon_ID <- rownames(tax_df)
  res$Taxon_name <- as.character(tax_df[[rank]][match(res$Taxon_ID, tax_df$Taxon_ID)])
  missing_name <- is.na(res$Taxon_name) | res$Taxon_name == ""
  res$Taxon_name[missing_name] <- paste0("Unclassified_", rank, "_", res$Taxon_ID[missing_name])
  deseq_results[[rank]] <- res
}

if (length(deseq_results) > 0) {
  deseq_all <- bind_rows(deseq_results)
  deseq_all$padj_global_BH <- p.adjust(deseq_all$pvalue, method = "BH")
  deseq_all <- deseq_all %>% arrange(padj_global_BH, pvalue)
  write.csv(deseq_all, file.path(outdir, "DESeq2", "DESeq2_all_interaction_results.csv"), row.names = FALSE)
  sig <- deseq_all %>% filter(!is.na(padj_global_BH), padj_global_BH < 0.05)
  write.csv(sig, file.path(outdir, "DESeq2", "DESeq2_significant_interaction_results.csv"), row.names = FALSE)
} 

################### DESEQ2: GRAPHS ###################
COLORS <- c("control" = "gray", "intervention" = "red")
GROUP_LEVELS <- c("control", "intervention")
TIME_LEVELS <- c("T0", "T1")
plot_dir <- file.path(outdir, "DESeq2", "Individual_Taxa_Plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

sig_taxa <- sig %>%
  filter(
    !is.na(padj_global_BH),
    padj_global_BH < 0.05,
    Taxonomic_rank %in% c("Phylum", "Class", "Order", "Family", "Genus")
  )

if (nrow(sig_taxa) == 0) {
  message("No significance after BH")
} else {
  
  for (rank in unique(sig_taxa$Taxonomic_rank)) {
    
    message("Graphs creation for: ", rank)
    
    # Aggregazione al livello tassonomico corrente
    ps_rank <- tax_glom(ps_dds, taxrank = rank, NArm = FALSE)
    
    # Conversione in abbondanze relative percentuali
    ps_rank_rel <- transform_sample_counts(
      ps_rank,
      function(x) if (sum(x) > 0) x / sum(x) * 100 else x
    )
    
    # Risultati significativi del rango corrente
    sig_rank <- sig_taxa %>%
      filter(Taxonomic_rank == rank)
    
    # Mantiene solo i taxa significativi presenti nella phyloseq aggregata
    taxa_keep <- intersect(taxa_names(ps_rank_rel), sig_rank$Taxon_ID)
    
    if (length(taxa_keep) == 0) next
    
    ps_rank_rel <- prune_taxa(taxa_keep, ps_rank_rel)
    
    # Tabella in formato long
    plot_data <- psmelt(ps_rank_rel) %>%
      mutate(
        Time = factor(as.character(Time), levels = TIME_LEVELS),
        Group = factor(as.character(Group), levels = GROUP_LEVELS),
        Sample_Name = as.character(Sample_Name),
        Time_num = as.numeric(Time),
        x_plot = Time_num + ifelse(Group == "control", -0.09, 0.09)
      ) %>%
      filter(!is.na(Time), !is.na(Group))
    
    # Un grafico per ciascun taxon significativo
    for (taxon_id in taxa_keep) {
      
      result_taxon <- sig_rank %>%
        filter(Taxon_ID == taxon_id)
      
      taxon_name <- result_taxon$Taxon_name[1]
      fdr <- result_taxon$padj_global_BH[1]
      log2fc <- result_taxon$log2FoldChange[1]
      
      dat_taxon <- plot_data %>%
        filter(OTU == taxon_id)
      
      if (nrow(dat_taxon) == 0) next
      
      p <- ggplot(dat_taxon, aes(x = x_plot, y = Abundance)) +
        geom_line(
          aes(group = Sample_Name),
          color = "grey65",
          alpha = 0.55,
          linewidth = 0.45
        ) +
        geom_boxplot(
          aes(group = interaction(Group, Time), fill = Group),
          width = 0.13,
          alpha = 0.22,
          outlier.shape = NA,
          linewidth = 0.45
        ) +
        geom_point(
          aes(color = Group),
          position = position_jitter(width = 0.025, height = 0, seed = 1994),
          size = 2.1,
          alpha = 0.85
        ) +

        stat_summary(
          aes(group = Group, color = Group),
          fun = mean,
          geom = "line",
          linewidth = 1.05
        ) +
        
        stat_summary(
          aes(group = Group, color = Group),
          fun = mean,
          geom = "point",
          size = 3.1,
          shape = 18
        ) +
        
        scale_color_manual(values = COLORS, breaks = GROUP_LEVELS) +
        scale_fill_manual(values = COLORS, breaks = GROUP_LEVELS) +
        
        scale_x_continuous(
          breaks = c(1, 2),
          labels = TIME_LEVELS,
          expand = expansion(mult = c(0.12, 0.12))
        ) +
        
        scale_y_continuous(
          expand = expansion(mult = c(0.04, 0.12))
        ) +
        
        labs(
          title = taxon_name,
          subtitle = paste0(
            rank,
            " | Group × time interaction: FDR = ", signif(fdr, 3),
            " | log2FC interaction = ", signif(log2fc, 3)
          ),
          x = NULL,
          y = "Relative abundance (%)",
          color = NULL,
          fill = NULL
        ) +
        
        theme_classic(base_size = 12) +
        theme(
          plot.title = element_text(face = "italic", size = 14),
          plot.subtitle = element_text(size = 9.5),
          legend.position = "top",
          axis.title.y = element_text(margin = margin(r = 8))
        )
      
      print(p)
      safe_name <- gsub("[^A-Za-z0-9_-]", "_", taxon_name)
      safe_name <- substr(safe_name, 1, 80)
      
      ggsave(
        filename = file.path(
          plot_dir,
          paste0("DESeq2_", rank, "_", safe_name, ".png")
        ),
        plot = p,
        width = 7.2,
        height = 5.2,
        dpi = 320
      )
    }
  }
}




