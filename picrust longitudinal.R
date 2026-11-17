############################
### PICRUSt2 LONGITUDINALE
############################
library(lme4)
packageVersion("lme4")
# Pacchetti
packages <- c(
  "Maaslin2",
  "dplyr",
  "stringr",
  "ggplot2",
  "readr"
)


library(Maaslin2)
library(dplyr)
library(stringr)
library(ggplot2)
library(readr)

############################
### IMPORT PICRUST2
############################

a <- read.delim(
  "QIIME_output/exported-feature-table/picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz",
  check.names = TRUE
)

# Ispeziona i primi nomi delle colonne
head(colnames(a))

# Mantieni la correzione dei nomi già usata nel tuo script
idx_X <- stringr::str_starts(colnames(a), pattern = "X")

colnames(a)[idx_X] <- substring(
  colnames(a)[idx_X],
  first = 10
)

# Controlla i nomi risultanti
head(colnames(a))

# Salva gli ID e le descrizioni dei pathway
Descriptions <- a[, c("pathway", "description")]

# Controlla che gli ID siano univoci
stopifnot(!anyDuplicated(Descriptions$pathway))

# Rimuovi le colonne descrittive dalla matrice delle abbondanze
a <- a[, !colnames(a) %in% c("pathway", "description"),
       drop = FALSE]

# Assicura che i valori siano numerici
a[] <- lapply(a, function(x) as.numeric(as.character(x)))

# Imposta gli ID dei pathway come nomi delle righe
rownames(a) <- Descriptions$pathway

# Controlli preliminari
stopifnot(!anyDuplicated(colnames(a)))
stopifnot(all(is.finite(as.matrix(a))))
stopifnot(all(as.matrix(a) >= 0))

dim(a)
head(a[, 1:min(5, ncol(a)), drop = FALSE])


############################
### PREPARAZIONE METADATA
############################
metadata<-read.xlsx('Metadata.xlsx', sheetIndex = 2)
metadata<-metadata[,1:5]
metadata<-na.omit(metadata)
metadata$ID <- gsub(" ", "", metadata$ID)
metadata<-metadata [1:217,1:5]
rows<-metadata$ID
row.names(metadata)<-rows
metadata$Sample_Name <- as.character(metadata$Sample_Name)

# Gli Sample_Name dei campioni devono essere univoci
# L'ID del campione deve essere univoco
stopifnot(!anyDuplicated(metadata$ID))

# Sample_Name identifica il partecipante e può ripetersi
# tra T0 e T1
subject_col <- "Sample_Name"
library(stringr)

colnames(a) <- str_extract(colnames(a), "SL[0-9]+")

# Controlla che gli ID dei campioni corrispondano
if (!setequal(colnames(a), metadata$ID)) {
  stop("Gli ID dei campioni in a non corrispondono a metadata$ID.")
}

# Riordina i metadati secondo l'ordine delle colonne di a
Metadata <- metadata[match(colnames(a), metadata$ID), , drop = FALSE]

# Verifica che l'ordine sia corretto
stopifnot(identical(colnames(a), Metadata$ID))

# Identificativo del partecipante per il modello longitudinale
subject_col <- "Sample_Name"

# Controlla l'allineamento
stopifnot(identical(colnames(a), Metadata$ID))

# Controlla le variabili disponibili
colnames(Metadata)
table(Metadata$Time, useNA = "ifany")
table(Metadata$Group, useNA = "ifany")
table(Metadata$Type_diet, useNA = "ifany")

############################
### CONFIGURAZIONE MODELLO
############################

# Cambia questo nome se il tuo Sample_Name del partecipante
# è registrato in una colonna diversa.
subject_col <- "Sample_Name"

# Verifica che tutte le variabili necessarie esistano
required_vars <- c(
  "Sample_Name",
  "Time",
  "Group",
  "Type_diet",
  subject_col
)

missing_vars <- setdiff(required_vars, colnames(Metadata))

if (length(missing_vars) > 0) {
  stop(
    "Mancano nel metadata queste colonne: ",
    paste(missing_vars, collapse = ", ")
  )
}

# Mantieni T0 e T1 e definisci T0 come riferimento
Metadata$Time <- factor(
  Metadata$Time,
  levels = c("T0", "T1")
)

if (anyNA(Metadata$Time)) {
  stop("Time deve contenere i valori T0 e T1 esattamente come specificati.")
}

# Identificativo del partecipante
Metadata[[subject_col]] <- factor(Metadata[[subject_col]])

# Verifica che i partecipanti abbiano misure ripetute
n_timepoints <- Metadata %>%
  group_by(.data[[subject_col]]) %>%
  summarise(
    n_time = n_distinct(Time),
    .groups = "drop"
  )

table(n_timepoints$n_time)

if (!any(n_timepoints$n_time == 2)) {
  stop(
    "Non risultano partecipanti con misure sia a T0 sia a T1. ",
    "Verifica subject_col e il metadata."
  )
}
############################
### DEFINIZIONE CONFRONTI
############################

comparisons <- list(
  
  intervention_vs_control = list(
    variable = "Group",
    reference = "control",
    comparison = "intervention"
  ),
  
  morning_vs_evening = list(
    variable = "Type_diet",
    reference = "evening",
    comparison = "morning"
  ),
  
  morning_vs_control = list(
    variable = "Type_diet",
    reference = "control",
    comparison = "morning"
  ),
  
  evening_vs_control = list(
    variable = "Type_diet",
    reference = "control",
    comparison = "evening"
  )
  
)

############################
### FUNZIONE DI ANALISI
############################

run_picrust_longitudinal <- function(
    comparison_name,
    comparison_info,
    abundance_matrix,
    metadata,
    descriptions,
    subject_col
) {
  
  message("\n-----------------------------------")
  message("Analisi: ", comparison_name)
  message("-----------------------------------")
  
  group_var <- comparison_info$variable
  ref_group <- comparison_info$reference
  comp_group <- comparison_info$comparison
  
  # Seleziona i due gruppi del confronto
  keep <- metadata[[group_var]] %in% c(ref_group, comp_group)
  
  meta_sub <- droplevels(
    metadata[keep, , drop = FALSE]
  )
  
  # Mantieni i campioni presenti nella matrice
  meta_sub <- meta_sub[
    meta_sub$ID %in% colnames(abundance_matrix),
    ,
    drop = FALSE
  ]
  
  # Controlla che siano presenti entrambi i gruppi
  observed_groups <- unique(as.character(meta_sub[[group_var]]))
  
  if (!all(c(ref_group, comp_group) %in% observed_groups)) {
    warning(
      "Confronto saltato: non sono presenti entrambi i gruppi per ",
      comparison_name
    )
    return(NULL)
  }
  
  # Mantieni T0 e T1
  meta_sub <- meta_sub[
    meta_sub$Time %in% c("T0", "T1"),
    ,
    drop = FALSE
  ]
 
  # Mantieni soltanto i partecipanti con un campione a T0 e uno a T1
  paired_subjects <- meta_sub %>%
    group_by(.data[[subject_col]]) %>%
    summarise(
      has_T0 = any(Time == "T0"),
      has_T1 = any(Time == "T1"),
      n_T0 = sum(Time == "T0"),
      n_T1 = sum(Time == "T1"),
      .groups = "drop"
    ) %>%
    filter(has_T0, has_T1, n_T0 == 1, n_T1 == 1) %>%
    pull(.data[[subject_col]])
  
  meta_sub <- meta_sub[
    meta_sub[[subject_col]] %in% paired_subjects,
    ,
    drop = FALSE
  ]
  
  message(
    "Partecipanti con coppia T0-T1: ",
    length(unique(meta_sub[[subject_col]]))
  )

  # Ordina i metadata secondo l'ordine delle colonne della matrice
  meta_sub <- meta_sub[
    match(meta_sub$ID, colnames(abundance_matrix)),
    ,
    drop = FALSE
  ]
  
  # Verifica l'allineamento
  stopifnot(
    identical(meta_sub$ID, colnames(abundance_matrix)[
      match(meta_sub$ID, colnames(abundance_matrix))
    ])
  )
  # Tabella di controllo del disegno
  print(table(meta_sub[[group_var]], meta_sub$Time))
  
  # Escludi campioni con dati mancanti nelle variabili del modello
  complete <- complete.cases(
    meta_sub[, c(group_var, "Time", subject_col), drop = FALSE]
  )
  
  meta_sub <- meta_sub[complete, , drop = FALSE]
  
  # Definisci il gruppo di riferimento
  meta_sub[[group_var]] <- factor(
    as.character(meta_sub[[group_var]]),
    levels = c(ref_group, comp_group)
  )
  
  # Mantieni T0 come riferimento
  meta_sub$Time <- factor(
    as.character(meta_sub$Time),
    levels = c("T0", "T1")
  )
  
  # Imposta il partecipante come fattore
  meta_sub[[subject_col]] <- factor(meta_sub[[subject_col]])
  
  # Controllo delle misure ripetute
  repeated <- meta_sub %>%
    group_by(.data[[subject_col]]) %>%
    summarise(n_time = n_distinct(Time), .groups = "drop")
  
  if (!any(repeated$n_time == 2)) {
    warning(
      "Attenzione: nessun partecipante con entrambi i tempi in ",
      comparison_name
    )
  }
  
  # Matrice pathway x campioni
  mat_sub <- abundance_matrix[
    ,
    meta_sub$ID,
    drop = FALSE
  ]
  
  # MaAsLin2 richiede campioni nelle righe e pathway nelle colonne
  abundance_df <- as.data.frame(
    t(mat_sub),
    check.names = FALSE
  )
  
  # Usa l'ID univoco del campione come nome di riga
  rownames(meta_sub) <- meta_sub$ID
  
  # Verifica che matrice e metadata siano allineati
  stopifnot(
    identical(rownames(abundance_df), rownames(meta_sub))
  )
  
  # Crea la cartella dei risultati
  outdir <- file.path("PICRUSt2_longitudinal", comparison_name)
  
  dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # Salva i dati usati nel modello per tracciabilità
  write.csv(
    meta_sub,
    file = file.path(outdir, "metadata_used.csv"),
    row.names = FALSE
  )
  
  write.table(
    abundance_df,
    file = file.path(outdir, "pathway_abundances_used.tsv"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  
  # Modello:
  # gruppo + tempo + gruppo:tempo
  # con intercetta casuale per partecipante
  # Variabile di interazione: gruppo di confronto a T1
  meta_sub$GroupTime_interaction <- as.numeric(
    meta_sub[[group_var]] == comp_group &
      meta_sub$Time == "T1"
  )
  fit <- Maaslin2(
    input_data = abundance_df,
    input_metadata = meta_sub,
    output = file.path(outdir, "MaAsLin2_results"),
    fixed_effects = c(
      group_var,
      "Time",
      "GroupTime_interaction"
    ),
    random_effects = c(subject_col),
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM",
    max_significance = 1.0,
    min_prevalence = 0.10,
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  )
  
  # Importa la tabella completa dei risultati
  results_file <- file.path(
    outdir,
    "MaAsLin2_results",
    "all_results.tsv"
  )
  
  if (!file.exists(results_file)) {
    warning("Non trovo all_results.tsv per ", comparison_name)
    return(fit)
  }
  
  res <- read.delim(
    results_file,
    check.names = FALSE
  )
  
  # Estrai i coefficienti di interazione gruppo x tempo.
  # I nomi possono essere riportati in uno dei due ordini.
  interaction_results <- res %>%
    filter(metadata == "GroupTime_interaction")
  
  # Aggiungi le descrizioni MetaCyc ai risultati
  descriptions_clean <- descriptions %>%
    transmute(
      pathway_key = str_replace_all(
        trimws(as.character(pathway)),
        "-",
        "."
      ),
      Pathway_description = trimws(as.character(description))
    ) %>%
    distinct(pathway_key, .keep_all = TRUE)
  
  interaction_results <- interaction_results %>%
    mutate(
      feature = trimws(as.character(feature)),
      pathway_key = feature
    ) %>%
    left_join(
      descriptions_clean,
      by = "pathway_key",
      suffix = c(".old", "")
    ) %>%
    mutate(
      Pathway_description = coalesce(
        na_if(Pathway_description, ""),
        feature
      )
    ) %>%
    select(-pathway_key, -any_of("Pathway_description.old"))
  
  # Controlla il risultato
  message(
    "Pathway con descrizione recuperata: ",
    sum(
      interaction_results$Pathway_description !=
        interaction_results$feature,
      na.rm = TRUE
    ),
    " su ",
    nrow(interaction_results)
  )
  
  # Ordina i risultati per significatività
  interaction_results <- interaction_results %>%
    arrange(qval, pval)
  
  # Salva tutti i risultati dell'interazione
  write.csv(
    interaction_results,
    file = file.path(
      outdir,
      paste0(comparison_name, "_interaction_all.csv")
    ),
    row.names = FALSE
  )
  
  # Seleziona quelli significativi con FDR < 0.05
  significant <- interaction_results %>%
    filter(!is.na(qval), qval < 0.05)
  
  write.csv(
    significant,
    file = file.path(
      outdir,
      paste0(comparison_name, "_interaction_FDR05.csv")
    ),
    row.names = FALSE
  )
  
  # Prepara il grafico: primi 20 pathway per q-value.
  # Se non sono significativi, il grafico rimane esplorativo.
  plot_data <- interaction_results %>%
    filter(!is.na(coef), !is.na(qval)) %>%
    arrange(qval, pval) %>%
    slice_head(n = 20)
  
  if (nrow(plot_data) > 0) {
    
    plot_data$Pathway_label <- ifelse(
      is.na(plot_data$Pathway_description) |
        plot_data$Pathway_description == "",
      plot_data$feature,
      plot_data$Pathway_description
    )
    
    # Accorcia le descrizioni molto lunghe
    plot_data$Pathway_label <- ifelse(
      nchar(plot_data$Pathway_label) > 75,
      paste0(
        substr(plot_data$Pathway_label, 1, 72),
        "..."
      ),
      plot_data$Pathway_label
    )
    
    # Ordina per coefficiente per visualizzare meglio il grafico
    plot_data <- plot_data %>%
      arrange(coef)
    
    plot_data$Pathway_label <- factor(
      plot_data$Pathway_label,
      levels = plot_data$Pathway_label
    )
    
    p <- ggplot(
      plot_data,
      aes(
        x = coef,
        y = Pathway_label,
        color = coef > 0
      )
    ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.5,
        color = "grey50"
      ) +
      geom_segment(
        aes(
          x = 0,
          xend = coef,
          y = Pathway_label,
          yend = Pathway_label
        ),
        linewidth = 0.7,
        alpha = 0.75
      ) +
      geom_point(size = 3) +
      scale_color_manual(
        values = c(
          "TRUE" = "#0072B2",
          "FALSE" = "#D55E00"
        ),
        labels = c(
          "TRUE" = "Interaction coefficient > 0",
          "FALSE" = "Interaction coefficient < 0"
        ),
        name = NULL
      ) +
      labs(
        title = gsub("_", " ", comparison_name),
        subtitle = paste(
          "Longitudinal interaction:",
          comp_group, "vs", ref_group,
          "| T1 vs T0"
        ),
        x = "Interaction coefficient (MaAsLin2)",
        y = NULL,
        caption = paste(
          "Top 20 pathways by q-value.",
          "Positive coefficients indicate a larger T0-to-T1 change",
          "in the comparison group relative to the reference group."
        )
      ) +
      theme_classic(base_size = 11) +
      theme(
        plot.title = element_text(face = "plain", size = 13),
        plot.subtitle = element_text(size = 10),
        axis.text.y = element_text(size = 8.5, color = "black"),
        axis.text.x = element_text(color = "black"),
        legend.position = "bottom",
        plot.caption = element_text(size = 8, hjust = 0),
        plot.margin = margin(8, 15, 8, 8)
      )
    
    ggsave(
      filename = file.path(
        outdir,
        paste0(comparison_name, "_top20_interaction.png")
      ),
      plot = p,
      width = 11,
      height = max(5, 0.32 * nrow(plot_data)),
      dpi = 300,
      bg = "white"
    )
    
    print(p)
    
  } else {
    message(
      "Nessun risultato con coefficiente e q-value disponibili per il grafico."
    )
  }
  
  message(
    "Analisi completata: ", comparison_name,
    "\nRisultati: ", outdir,
    "\nPathway con FDR < 0.05: ", nrow(significant)
  )
  
  return(
    list(
      fit = fit,
      all_results = res,
      interaction = interaction_results,
      significant = significant
    )
  )
}

############################
### ESECUZIONE ANALISI
############################

results_list <- list()

for (nm in names(comparisons)) {
  
  results_list[[nm]] <- run_picrust_longitudinal(
    comparison_name = nm,
    comparison_info = comparisons[[nm]],
    abundance_matrix = a,
    metadata = Metadata,
    descriptions = Descriptions,
    subject_col = subject_col
  )
  
}

significant_pathways <- res_interaction %>%
  filter(qval < 0.05) %>%
  arrange(qval)

write.csv(
  significant_pathways,
  "Significant_longitudinal_pathways.csv",
  row.names = FALSE
)
