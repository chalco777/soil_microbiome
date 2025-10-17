---
title: "differential_COG_GO_KEGG"
output: html_document
date: "2024-11-15"
---

## Setup working directory and invoke libraries

We read the eggnog mapper results with the functional assignments for each coding sequence

```{r setup, include=FALSE}
library(knitr)
knitr::opts_knit$set(root.dir = "C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/alen-belem/eggnog_abundance_new/results_differential_2.0")

library(tidyverse)
library(Maaslin2)
library(KEGGREST)
library(httr)
library(jsonlite)
library(rvest)
library(ggpubr)

setwd("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/alen-belem/eggnog_abundance_new")
cog_categories <- c(
  "A" = "RNA processing and modification",
  "B" = "Chromatin structure and dynamics",
  "C" = "Energy production and conversion",
  "D" = "Cell cycle control, cell division,\nchromosome partitioning",
  "E" = "Amino acid transport and metabolism",
  "F" = "Nucleotide transport and metabolism",
  "G" = "Carbohydrate transport and metabolism",
  "H" = "Coenzyme transport and metabolism",
  "I" = "Lipid transport and metabolism",
  "J" = "Translation, ribosomal structure and biogenesis",
  "K" = "Transcription",
  "L" = "Replication, recombination and repair",
  "M" = "Cell wall/membrane/envelope biogenesis",
  "N" = "Cell motility",
  "O" = "Posttranslational modification, protein turnover,\nchaperones",
  "P" = "Inorganic ion transport and metabolism",
  "Q" = "Secondary metabolites biosynthesis, transport\nand catabolism",
  "R" = "General function prediction only",
  "S" = "Function unknown",
  "T" = "Signal transduction mechanisms",
  "U" = "Intracellular trafficking, secretion, and\nvesicular transport",
  "V" = "Defense mechanisms",
  "W" = "Extracellular structures",
  "Y" = "Nuclear structure",
  "Z" = "Cytoskeleton", "EU" = "Amino acid transport and metabolism+\nIntracellular trafficking, secretion, and vesicular transport",
  "KT" = "Transcription+Signal transduction mechanisms",
  "NT" = "Cell motility+Signal transduction mechanisms",
  "CG" = "Energy production and conversion+\nCarbohydrate transport and metabolism",
  "HL" = "Coenzyme transport and metabolism+\nReplication, recombination and repair",
  "IT" = "Lipid transport and metabolism+\nSignal transduction mechanisms",
  "EGP" = "Amino acid + Carbohydrate + Inorganic ion\ntransport and metabolism",
  "MT" = "Cell wall/membrane/envelope biogenesis+\nSignal transduction mechanisms",
  "CEH" = "Energy production and conversion+Amino acid+ \nCoenzyme transport and metabolism",
  "GK" = "Carbohydrate transport and metabolism+Transcription",
  "KP" = "Transcription+Inorganic ion transport and metabolism",
  "IM" = "Lipid transport and metabolism+\nCell wall/membrane/envelope biogenesis",
  "CJ" = "Energy production and conversion+\nTranslation, ribosomal structure and biogenesis",
  "MU" = "Cell wall/membrane/envelope biogenesis+\nIntracellular trafficking, secretion, and vesicular transport",
  "FH" = "Nucleotide transport and metabolism+\nCoenzyme transport and metabolism",
  "CO" = "Energy production and conversion+\nPosttranslational modification, protein turnover, chaperones",
  "EH" = "Amino acid transport and metabolism+\nCoenzyme transport and metabolism",
  "KNT" = "Transcription+Cell motility+\nSignal transduction mechanisms",
  "KLT" = "Transcription+Replication, recombination and repair+\nSignal transduction mechanisms",
  "JM" = "Translation, ribosomal structure and biogenesis+\nCell wall/membrane/envelope biogenesis",
  "EK" = "Amino acid transport and metabolism+Transcription",
  "CI" = "Energy production and conversion+\nLipid transport and metabolism",
  "IK" = "Lipid transport and metabolism+Transcription",
  "JKL" = "Translation, ribosomal structure and biogenesis+\nTranscription+Lipid transport and metabolism",
  "HJ" = "Coenzyme transport and metabolism+\nTranslation, ribosomal structure and biogenesis",
  "EJ" = "Amino acid transport and metabolism+\nTranslation, ribosomal structure and biogenesis",
  "KQ" = "Transcription+Secondary metabolites biosynthesis,\ntransport and catabolism",
  "HK" = "Coenzyme transport and metabolism+Transcription",
  "PQ" = "Inorganic ion transport and metabolism+\nSecondary metabolites biosynthesis, transport and catabolism",
  "CQ" = "Energy production and conversion+\nSecondary metabolites biosynthesis, transport and catabolism",
  "ET" = "Amino acid transport and metabolism+\nSignal transduction mechanisms",
  "GM" = "Carbohydrate transport and metabolism+\nCell wall/membrane/envelope biogenesis",
  "UW" = "Intracellular trafficking, secretion, and vesicular transport+\nExtracellular structures",
  "FP" = "Nucleotide transport and metabolism+\nInorganic ion transport and metabolism",
  "GT" = "Carbohydrate transport and metabolism+\nSignal transduction mechanisms")
df<-read_tsv("combined_gff_final_corrected.tsv",col_names = c("sample","contig","program","type","start","end",                 "score","strand","phase","attributes","mge")) %>% filter(!sample%in%c("barcode03",paste0("barcode",30:33))) 

libra<-read_table("wcline.tsv",col_names = c("size","sample")) %>%
  mutate(sample=gsub("^./(.*)/.*fastq","\\1",sample)) %>% mutate(size=as.numeric(size)) %>% filter(!sample%in%c("barcode03",paste0("barcode",30:33))) %>% mutate(size=size/4)
```

## OBTAINING MAIN TABLES
```{r}
df_main<-df %>%
  # First, separate the basic fields into ID, name, locus_tag, product, rest
  separate(attributes,
           into = c("ID", "name","locus_tag","product","rest"),
           sep = ";",
           extra = "merge") %>% 
  mutate(
    cog    = str_extract(rest, "(?<=em_COG_cat=)[^;]+"),
    GO       = str_extract(rest, "GOs=[^;]+"),
    pathways = str_extract(rest, "KEGG_Pathway=[^;]+"),
    brite    = str_extract(rest, "BRITE=[^;]+"),
    pfams    = str_extract(rest, "PFAMs=[^;]+"),
    conteo   = str_extract(rest, "CONTEO=[^;]+")
  ) %>% 
    mutate(across(c(17:21, 10:13), ~ str_remove(., ".*="))) %>% 
    select(-"rest")
df_main<-df_main %>%filter(!(sample %in% c(paste0("barcode",30:33),"barcode07","barcode11","barcode13","barcode14","barcode18",paste0("barcode",20:22)))) %>% 
  mutate(status = case_when(
    sample %in% c('barcode01', 'barcode02', 'barcode03', 'barcode04', 'barcode05','barcode06') ~ 'barcode30',
    sample %in% c('barcode08', 'barcode09', 'barcode10', 'barcode12') ~ 'barcode32'  ,
    TRUE ~ 'barcode33'
  ))%>% filter(!is.na(conteo))
df_main<-df_main%>% 
  mutate(conteo = replace_na(conteo, "0")) %>% 
  mutate(across(
    c(15:19), 
    ~ replace_na(.x, "Not_assigned")),
    cog = str_replace(cog, "None", "Not_assigned")
  ) %>% mutate(conteo=as.numeric(conteo))
df_main$status <- factor(df_main$status, levels = c("barcode30", "barcode32", "barcode33"))
```

# UTILITY FUNCTIONS

## PARSING FROM API FUNCTIONS

These are used to decorate plot facet titles/heatmaps with human-readable descriptions

get_kegg_description(code) → KEGGREST lookup.
get_go_description_amigo(go_term) → GO API, with AmiGO fallback if a term was replaced by another.
get_brite_description(code) → KEGG list lookup.

```{r}
get_kegg_description <- function(code) {
  tryCatch({
    # Try to get information from KEGG
    kegg_info <- keggGet(code)
    
    # Check if the response contains valid information
    if (length(kegg_info) > 0 && !is.null(kegg_info[[1]]$NAME)) {
      return(kegg_info[[1]]$NAME)
    } else {
      return("Not found")
    }
  }, error = function(e) {
    # Catch errors, including the case of 404 code
    warning(paste("Error fetching KEGG information for code:", code, "-", e$message))
    return("Not found")
  })
}

##TEST
kegg_codes <- c("map00010", "map00020", "map00270")  # Replace with your codes

# Get descriptions of KEGG codes
kegg_descriptions <- sapply(kegg_codes, get_kegg_description)

get_go_description_amigo <- function(go_term) {

# 1) Query the API with the original term
  base_url <- "http://api.geneontology.org/api/ontology/term/"
  query_url <- paste0(base_url, go_term)
  response <- GET(query_url)
  
  if (status_code(response) == 200) {
    content_json <- fromJSON(content(response, "text", encoding = "UTF-8"))
    if (!is.null(content_json$label) && nzchar(content_json$label)) {
      # If we find a valid label, return it
      return(content_json$label)
    }
  }
  
  # 2) No success with the API or there was no 'label'.
  #    Scrape AmiGO to see if there is a "replaced by".
  
  web_url <- paste0("https://amigo.geneontology.org/amigo/term/", go_term)
  page <- tryCatch(read_html(web_url), error = function(e) NULL)
  
  # If we can't access the page or page is NULL, return "Not found".
  if (is.null(page)) {
    return("Not found")
  }
  replaced_by_term <- page %>%
    html_nodes(xpath = "//dd[strong[text()='replaced by']]") %>%
    html_node("a") %>%
    html_text(trim = TRUE)
  
  # If nothing, return "Not found".
  if (is.na(replaced_by_term) || is.null(replaced_by_term) || replaced_by_term == "") {
    return("Not found")
  }
  
  # 3) If we find a replacement term, query the API with that new term.
  new_query_url <- paste0(base_url, replaced_by_term)
  new_response <- GET(new_query_url)
  
  if (status_code(new_response) == 200) {
    new_content <- fromJSON(content(new_response, "text", encoding = "UTF-8"))
    if (!is.null(new_content$label) && nzchar(new_content$label)) {
      return(new_content$label)
    }
  }
  
  # 4) If we get here, nothing useful was found. Return "Not found".
  return("Not found")
}
##TEST
go_terms <- c("GO:0008150", "GO:0000004", "GO:0017144"

)

get_brite_description <- function(brite_code) {
  # Check if the code is among the 'names()' of the brite_data vector
  brite_data <- keggList("brite")
  if (brite_code %in% names(brite_data)) {
    return(brite_data[[brite_code]])
  } else {
    return(paste("The code", brite_code, "was not found in brite_data."))
  }
}
brite_codes <- c("br08303", "ko01000", "ko01004","ko03110" )  # Examples

brite_descriptions <- sapply(as.character(brite_codes), get_brite_description)
brite_descriptions
kegg_descriptions

```

## RECODE FUNCTIONS

```{r}
recode_feature_names <- function(features, annotation_col, max_width = 300) {
  
  # Internal function to add line breaks
  wrap_text <- function(txt, width = 300) {
    paste(strwrap(txt, width = width), collapse = "\n")
  }
  
  recoded <- sapply(features, function(ft) {
    
    # For 'cog' or 'pfam', return the feature unchanged
    if (tolower(annotation_col) %in% c("cog", "pfams")) {
      return(ft)
    }
    
    # Otherwise, try to recode
    if (tolower(annotation_col) == "pathways") {
      desc <- get_kegg_description(ft)
      print("IN KEGGGGGGGG IN THE RECODE FUNCTION")
      print(desc)
    } else if (tolower(annotation_col) == "brite") {
      desc <- get_brite_description(ft)
    } else if (tolower(annotation_col) == "go") {
      desc <- get_go_description_amigo(ft)  # Example
    } else {
      desc <- "No Description"
      print("NO DESCRIPTION INN RECODE FUNCTIONN")
    }
    
    # Adjust the description if it is empty
    if (is.null(desc) || desc == "") desc <- "No Description"
    
    # Generate the final string
    out <- paste0(ft, ": ", desc)
    wrap_text(out, width = max_width)
    
  }, USE.NAMES = FALSE)
  
  # Return a named vector (original names -> recoded)
  names(recoded) <- features
  return(recoded)
}

```

### FUNCTION ANALYSIS II

This is a full pipeline function: pivots to wide CPM, builds metadata, runs MaAsLin2, ranks features by KW p-value, picks a slice (e.g., rows 1–9), recodes labels (COG/GO/KEGG), draws & saves a faceted boxplot PNG, and returns a list of useful objects (wide matrix, metadata, results, KW tables, filtered data, and the ggplot object).
```{r}
analysis_maaslin2_boxplot <- function(
  input_data,          # Data frame with columns: sample, status, conteo_total and the annotation column
  annotation_col,      # Column name: "GO", "KEGG", "BRITE", "PFAM", etc.
  origin = c("virus", "plasmid", "chromosome"),  # For naming files and/or plots
  slice_range = c(1, 9), # Row range to take after sorting by p-value (e.g. c(1,9) or c(10,18))
  vector_recode = NULL # Named vector to recode categories
) {
  # --- Internal function to format the title of each facet --- #
  # Combines the "code" (e.g. GO:XXXX) with the "description" and,
  # if it exceeds 25 characters, inserts line breaks.
  wrap_title <- function(code, desc, max_width = 25) {
    if (is.null(desc) || desc == "") desc <- "No Description"
    # Add ": " between the code and the description
    combined <- paste0(code, ": ", desc)
    # Adjust with line breaks if it exceeds max_width
    paste(strwrap(combined, width = max_width), collapse = "\n")
  }

  # Titles for "origin" (virus, plasmid, chromosome)
  origin_titles <- c(
    virus = "Virus",
    plasmid = "Plasmids",
    chromosome = "Chromosomal Contigs"
  )

  # Titles for each annotation
  annotation_titles <- c(
    cog = "Clusters of Orthologous Groups (COGs)",
    GO = "Gene Ontology (GO) Classes",
    pathways = "KEGG Pathways",
    brite = "BRITE Hierarchies",
    pfams = "PFAM Domains"
  )

  origin <- match.arg(origin)  # Validate 'origin' argument
  annotation_col_sym <- rlang::sym(annotation_col)  # Convert the column name to symbol

  # 1. pivot_wider to create the "cpm" matrix
  df_cpm <- input_data %>%
    tidyr::pivot_wider(
      names_from  = !!annotation_col_sym,   # Column with GO/KEGG/BRITE/PFAM
      values_from = conteo_total,
      values_fill = 0
    ) %>% 
    dplyr::select(-1) %>% 
    as.data.frame()

  # Assume the first column is "sample" and use it as rownames
  rownames(df_cpm) <- df_cpm[, 1]
  df_cpm <- df_cpm[, -1]

  # 2. Create the metadata
  metadata_df_cpm <- input_data %>%
    dplyr::select(sample, status) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      name = dplyr::case_when(
        status == "barcode30" ~ "Baseline",
        status == "barcode32" ~ "24 hours",
        status == "barcode33" ~ "10 days",
        TRUE                  ~ status
      )
    ) %>%
    as.data.frame()

  rownames(metadata_df_cpm) <- metadata_df_cpm[, 1]
  metadata_df_cpm <- metadata_df_cpm[, -1]

  # 3. Run Maaslin2 with dynamic output folder
  output_folder <- paste0("Maaslin2_", origin, "_", annotation_col, "_cpm")

  fit_data <- Maaslin2(
    input_data      = df_cpm,
    input_metadata  = metadata_df_cpm,
    output          = output_folder,
    fixed_effects   = c("name"),
    reference       = "name,Baseline",
    min_prevalence  = 0.2,
    normalization   = "TSS",
    transform       = "LOG",
    analysis_method = "LM"
  )

  # 4. Kruskal-Wallis test for each category
  group_by_string <- annotation_col
  k_cpm <- input_data %>%
    dplyr::group_by(.data[[group_by_string]]) %>%
    dplyr::filter(dplyr::n_distinct(status) >= 2, dplyr::n_distinct(sample) >= 4) %>%
    dplyr::summarise(
      k          = stats::kruskal.test(conteo_total ~ status)$p.value,
      wilcoxon_p = list(
        pairwise.wilcox.test(
          conteo_total, status, 
          exact = FALSE,
          p.adjust.method = "BH"
        )$p.value
      )
    ) %>%
    dplyr::arrange(k)

  # 5. Select the subset of rows by range
  start_slice <- slice_range[1]
  end_slice   <- slice_range[2]
  k_subset <- k_cpm %>%
    dplyr::slice(start_slice:end_slice)

  # 6. Join to the original dataframe
  filtered_cpm <- input_data %>%
    dplyr::filter(.data[[group_by_string]] %in% k_subset[[group_by_string]]) %>%
    dplyr::left_join(k_subset, by = group_by_string) %>%
    dplyr::mutate(
      !!group_by_string := factor(
        .data[[group_by_string]],
        levels = k_subset[[group_by_string]]
      )
    )

  # 7. Recoding for boxplot titles
  #    Add code + description with line break and
  #    for GO, force the first letter to uppercase.
  if (tolower(annotation_col) == "go") {
    all_go_codes <- unique(filtered_cpm[[group_by_string]])
    
    go_wrapped <- sapply(all_go_codes, function(go_code) {
      desc <- get_go_description_amigo(go_code)
      # First letter uppercase
      if (!is.null(desc) && nzchar(desc)) {
        desc <- paste0(toupper(substr(desc, 1, 1)), substr(desc, 2,nchar(desc)))
      } else {
        desc <- "No Description"
      }
      wrap_title(go_code, desc, max_width = 25)
    })
    
    names(go_wrapped) <- all_go_codes
    
    filtered_cpm <- filtered_cpm %>%
      dplyr::mutate(
        !!group_by_string := dplyr::recode(
          .data[[group_by_string]],
          !!!go_wrapped
        )
      )
    
  } else if (tolower(annotation_col) == "pathways") {
    all_kegg_codes <- unique(filtered_cpm[[group_by_string]])
    
    kegg_wrapped <- sapply(all_kegg_codes, function(kegg_code) {
      desc <- get_kegg_description(kegg_code)
      wrap_title(kegg_code, desc, max_width = 25)
    })
    names(kegg_wrapped) <- all_kegg_codes
    filtered_cpm <- filtered_cpm %>%
      dplyr::mutate(
        !!group_by_string := dplyr::recode(
          .data[[group_by_string]],
          !!!kegg_wrapped
        )
      )
    
  } else if (tolower(annotation_col) == "cog" && !is.null(vector_recode)) {
    # For COG it is assumed that 'vector_recode' already contains the recoding
    filtered_cpm <- filtered_cpm %>%
      dplyr::mutate(
        !!group_by_string := dplyr::recode(
          .data[[group_by_string]],
          !!!vector_recode
        )
      )
  }

  # 8. Create the boxplot for this subset
  origin_title <- origin_titles[[origin]]
  annotation_title <- annotation_titles[[annotation_col]]
  combined_title <- paste0(annotation_title, " Abundance in ", origin_title)

  plot_box <- filtered_cpm %>%
    ggplot2::ggplot(aes(x = status, y = conteo_total, fill = status)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.5) +
    ggplot2::facet_wrap(
      stats::as.formula(paste0("~ ", group_by_string)),
      scales = "free_y"
    ) +
    ggplot2::geom_text(
      aes(label = paste("p =", round(k, 4))),
      x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
      size = 2.5, color = "black"
    ) +
    ggplot2::labs(
      title = combined_title, x = "Status", y = "CPM"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_discrete(labels = function(x) dplyr::case_when(
      x == "barcode30" ~ "Baseline",
      x == "barcode32" ~ "24 hours",
      x == "barcode33" ~ "10 days",
      TRUE             ~ x
    )) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title.y     = ggplot2::element_text(size = 10),
      axis.text.x      = ggplot2::element_text(size = 9),
      axis.text.y      = ggplot2::element_text(size = 8),
      strip.text       = ggplot2::element_text(size = 7, color = "black"),
      panel.grid.major = ggplot2::element_line(color = "grey90", size = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.placement  = "outside"
    ) +
    ggplot2::guides(fill = "none")

  filename <- paste0(
    tolower(annotation_col), "_",
    tolower(origin), "_cpm_",
    start_slice, "-", end_slice, ".png"
  )

  ggplot2::ggsave(
    filename = filename,
    plot     = plot_box,
    width    = 8,
    height   = 6,
    dpi      = 300,
    bg       = "white"
  )

  # 9. Return everything needed
  list(
    wide_data       = df_cpm,         # Resulting pivot matrix
    metadata        = metadata_df_cpm,
    maaslin_results = fit_data,       # Maaslin2
    kruskal_full    = k_cpm,          # Full table with p-values
    kruskal_subset  = k_subset,       # Selected subset
    filtered_data   = filtered_cpm,   # Data frame for the boxplot
    boxplot         = plot_box        # ggplot2 object
  )
}

```

### FUNCTION ANALYSIS III

This functions is richer than function analysis 2 and overrides it.

- Uses an expanded_data long table internally (safer for KW across all features).
- Adds ggpubr::stat_compare_means to draw global KW p-values and pairwise Wilcoxon asterisks on the boxplots.
- Writes output to a different folder suffix (_cpm_test5) and filename suffix (test5.png).
-If MaAsLin2 finds significant features (q < 0.25), it builds two heatmaps:

  “significance” heatmap using −log10(q)·sign(coef), after recoding feature names with descriptions;

  mean CPM-based log2 fold change heatmap vs. Baseline (also recoded).

- Returns, besides what vII returns, heat_mat, expanded_data, and expanded_filt for downstream reuse

```{r}
analysis_maaslin2_boxplot <- function(
  input_data,          # Data frame with columns: sample, status, conteo_total and the annotation column
  annotation_col,      # Column name: "GO", "KEGG", "BRITE", "PFAM", etc.
  origin = c("virus", "plasmid", "chromosome"),  # For naming files and/or plots
  slice_range = c(1, 9), # Row range to take after sorting by p-value (e.g. c(1,9) or c(10,18))
  vector_recode = NULL # Named vector to recode categories (e.g. for COG)
) {
  # --- Internal function to format the title of each facet --- #
  wrap_title <- function(code, desc, max_width = 25) {
    if (is.null(desc) || desc == "") desc <- "No Description"
    # Add ": " between the code and the description
    combined <- paste0(code, ": ", desc)
    # Adjust with line breaks if it exceeds max_width
    paste(strwrap(combined, width = max_width), collapse = "\n")
  }
  # Titles for "origin" (virus, plasmid, chromosome)
  origin_titles <- c(
    virus     = "Virus",
    plasmid   = "Plasmids",
    chromosome = "Chromosomal Contigs"
  )

  # Titles for each annotation
  annotation_titles <- c(
    cog      = "Clusters of Orthologous Groups (COGs)",
    GO       = "Gene Ontology (GO) Classes",
    pathways = "KEGG Pathways",
    brite    = "BRITE Hierarchies",
    pfams    = "PFAM Domains"
  )

  origin <- match.arg(origin)  # Validate 'origin' argument
  annotation_col_sym <- rlang::sym(annotation_col)  # Convert the column name to symbol

  # 1. pivot_wider to create the "cpm" matrix
  df_cpm <- input_data %>%
    tidyr::pivot_wider(
      names_from  = !!annotation_col_sym,   # Column with GO/KEGG/BRITE/PFAM/COG
      values_from = conteo_total,
      values_fill = 0
    ) %>% 
        dplyr::select(-1) %>% 
    as.data.frame()

  # Assume the first column is "sample" and use it as rownames
  rownames(df_cpm) <- df_cpm[, 1]
  df_cpm <- df_cpm[, -1]

  # 2. Create the metadata
  metadata_df_cpm <- input_data %>%
    dplyr::select(sample, status) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      name = dplyr::case_when(
        status == "barcode30" ~ "Baseline",
        status == "barcode32" ~ "24 hours",
        status == "barcode33" ~ "10 days",
        TRUE                  ~ status
      )
    ) %>%
    as.data.frame()

  rownames(metadata_df_cpm) <- metadata_df_cpm[, 1]
  metadata_df_cpm <- metadata_df_cpm[, -1]

  # 3. Run Maaslin2 with dynamic output folder
  output_folder <- paste0("Maaslin2_", origin, "_", annotation_col, "_cpm_test5")

  fit_data <- Maaslin2(
    input_data      = df_cpm,
    input_metadata  = metadata_df_cpm,
    output          = output_folder,
    fixed_effects   = c("name"),
    reference       = "name,Baseline",
    min_prevalence  = 0.2,
    normalization   = "TSS",
    transform       = "LOG",
    analysis_method = "LM"
  )

  # 4. Kruskal-Wallis test for each category
  group_by_string <- annotation_col
  
  expanded_data <- input_data %>%
  tidyr::pivot_wider(
    names_from  = !!annotation_col_sym,   # Column with GO/KEGG/BRITE/PFAM/COG
    values_from = conteo_total,
    values_fill = 0
  ) %>%
  tidyr::pivot_longer(
    cols      = -c(sample, status),  # Keep 'sample' and 'status', the rest becomes "long"
    names_to  = annotation_col,      # Reuse the same variable
    values_to = "conteo_total"
  )
  
k_cpm <- expanded_data %>%
  dplyr::group_by(.data[[annotation_col]]) %>%
  #dplyr::filter(dplyr::n_distinct(status) >= 2, dplyr::n_distinct(sample) >= 4) %>%
  dplyr::summarise(
    k          = stats::kruskal.test(conteo_total ~ status)$p.value,
    wilcoxon_p = list(
      pairwise.wilcox.test(
        conteo_total, status, 
        exact = FALSE,
        p.adjust.method = "BH"
      )$p.value
    )
  ) %>%
  dplyr::arrange(k)

  # 5. Select the subset of rows by range
  start_slice <- slice_range[1]
  end_slice   <- slice_range[2]
  k_subset <- k_cpm %>%
    dplyr::slice(start_slice:end_slice)

  # 6. Join to the original dataframe with p-values
filtered_cpm <- expanded_data %>%
  dplyr::filter(.data[[annotation_col]] %in% k_subset[[annotation_col]]) %>%
  dplyr::left_join(k_subset, by = annotation_col) %>%
  dplyr::mutate(
    # Adjust the order of the facet according to k_subset
    !!annotation_col := factor(
      .data[[annotation_col]],
      levels = k_subset[[annotation_col]]
    )
  )

  # 7. Recoding for boxplot titles
  if (tolower(annotation_col) == "go") {
    all_go_codes <- unique(filtered_cpm[[group_by_string]])
    
    go_wrapped <- sapply(all_go_codes, function(go_code) {
      desc <- get_go_description_amigo(go_code)
      # First letter uppercase
      if (!is.null(desc) && nzchar(desc)) {
        desc <- paste0(toupper(substr(desc, 1, 1)), substr(desc, 2,nchar(desc)))
      } else {
        desc <- "No Description"
      }
      wrap_title(go_code, desc, max_width = 25)
    })
    
    names(go_wrapped) <- all_go_codes
    
    filtered_cpm <- filtered_cpm %>%
      dplyr::mutate(
        !!group_by_string := dplyr::recode(
          .data[[group_by_string]],
          !!!go_wrapped
        )
      )
    
  } else if (tolower(annotation_col) == "pathways") {
    all_kegg_codes <- unique(filtered_cpm[[group_by_string]])
    kegg_wrapped <- sapply(all_kegg_codes, function(kegg_code) {
      desc <- get_kegg_description(kegg_code)
      wrap_title(kegg_code, desc, max_width = 25)
    })
    names(kegg_wrapped) <- all_kegg_codes
    filtered_cpm <- filtered_cpm %>%
      dplyr::mutate(
        !!group_by_string := dplyr::recode(
          .data[[group_by_string]],
          !!!kegg_wrapped
        ))
  } else if (tolower(annotation_col) == "brite") {
    all_brite_codes <- unique(filtered_cpm[[group_by_string]])
    brite_wrapped <- sapply(as.character(all_brite_codes), function(brite_code) {
      desc <- get_brite_description(brite_code)   # NEW FUNCTION
      wrap_title(brite_code, desc, max_width = 25)
    })
    gg<-brite_wrapped
    names(brite_wrapped) <- all_brite_codes
    filtered_cpm <- filtered_cpm %>%
      dplyr::mutate(
        !!group_by_string := dplyr::recode(
          .data[[group_by_string]],
          !!!brite_wrapped
        )
      )

  } else if (tolower(annotation_col) == "cog" && !is.null(vector_recode)) {
    filtered_cpm <- filtered_cpm %>%
      dplyr::mutate(
        !!group_by_string := dplyr::recode(
          .data[[group_by_string]],
          !!!vector_recode
        )
      )
  }
  # PFAM -> No recoding

  # 8. Create the boxplot for this subset
  origin_title <- origin_titles[[origin]]
  annotation_title <- annotation_titles[[annotation_col]]
  combined_title <- paste0(annotation_title, " Abundance in ", origin_title)

  # To add KW asterisks, generate column with 'k_stars'
  plot_box <- filtered_cpm %>% 
  #     dplyr::group_by(!!rlang::sym(group_by_string),status) %>%
  # dplyr::filter(sd(conteo_total) > 0) %>% 
  # dplyr::ungroup()%>%
    ggplot2::ggplot(aes(x = status, y = conteo_total, fill = status)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.5) +
    ggplot2::facet_wrap(
      stats::as.formula(paste0("~ ", group_by_string)),
      scales = "free_y"
    ) +
    # Show KW p-value and stars in the upper right corner +
    #    # ggplot2::geom_text(
    # 
    #   aes(label = paste("p =", round(k, 4))),
    # 
    #   x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
    # 
    #   size = 2.5, color = "black"
    # 
    # ) +
    ggplot2::labs(
      title = combined_title, x = "Status", y = "CPM"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_discrete(labels = function(x) dplyr::case_when(
      x == "barcode30" ~ "Baseline",
      x == "barcode32" ~ "24 hours",
      x == "barcode33" ~ "10 days",
      TRUE             ~ x
    )) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title.y     = ggplot2::element_text(size = 10),
      axis.text.x      = ggplot2::element_text(size = 9),
      axis.text.y      = ggplot2::element_text(size = 8),
      strip.text       = ggplot2::element_text(size = 7, color = "black"),
      panel.grid.major = ggplot2::element_line(color = "grey90", size = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.placement  = "outside"
    ) +
    ggplot2::guides(fill = "none")
  
  #  # 8.1. Add statistics with ggpubr::stat_compare_means
  # # Define pairwise comparisons for Wilcoxon.
  # # Adjust according to your real groups in 'status'.
  comparisons_list <- list(
    c("barcode30", "barcode32"),  # Baseline vs 24h
    c("barcode30", "barcode33"),  # Baseline vs 10d
    c("barcode32", "barcode33")   # 24h vs 10d
  )

# (B) Add Kruskal-Wallis p-value (global)
plot_box <- plot_box +
  stat_compare_means(
    aes(group = status),
    method = "kruskal.test",
    label = "p.format",    # Shows *** if p<0.001, etc. and ns if not
    hide.ns = TRUE,        # Hides "ns" and its bracket if not significant
    label.x.npc = "right",  # Adjust X position
    label.y.npc = "top",   # Adjust Y position
    size = 2.5,
    vjust = 1.2,
  )

# (C) Add pairwise comparisons (Wilcoxon)
plot_box<-plot_box +
  stat_compare_means(
    aes(group = status),
    method = "wilcox.test",
    p.adjust.method = "BH",
    comparisons = comparisons_list,
    label = "p.signif",    
    hide.ns = TRUE,        # Hides brackets and label if p>=0.05
    step.increase = 0.1,   # Additional space between successive brackets
    bracket.size = 0.5,    
    tip.length = 0.02,     
    size = 4,              
    geom = "text",         # Text with asterisks
    vjust = 0.5
  )

  filename <- paste0(
    tolower(annotation_col), "_",
    tolower(origin), "_cpm_",
    start_slice, "-", end_slice, "test5.png"
  )

  ggplot2::ggsave(
    filename = filename,
    plot     = plot_box,
    width    = 8,
    height   = 6,
    dpi      = 300,
    bg       = "white"
  )
  print(head(expanded_data))

  res_sign <- fit_data$results %>%
  dplyr::filter(qval < 0.25)
  # --- HEATMAP CREATION ---
if (nrow(res_sign) == 0) {
  message("No significant features (qval < 0.25) found. Skipping heatmaps.")
} else {
  # ---------------------------
  # A) SIGNIFICANCE HEATMAP (-log10(qval)*sign(coef))
  # ---------------------------
  res_sign <- res_sign %>%
    dplyr::mutate(
      qval_adj   = dplyr::if_else(qval == 0, 1e-300, qval),
      heat_value = -log10(qval_adj) * sign(coef)
    )
  heat_data <- res_sign %>%
    dplyr::select(feature, value, heat_value) %>%
    tidyr::pivot_wider(
      names_from  = value, 
      values_from = heat_value
    )
  
  heat_mat <- as.data.frame(heat_data[,-1])
  rownames(heat_mat) <- heat_data$feature
  print(heat_mat)
  # 2A. Recode (in case there is a KEGG, BRITE, GO, etc. description)
       # "kegg", "brite", "go", etc.
  sig_features    <- rownames(heat_mat)  # significant features
  recoded_vec     <-recode_feature_names(sig_features,group_by_string)
  print("BEFORE")
  print(sig_features)
  print("AFTER")
  print(recoded_vec)
  print("ADDING TO HEAT_MAT")
  rownames(heat_mat) <- recoded_vec
  print(heat_mat)
  class(heat_mat)
# 1. Open the device (PNG, PDF, etc.)
heat_mat$feature <- rownames(heat_mat)
df_long <- heat_mat %>%
  pivot_longer(
    cols      = -feature,
    names_to  = "status",
    values_to = "value"
  )

p <- ggplot(df_long, aes(x = status, y = feature, fill = value)) +
  geom_tile(color = "black") +  # Black border on each cell
  scale_fill_gradient2(
    midpoint = 0,   # Center the scale at 0
    low = "blue", mid = "white", high = "red",
    name = "log2FC"
  ) +
  labs(
    title = "Heatmap of Log-Scaled Significance of Features with Baseline as Reference" ,
    x = "Condition",
    y = "Feature"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    # Center and format the title
    plot.title = element_text(
      hjust  = 0.5,  # Center the title horizontally
      face   = "bold",
      size   = 14
    ),
    # Adjust the size of the X axis (columns) and Y axis (rows) labels
    axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x =element_text(size = 14),
    axis.title.y=element_text(size = 14),
    # Remove grid lines
    panel.grid  = element_blank()
  ) +
  # Reduce empty space around cells
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  # Keep fixed aspect ratio (makes cells square).
  # If you want rectangular, adjust 'ratio'
  coord_fixed()

# Now save with a smaller size so the cells are smaller
ggsave(
  filename = paste0("heatmap_significance_", annotation_col, ".png"),
  plot     = p,
  width    = 8,   # Adjust width (in inches if units not specified)
  height   = 6,   # Adjust height
  dpi      = 300,
  bg       = "white"
)
  
  # ---------------------------
  # B) HEATMAP OF LOG2(ABUNDANCE)
  # ---------------------------
  # # 3A. Filter 'expanded_data' only with those significant features
expanded_filt <- gt$expanded_data %>%
  dplyr::filter(!!annotation_col_sym%in% sig_features) %>%
  mutate(name=ifelse(status == "barcode30", "Baseline", 
                     ifelse(status == "barcode32", "24 hours",
                            ifelse(status=="barcode33", "10 days",status)))) %>% 
  dplyr::filter(name%in% c(res_sign$value,"Baseline"))  %>% group_by(!!annotation_col_sym,name) %>%
  summarize(mean_value=mean(conteo_total)) %>% 
  pivot_wider(names_from=name,values_from = mean_value) %>%  ungroup() %>% 
  mutate(across(-c(!!annotation_col_sym,Baseline),~ log2( (.x + 1) / (Baseline + 1) ), .names = "{.col}")) 
  old_feats_2      = unique(expanded_filt %>% dplyr::pull(!!annotation_col_sym))
  print(expanded_filt)
  print("BEFORE")
  print(old_feats_2)
  recoded_vec_2 <- recode_feature_names(old_feats_2,group_by_string)
  # 
  print("AFTER")
  print(recoded_vec_2)
  expanded_filt <- expanded_filt %>%
    dplyr::mutate(
      !!group_by_string := dplyr::recode(
        .data[[group_by_string]],
        !!!recoded_vec_2
      )
    )
  expanded_filt<-expanded_filt %>%
  column_to_rownames(var = rlang::as_string(annotation_col_sym)) %>% select(-Baseline)
  print(head(expanded_filt))
expanded_filt$feature <- rownames(expanded_filt)
df_long <- expanded_filt %>%
  pivot_longer(
    cols      = -feature,
    names_to  = "status",
    values_to = "value"
  )

p <- ggplot(df_long, aes(x = status, y = feature, fill = value)) +
  geom_tile(color = "black") +  # Black border on each cell
  scale_fill_gradient2(
    midpoint = 0,   # Center the scale at 0
    low = "blue", mid = "white", high = "red",
    name = "log2FC"
  ) +
  labs(
    title = "Heatmap of log2(Abundance) of Features (Qval < 0.25) with Baseline as Reference",
    x = "Condition",
    y = "Feature"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    # Center and format the title
    plot.title = element_text(
      hjust  = 0.5,  # Center the title horizontally
      face   = "bold",
      size   = 14
    ),
    # Adjust the size of the X axis (columns) and Y axis (rows) labels
    axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x =element_text(size = 14),
    axis.title.y=element_text(size = 14),
    # Remove grid lines
    panel.grid  = element_blank()
  ) +
  # Reduce empty space around cells
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  # Keep fixed aspect ratio (makes cells square).
  # If you want rectangular, adjust 'ratio'
  coord_fixed()

# Now save with a smaller size so the cells are smaller
ggsave(
  filename = paste0("heatmap_log2abundance_", annotation_col, ".png"),
  plot     = p,
  width    = 8,   # Adjust width (in inches if units not specified)
  height   = 6,   # Adjust height
  dpi      = 300,
  bg       = "white"
)





  
}
  # 9. Return everything needed
  list(
    wide_data       = df_cpm,         # Resulting pivot matrix
    metadata        = metadata_df_cpm,
    maaslin_results = fit_data,       # Maaslin2
    kruskal_full    = k_cpm,          # Full table with p-values
    kruskal_subset  = k_subset,       # Selected subset
    filtered_data   = filtered_cpm,   # Data frame for the boxplot
    boxplot         = plot_box,        # ggplot2 object
  #log2_mat=log2_mat,
  heat_mat=heat_mat,
  expanded_data=expanded_data,
  expanded_filt=expanded_filt
    )  
}

```



# Clusters of Orthologous Genes analysis
MAASLIN2 DIFFERENTIAL ANALYSIS AND BOXPLOTS

Raw count of reads assigned to plasmid CDS

```{r}
cat_plasmid<-df_main %>% 
    filter(str_detect(mge, "PLASMID")) %>% filter(cog!="Not_assigned") %>% 
  group_by(status, sample, cog) %>% summarise(conteo_total=sum(conteo)) %>% ungroup()

plasmid_df <- cat_plasmid %>% 
  pivot_wider(names_from = cog, values_from = conteo_total, values_fill = 0) %>% select(-1) %>% as.data.frame() 
rownames(plasmid_df)<-plasmid_df[,1]
plasmid_df<-plasmid_df[,-1]

# Create the metadata dataframe
metadata <- cat_plasmid %>%
  select(sample, status) %>% distinct() %>% 
  mutate(name=ifelse(status == "barcode30", "Baseline", 
      ifelse(status == "barcode32", "24 hours",
      ifelse(status=="barcode33", "10 days",status)))) %>% as.data.frame()
rownames(metadata)<-metadata[,1]
metadata<-metadata[,-1]

fit_data=Maaslin2(input_data = plasmid_df,
                  input_metadata = metadata,
                      output = "Maaslin2_plasmid_cog_none",       # Output folder for results
                  fixed_effects = c("name"), 
                  reference = "name,Baseline",
                  min_prevalence = 0.2,
                  normalization="NONE",
                  transform = "NONE",
                  analysis_method = "NEGBIN")

# FILTER BY LOWEST KRUSKAL WALLIS P-VALUE (data is not normal according to shapiro test)
k<-cat_plasmid %>% group_by(cog) %>% filter(n_distinct(status) >= 2) %>% summarise(k=kruskal.test(conteo_total ~ status)$p.value,                                                  wilcoxon_p=list(pairwise.wilcox.test(conteo_total,status, exact=FALSE, p.adjust.method="BH")$p.value))  %>% arrange(k) %>% 
  slice_head(n = 9) %>%
  pull(cog)

test<-cat_plasmid %>% filter(cog=="C")
p<-pairwise.wilcox.test(test$conteo_total,test$status, exact=FALSE, p.adjust.method = "BH")

# Filter the dataframe for only the 9 categories with the lowest p value in Kruskal
filtered_cat_plasmid <- cat_plasmid %>%
  filter(cog %in% k)

# Create the boxplot graphs
plots <- filtered_cat_plasmid %>%
  ggplot(aes(x = status, y = conteo_total, fill = status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ cog, scales = "free_y") +
  labs(title = "Boxplots by Status for the 9 Categories with Highest Median",
       x = "Status",
       y = "Total Count") +
  theme_minimal()+
    scale_x_discrete(labels = function(x) ifelse(x == "barcode30", "Baseline", 
                                               ifelse(x == "barcode32", "24 hours",
                                                      ifelse(x=="barcode33", "10 days",x))))
```

### COG PLASMIDS CPM 

```{r}
cat_plasmid_cpm <- df_main  %>% filter(str_detect(mge, "PLASMID"))%>% filter(cog!="Not_assigned") %>% inner_join(libra, by="sample") %>%mutate(cpm=(conteo/size)*1e6) %>% 
  group_by(status, sample, cog) %>% summarise(conteo_total=sum(cpm)) %>%ungroup()

## Perform Maaslin2 analysis

analysis_maaslin2_boxplot(input_data = cat_plasmid_cpm,"cog","plasmid",vector_recode =cog_categories)

```

### COG VIRUS CPM 
```{r}
cat_virus_cpm <- df_main  %>% filter(str_detect(mge, "VIRUS")) %>% filter(cog!="Not_assigned")%>% inner_join(libra, by="sample") %>% separate(mge, into = c("type_mge","length","rango"), sep=";",fill="right") %>% #now keep genes within provirus
  mutate(
    in_range = case_when(
      str_detect(type_mge, "PROVIRUS") ~ if_else(!is.na(rango), {
                  rango_parts <- str_split(rango, "-", simplify = TRUE)
                  rango_start <- as.numeric(rango_parts[, 1])
                  rango_end <- as.numeric(rango_parts[, 2])
                  start >= rango_start & end <= rango_end
                }, FALSE,FALSE #else is false
                ),TRUE ~ TRUE)
  ) %>%
  filter(in_range)%>%mutate(cpm=(conteo/size)*1e6) %>% 
  group_by(status, sample, cog) %>% summarise(conteo_total=sum(cpm)) %>%ungroup()
analysis_maaslin2_boxplot(input_data = cat_virus_cpm,"cog","virus", vector_recode = cog_categories)

```

### COG CHROMOSOMAL CONTIGS (NOT MGE)
```{r}
cat_main_cpm <- df_main  %>% filter(str_detect(mge, "NOT_MGE"))%>% filter(cog!="Not_assigned") %>% inner_join(libra, by="sample") %>%mutate(cpm=(conteo/size)*1e6) %>% 
  group_by(status, sample, cog) %>% summarise(conteo_total=sum(cpm)) %>%ungroup()
analysis_maaslin2_boxplot(input_data = cat_main_cpm,"cog","chromosome", vector_recode = cog_categories)

```

# Gene Ontology terms analysis
 MAASLIN2 AND BOXPLOTS 

### GO RAW COUNT
```{r}
### PLASMID RAW COUNT
cat_plasmid_go<-df_main %>% 
    filter(str_detect(mge, "PLASMID"),GO!="Not_assigned") %>% select(status,sample,GO, conteo)%>%separate_rows(GO, sep = ",")%>% 
  group_by(status, sample, GO) %>% summarise(conteo_total=sum(conteo)) %>% ungroup()  

plasmid_df <- cat_plasmid_go %>% 
  pivot_wider(names_from = GO, values_from = conteo_total, values_fill = 0) %>% select(-1) %>% as.data.frame()
rownames(plasmid_df)<-plasmid_df[,1]
plasmid_df<-plasmid_df[,-1]

# Create the metadata dataframe
metadata <- cat_plasmid_go %>%
  select(sample, status) %>% distinct() %>% 
  mutate(name=ifelse(status == "barcode30", "Baseline", 
      ifelse(status == "barcode32", "24 hours",
      ifelse(status=="barcode33", "10 days",status)))) %>% as.data.frame()
rownames(metadata)<-metadata[,1]
metadata<-metadata[,-1]

fit_data=Maaslin2(input_data = plasmid_df,
                  input_metadata = metadata,
                      output = "Maaslin2_plasmid_GO",       # Output folder for results
                  fixed_effects = c("name"), 
                  reference = "name,Baseline",
                  min_prevalence = 0.2,
                  normalization="TSS",
                  transform = "LOG",
                  analysis_method = "LM")

# FILTER BY LOWEST KRUSKAL WALLIS P-VALUE (data is not normal according to shapiro test)
k<-cat_plasmid_go %>% group_by(GO) %>% filter(n_distinct(status) >= 2) %>% summarise(k=kruskal.test(conteo_total ~ status)$p.value,                                     wilcoxon_p=list(pairwise.wilcox.test(conteo_total,status, exact=FALSE, p.adjust.method="BH")$p.value))  %>% arrange(k) %>% 
  slice_head(n = 9) %>%
  pull(GO)
##INDIVIDUAL TEST SEPARATELY
test<-cat_plasmid_go %>% filter(GO=="GO:0017144")
p<-pairwise.wilcox.test(test$conteo_total,test$status, exact=FALSE, p.adjust.method = "BH")

# Filter the dataframe for only the 9 categories with the lowest p value in Kruskal
filtered_cat_plasmid <- cat_plasmid_go %>%
  filter(GO %in% k)

# Create the boxplot graphs
plots <- filtered_cat_plasmid %>%
  ggplot(aes(x = status, y = conteo_total, fill = status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ GO, scales = "free_y") +
  labs(title = "Boxplots by Status for the 9 Categories with Highest Median",
       x = "Status",
       y = "Total Count") +
  theme_minimal()+
    scale_x_discrete(labels = function(x) ifelse(x == "barcode30", "Baseline", 
                                               ifelse(x == "barcode32", "24 hours",
                                                      ifelse(x=="barcode33", "10 days",x))))
```

### GO PLASMIDS CPM 
```{r}
cat_plasmid_cpm <- df_main  %>% filter(str_detect(mge, "PLASMID")) %>% inner_join(libra, by="sample") %>%mutate(cpm=(conteo/size)*1e6) %>% select(status,sample,GO, cpm)%>%separate_rows(GO, sep = ",") %>% 
  group_by(status, sample, GO) %>% summarise(conteo_total=sum(cpm)) %>%ungroup()

  analysis_maaslin2_boxplot(cat_plasmid_cpm,"GO",origin =  "plasmid")
```

### GO VIRUS CPM 
```{r}
cat_virus_cpm <- df_main %>% filter(str_detect(mge, "VIRUS")) %>% inner_join(libra, by="sample") %>% separate(mge, into = c("type_mge","length","rango"), sep=";",fill="right") %>% #now keep genes within provirus
  mutate(
    in_range = case_when(
      str_detect(type_mge, "PROVIRUS") ~ if_else(!is.na(rango), {
                  rango_parts <- str_split(rango, "-", simplify = TRUE)
                  rango_start <- as.numeric(rango_parts[, 1])
                  rango_end <- as.numeric(rango_parts[, 2])
                  start >= rango_start & end <= rango_end
                }, FALSE,FALSE #else is false
                ),TRUE ~ TRUE)
  ) %>%
  filter(in_range)%>%mutate(cpm=(conteo/size)*1e6) %>% select(status,sample,GO, cpm)%>%separate_rows(GO, sep = ",") %>% 
  group_by(status, sample, GO) %>% summarise(conteo_total=sum(cpm)) %>%ungroup()

analysis_maaslin2_boxplot(cat_virus_cpm,"GO","virus")
```

### GO Chromosomal contigs (NOT MGE)
```{r}
cat_main_cpm <- df_main  %>% filter(str_detect(mge, "NOT_MGE")) %>% inner_join(libra, by="sample") %>%mutate(cpm=(conteo/size)*1e6)%>% select(status,sample,GO, cpm)%>%separate_rows(GO, sep = ",") %>% 
  group_by(status, sample, GO) %>% summarise(conteo_total=sum(cpm)) %>%ungroup()
analysis_maaslin2_boxplot(cat_main_cpm,"GO","chromosome",  slice_range = c(10, 18))

```

# KEGG MAASLIN2 AND BOXPLOTS 

### KEGG RAW COUNT
```{r}
##PLASMID RAW COUNT
cat_plasmid_pathways<-df_main %>% 
    filter(str_detect(mge, "PLASMID"),pathways!="Not_assigned") %>% select(status,sample,pathways, conteo)%>%separate_rows(pathways, sep = ",")%>% filter(str_detect(pathways,"map")) %>% 
  group_by(status, sample, pathways) %>% summarise(conteo_total=sum(conteo)) %>% ungroup()  

plasmid_df <- cat_plasmid_pathways %>% 
  pivot_wider(names_from = pathways, values_from = conteo_total, values_fill = 0) %>% select(-1) %>% as.data.frame()
rownames(plasmid_df)<-plasmid_df[,1]
plasmid_df<-plasmid_df[,-1]

# Create the metadata dataframe
metadata <- cat_plasmid_pathways %>%
  select(sample, status) %>% distinct() %>% 
  mutate(name=ifelse(status == "barcode30", "Baseline", 
      ifelse(status == "barcode32", "24 hours",
      ifelse(status=="barcode33", "10 days",status)))) %>% as.data.frame()
rownames(metadata)<-metadata[,1]
metadata<-metadata[,-1]

fit_data=Maaslin2(input_data = plasmid_df,
                  input_metadata = metadata,
                      output = "Maaslin2_plasmid_pathways",       # Output folder for results
                  fixed_effects = c("name"), 
                  reference = "name,Baseline",
                  min_prevalence = 0.2,
                  normalization="TSS",
                  transform = "LOG",
                  analysis_method = "LM")

# FILTER BY LOWEST KRUSKAL WALLIS P-VALUE (data is not normal according to shapiro test)
k<-cat_plasmid_pathways %>% group_by(pathways) %>% filter(n_distinct(status) >= 2) %>% summarise(k=kruskal.test(conteo_total ~ status)$p.value,                                     wilcoxon_p=list(pairwise.wilcox.test(conteo_total,status, exact=FALSE, p.adjust.method="BH")$p.value))  %>% arrange(k) %>% 
  slice_head(n = 9) %>%
  pull(pathways)

# Filter the dataframe for only the 9 categories with the lowest p value in Kruskal
filtered_cat_plasmid <- cat_plasmid_pathways %>%
  filter(pathways %in% k)

# Create the boxplot graphs
plots <- filtered_cat_plasmid %>%
  ggplot(aes(x = status, y = conteo_total, fill = status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ pathways, scales = "free_y") +
  labs(title = "Boxplots by Status for the 9 Categories with Highest Median",
       x = "Status",
       y = "Total Count") +
  theme_minimal()+
    scale_x_discrete(labels = function(x) ifelse(x == "barcode30", "Baseline", 
                                               ifelse(x == "barcode32", "24 hours",
                                                      ifelse(x=="barcode33", "10 days",x))))
```

### KEGG PLASMIDS CPM 
```{r}
cat_plasmid_cpm <- df_main  %>% filter(str_detect(mge, "PLASMID"),pathways!="Not_assigned") %>% inner_join(libra, by="sample") %>%mutate(cpm=(conteo/size)*1e6) %>% select(status,sample,pathways, cpm)%>%separate_rows(pathways, sep = ",")%>% filter(str_detect(pathways,"map")) %>% 
  group_by(status, sample, pathways) %>% summarise(conteo_total=sum(cpm)) %>%ungroup()

analysis_maaslin2_boxplot(cat_plasmid_cpm,"pathways","plasmid")
```

### KEGG VIRUS CPM 
```{r}
cat_virus_cpm <- df_main %>% filter(str_detect(mge, "VIRUS"),pathways!="Not_assigned") %>% inner_join(libra, by="sample") %>% separate(mge, into = c("type_mge","length","rango"), sep=";",fill="right") %>% #now keep genes within provirus
  mutate(
    in_range = case_when(
      str_detect(type_mge, "PROVIRUS") ~ if_else(!is.na(rango), {
                  rango_parts <- str_split(rango, "-", simplify = TRUE)
                  rango_start <- as.numeric(rango_parts[, 1])
                  rango_end <- as.numeric(rango_parts[, 2])
                  start >= rango_start & end <= rango_end
                }, FALSE,FALSE #else is false
                ),TRUE ~ TRUE)
  ) %>%
  filter(in_range)%>%mutate(cpm=(conteo/size)*1e6) %>% select(status,sample,pathways, cpm)%>%separate_rows(pathways, sep = ",")%>% filter(str_detect(pathways,"map")) %>% 
  group_by(status, sample, pathways) %>% summarise(conteo_total=sum(cpm)) %>%ungroup()

analysis_maaslin2_boxplot(cat_virus_cpm,"pathways","virus")
```

### KEGG Chromosomal contigs (NOT MGE)
```{r}
cat_main_cpm <- df_main  %>% filter(str_detect(mge, "NOT_MGE"),pathways!="Not_assigned") %>% inner_join(libra, by="sample") %>%mutate(cpm=(conteo/size)*1e6)%>% select(status,sample,pathways, cpm)%>%separate_rows(pathways, sep = ",") %>% filter(str_detect(pathways,"map"))%>% 
  group_by(status, sample, pathways) %>% summarise(conteo_total=sum(cpm)) %>%ungroup()
gt<-analysis_maaslin2_boxplot(cat_main_cpm,"pathways","chromosome")
```

# BRITE ontology databse from KEGG MAASLIN2 AND BOXPLOTS 
```{r}

##PLASMID
cat_plasmid_cpm <- df_main  %>% filter(str_detect(mge, "PLASMID"),brite!="Not_assigned") %>% inner_join(libra, by="sample") %>%mutate(cpm=(conteo/size)*1e6) %>% select(status,sample,brite, cpm)%>%separate_rows(brite, sep = ",")%>% 
  group_by(status, sample, brite) %>% summarise(conteo_total=sum(cpm)) %>%ungroup()

g<-analysis_maaslin2_boxplot(cat_plasmid_cpm,"brite","plasmid")

###VIRUS
cat_virus_cpm <- df_main %>% filter(str_detect(mge, "VIRUS"),brite!="Not_assigned") %>% inner_join(libra, by="sample") %>% separate(mge, into = c("type_mge","length","rango"), sep=";",fill="right") %>% #now keep genes within provirus
  mutate(
    in_range = case_when(
      str_detect(type_mge, "PROVIRUS") ~ if_else(!is.na(rango), {
                  rango_parts <- str_split(rango, "-", simplify = TRUE)
                  rango_start <- as.numeric(rango_parts[, 1])
                  rango_end <- as.numeric(rango_parts[, 2])
                  start >= rango_start & end <= rango_end
                }, FALSE,FALSE #else is false
                ),TRUE ~ TRUE)
  ) %>%
  filter(in_range)%>%mutate(cpm=(conteo/size)*1e6) %>% select(status,sample,brite, cpm)%>%separate_rows(brite, sep = ",")%>%
  group_by(status, sample, brite) %>% summarise(conteo_total=sum(cpm)) %>%ungroup()

analysis_maaslin2_boxplot(cat_virus_cpm,"brite","virus")
###CHROMOSOMAL CONTIGS
cat_main_cpm <- df_main  %>% filter(str_detect(mge, "NOT_MGE"),brite!="Not_assigned") %>% inner_join(libra, by="sample") %>%mutate(cpm=(conteo/size)*1e6)%>% select(status,sample,brite, cpm)%>%separate_rows(brite, sep = ",") %>% 
  group_by(status, sample, brite) %>% summarise(conteo_total=sum(cpm)) %>%ungroup()
analysis_maaslin2_boxplot(cat_main_cpm,"brite","chromosome")
g

```
# PFAM (identical to kegg_ko)
```{r}

##PLASMID
cat_plasmid_cpm <- df_main  %>% filter(str_detect(mge, "PLASMID"),pfams!="Not_assigned") %>% inner_join(libra, by="sample") %>%mutate(cpm=(conteo/size)*1e6) %>% select(status,sample,pfams, cpm)%>%separate_rows(pfams, sep = ",")%>%  
  group_by(status, sample, pfams) %>% summarise(conteo_total=sum(cpm)) %>%ungroup()

analysis_maaslin2_boxplot(cat_plasmid_cpm,"pfams","plasmid")

###VIRUS
cat_virus_cpm <- df_main %>% filter(str_detect(mge, "VIRUS"),pfams!="Not_assigned") %>% inner_join(libra, by="sample") %>% separate(mge, into = c("type_mge","length","rango"), sep=";",fill="right") %>% #now keep genes within provirus
  mutate(
    in_range = case_when(
      str_detect(type_mge, "PROVIRUS") ~ if_else(!is.na(rango), {
                  rango_parts <- str_split(rango, "-", simplify = TRUE)
                  rango_start <- as.numeric(rango_parts[, 1])
                  rango_end <- as.numeric(rango_parts[, 2])
                  start >= rango_start & end <= rango_end
                }, FALSE,FALSE #else is false
                ),TRUE ~ TRUE)
  ) %>%
  filter(in_range)%>%mutate(cpm=(conteo/size)*1e6) %>% select(status,sample,pfams, cpm)%>%separate_rows(pfams, sep = ",")%>%
  group_by(status, sample, pfams) %>% summarise(conteo_total=sum(cpm)) %>%ungroup()

analysis_maaslin2_boxplot(cat_virus_cpm,"pfams","virus")
###CHROMOSOMAL CONTIGS
cat_main_cpm <- df_main  %>% filter(str_detect(mge, "NOT_MGE"),pfams!="Not_assigned") %>% inner_join(libra, by="sample") %>%mutate(cpm=(conteo/size)*1e6)%>% select(status,sample,pfams, cpm)%>%separate_rows(pfams, sep = ",") %>% 
  group_by(status, sample, pfams) %>% summarise(conteo_total=sum(cpm)) %>%ungroup()
analysis_maaslin2_boxplot(cat_main_cpm,"pfams","chromosome")



```

