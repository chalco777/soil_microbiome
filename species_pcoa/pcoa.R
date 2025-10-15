# PCoA Analysis Script
# 
# This script performs Principal Coordinates Analysis (PCoA) on taxonomic abundance data
# to visualize beta diversity patterns across samples. As it is known, PCoA is an ordination method that
# reduces multidimensional data into a 2D space while preserving distances between samples.
# It helps identify sample clustering patterns and assess compositional dissimilarity based
# on Bray-Curtis distances. The 2nd PcoA filters taxa based on prevalence and abundance
# thresholds to focus on relevant biological signals.

setwd("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/alen-belem")
library(tidyverse)
library(vegan)
library(glue)
library(ggrepel)

df <- read_tsv("relative_abundances_per_sample.txt")

# TRANSPOSE: rows are now samples and columns are taxa
df <- t(df)
# Move first row to column names
colnames(df) <- df[1,]
df <- df[-1,]
# Convert to data frame
df <- as.data.frame(df)

# Remove extra barcodes, keep ONLY those used for coassembly groups
df <- df %>%
  filter(!(rownames(df) %in% c("barcode11", "barcode13","barcode18", paste0("barcode",c(20:22,30:33)))))

# Create a tibble to define coassembly groups for each barcode
df_variable <- tibble(
  samples = rownames(df[order(rownames(df)), ]),  # Assign sample names to 'samples' column
  group = c(rep("barcode30", 6), "barcode31", rep("barcode32", 3), "barcode32",
            "barcode31", rep("barcode33", 3), "barcode33")
)

# WITH ALL SAMPLES
# df_variable <- tibble(
#   samples = rownames(df[order(rownames(df)), ]),
#   group = c(rep("barcode30", 6), "barcode31", rep("barcode32", 3), "None", "barcode32", "None",
#             "barcode31", rep("barcode33", 3), "None", "barcode33", rep("None",3),paste0("barcode_off_",c(30:33))
# ))


### FILTERING: taxa must have a total sum > 0.001 AND be present in at least 3 barcodes
### Otherwise, discard the taxon and exclude it from PCoA analysis

# Convert to numeric to filter columns summing to 0 (taxa absent in all barcodes)
dim(df)
df[] <- lapply(df, as.numeric)
dim(df)
# Filter taxa with total abundance > 0.001
df <- df[, which(colSums(df) > 0.001)]
dim(df)
# Filter taxa present in at least 3 samples
df <- df[, apply(df, 2, function(col) sum(col > 0) >= 3)]

df[is.na(df)] <- 0

# Convert to matrix for vegan
df <- as.matrix(df)
# Calculate Bray-Curtis distance
dist <- vegdist(df, method = "bray")

# Perform PCoA with cmdscale on distance matrix, add correction and get eigenvalues
pcoa <- cmdscale(dist, k = 2, eig = TRUE, add = TRUE)
positions <- pcoa$points  # Extract vectors
colnames(positions) <- c("pcoa1", "pcoa2")  # Name coordinates
positions <- positions[order(rownames(positions)), ]  # Order by row name

# Calculate variance explained
expl <- (pcoa$eig / sum(pcoa$eig)) * 100
exp <- format(round(expl[1:2], digits = 1), nsmall = 1, trim = TRUE)
labs <- c(glue("PCo 1 ({exp[1]}%)"),
          glue("PCo 2 ({exp[2]}%)"))

# Convert ordered PCoA results to tibble and join with df_variable
# Modify group labels
df_variable <- df_variable %>%
  mutate(group = case_when(
    group == "barcode30" ~ "baseline",
    group == "barcode32" ~ "24 hours",
    group == "barcode33" ~ "10 days",
    TRUE ~ NA_character_  # Filter other values
  )) %>%
  filter(!is.na(group))  # Keep only selected groups

# Create PCoA plot
pcoa_plot <- positions %>%
  as_tibble(rownames = "samples") %>%
  filter(!(rownames(positions) %in% c("barcode14", "barcode07"))) %>% 
  left_join(df_variable, by = "samples") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = group)) +
  geom_point(size = 2.5, show.legend = TRUE) +
  geom_text_repel(aes(label = samples), size = 3, max.overlaps = 2, force = 2) +
  labs(color = "Group", x = labs[1], y = labs[2]) +  
  theme_minimal(base_size = 16) +  
  theme(
    axis.title = element_text(size = 16, face = "bold"),  
    axis.text = element_text(size = 14),  
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12),  
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  
  ) +
  scale_color_brewer(palette = "Set1") +  
  ggtitle("PCoA Analysis with selected barcodes")
  #+ xlim(c(-0.105, -0.095))


# Display the plot
print(pcoa_plot)

# Save plot in high resolution
ggsave(
  filename = "pcoa_analysis_selected_barcodes.png",
  plot = pcoa_plot,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

#### FILTERED ANALYSIS ####

df <- read_tsv("relative_abundances_per_sample.txt")

# TRANSPOSE: rows are now samples and columns are taxa
df <- t(df)
# Move first row to column names
colnames(df) <- df[1,]
df <- df[-1,]
# Convert to data frame
df <- as.data.frame(df)

# Remove extra barcodes, keep ONLY those used for coassembly groups
df <- df %>%
  filter(!(rownames(df) %in% c("barcode11", "barcode13","barcode18", paste0("barcode",c(20:22,30:33)))))

# Create a tibble to define coassembly groups for each barcode
df_variable <- tibble(
  samples = rownames(df[order(rownames(df)), ]),  # Assign sample names to 'samples' column
  group = c(rep("barcode30", 6), "barcode31", rep("barcode32", 2),"barcode32","barcode32_bio",
            "barcode31", rep("barcode33", 3), "barcode33_bio")
)

### FILTERING: taxa must be present with total sum > 0.001 AND in at least 3 barcodes
### Otherwise, discard the taxon and exclude it from PCoA analysis

# Convert to numeric to filter columns summing to 0 (taxa absent in all barcodes)
dim(df)
df[] <- lapply(df, as.numeric)

df[is.na(df)] <- 0

# Convert to matrix for vegan
df <- as.matrix(df)
# Calculate Bray-Curtis distance
dist <- vegdist(df, method = "bray")

# Perform PCoA with cmdscale on distance matrix, add correction and get eigenvalues
pcoa <- cmdscale(dist, k = 2, eig = TRUE, add = TRUE)
positions <- pcoa$points  # Extract vectors
colnames(positions) <- c("pcoa1", "pcoa2")  # Name coordinates
positions <- positions[order(rownames(positions)), ]  # Order by row name

# Calculate variance explained
expl <- (pcoa$eig / sum(pcoa$eig)) * 100
exp <- format(round(expl[1:2], digits = 1), nsmall = 1, trim = TRUE)
labs <- c(glue("PCo 1 ({exp[1]}%)"),
          glue("PCo 2 ({exp[2]}%)"))

# Convert ordered PCoA results to tibble and join with df_variable
# Modify group labels
df_variable <- df_variable %>%
  mutate(group = case_when(
    group == "barcode30" ~ "baseline",
    group == "barcode32" ~ "24 hours",
    group == "barcode33" ~ "10 days",
    group == "barcode32_bio" ~ "24 hours bio",
    group == "barcode33_bio" ~ "10 days bio",
    TRUE ~ NA_character_  # Filter other values
  )) %>%
  filter(!is.na(group))  # Keep only selected groups

# Create PCoA plot
pcoa_2 <- positions %>%
  as_tibble(rownames = "samples") %>%
  filter(!(rownames(positions) %in% c("barcode14", "barcode07"))) %>% 
  left_join(df_variable, by = "samples") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = group)) +
  geom_point(size = 2.5, show.legend = TRUE) +
  geom_text_repel(aes(label = samples), size = 3, max.overlaps = 2, force = 2) +
  labs(color = "Group", x = labs[1], y = labs[2]) +  
  theme_minimal(base_size = 16) +  
  theme(
    axis.title = element_text(size = 16, face = "bold"),  
    axis.text = element_text(size = 14),  
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12),  
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  
  ) +
  scale_color_brewer(palette = "Set1") +  
  ggtitle("PCoA Analysis with selected barcodes")
  #+ xlim(c(-0.105, -0.095))

# Save plot in high resolution
ggsave(
  filename = "pcoa_analysis_selected_barcodes_filtered_sp.png",
  plot = pcoa_2,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)