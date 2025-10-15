# This script analyzes plasmid gene abundance data from metagenomic assemblies
# It processes COG (Clusters of Orthologous Groups) annotations, creates LEfSe input files,
# generates boxplots comparing gene categories across different timepoints (baseline, 24h, 10 days),
# and analyzes file sizes across samples

# IMPORTANCE OF FILE SIZE ANALYSIS:
# Before comparing gene abundances across samples, it is important to evaluate the read file 
# sizes (sequencing depth) of each barcode within each group. Differences in sequencing depth 
# can significantly bias abundance comparisons, as samples with higher sequencing coverage 
# will naturally have higher raw gene counts. Thus, we can:
# 1) Identify samples with unusually low/high coverage that may need to be excluded 
# 2) Determine if normalization is necessary
# 3) Ensure that observed differences in gene abundance reflect true biological variation 
#    rather than technical artifacts coming from unequal sequencing depth
# 4) Select appropriate normalization strategies (e.g., relative abundance, TPM, or rarefaction)

setwd("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/alen-belem/eggnog_abundance")
library(tidyverse)
library(gt)
library(stringr)
library(tidytext)
plasmid<-read_tsv("merged_cog_results_plasmid.gff",
                  col_names = c("Barcode","Contig","Feature","Start","End","Name", "Eggnog","COG","Depth"))
## Assign definitive categories
plasmid_new <- plasmid %>%
  mutate(new_name = ifelse(!is.na(Eggnog) & Eggnog != "None", Eggnog, 
                           ifelse(!is.na(Name), Name, Feature))) %>%
  select(Barcode, new_name, Depth)
## Sum the depth of genes with the same name within the same barcode
plasmid_new<-plasmid_new %>% group_by(Barcode,new_name) %>% 
  summarize(sum_depth=sum(Depth))
## Convert to wide format
lefse<-plasmid_new %>% 
  pivot_wider(names_from = Barcode, values_from = sum_depth) %>%
  replace(is.na(.), 0)

status = c("status",rep("barcode30", 6), "barcode31", rep("barcode32", 3), NA, "barcode32", NA,
           "barcode31", rep("barcode33", 3), NA, "barcode33", rep(NA,5))
# Add status and select coassemblies of interest
lefse_new <- rbind(status, lefse)  %>%
  mutate(new_name = gsub("\\s", "_", new_name)) %>%  # Replace spaces with underscores in new_name column
  select(new_name, where(~ .[1] %in% c("barcode30", "barcode32")))

write.table(lefse_new, file = "for_lefse_plasmid.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)



## FOR COG

plasmid_cog<-read_tsv("merged_cog_results_plasmid.gff",
                      col_names = c("Barcode","Contig","Feature","Start","End","Name", "Eggnog","cog","Depth"))
## Assign definitive categories
plasmid_new_cog <- plasmid_cog %>%
  mutate(new_name = ifelse(!is.na(Eggnog) & Eggnog != "None", cog, 
                           ifelse(!is.na(Name), Name, Feature))) %>%
  select(Barcode, new_name, Depth)
## Sum the depth of genes with the same name within the same barcode
plasmid_new_cog<-plasmid_new_cog %>% group_by(Barcode,new_name) %>% 
  summarize(sum_depth=sum(Depth))
## Convert to wide format
lefse_cog<-plasmid_new_cog %>% 
  pivot_wider(names_from = Barcode, values_from = sum_depth) %>%
  replace(is.na(.), 0)

status = c("status",rep("barcode30", 6), "barcode31", rep("barcode32", 3), NA, "barcode32", NA,
           "barcode31", rep("barcode33", 3), NA, "barcode33", rep(NA,5))
# Add status and select coassemblies of interest
lefse_new_cog <- rbind(status, lefse_cog)  %>%
  mutate(new_name = gsub("\\s", "_", new_name)) %>%  # Replace spaces with underscores in new_name column
  select(new_name, where(~ .[1] %in% c("barcode30", "barcode32","barcode33")))
## Rows that are present in all samples
lefse_filtered_cog <- lefse_new_cog %>%
  filter(rowSums(select(., -new_name) != 0) > 0)

write.table(lefse_filtered_cog, file = "for_lefse_cog_plasmid.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)


### BOXPLOT

lefse_filtered_cog2<-lefse_filtered_cog[-1,]
data<-lefse_filtered_cog2 %>% mutate(across(-new_name,as.numeric))

boxlef<-lefse_filtered_cog2 %>% 
  pivot_longer(cols=!new_name,
               names_to = "sample",values_to = "depth") %>% 
  mutate(status = case_when(
    sample %in% c('barcode01', 'barcode02', 'barcode03', 'barcode04', 'barcode05','barcode06') ~ 'barcode30',
    sample %in% c('barcode08', 'barcode09', 'barcode10', 'barcode12') ~ 'barcode32'  ,
    TRUE ~ 'barcode33'
  ))  %>% 
  mutate(depth = as.numeric(depth))

p1 <- ggplot(boxlef, aes(x = status, y = depth, fill = status)) +
  # Add boxplot with black borders and color by status
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.6) +
  # Add jitter points without custom colors
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  # Adjust labels and title
  labs(title = "Depth by Coassembly for 'Energy production and conversion'",
       x = "Coassembly",
       y = "Depth") +
  # Apply minimalist theme
  theme_minimal() +
  # Remove legend
  guides(fill = "none") +  # Alternatively: theme(legend.position = "none")
  # Modify Y-axis labels to assign "baseline" and "24 hours"
  scale_x_discrete(labels = function(x) ifelse(x == "barcode30", "baseline", 
                                               ifelse(x == "barcode32", "24 hours",
                                                      ifelse(x=="barcode33","10 days",x)))) +
  # Customize title and labels
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank()
  )
ggsave("boxplot_energy_production_conversion.png", plot = p1, width = 8, height = 6, dpi = 300)

# Step 5: Create second boxplot (without 'barcode03')
data_C_no_barcode03 <- boxlef %>% filter(sample != 'barcode03')
ggplot(data_C_no_barcode03, aes(x = status, y = depth, fill = status)) +
  # Add boxplot with black borders and color by status
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.6) +
  # Add jitter points without custom colors
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  # Adjust labels and title
  labs(title = "Depth by Coassembly without 'barcode03' for 'Energy production and conversion'",
       x = "Coassembly",
       y = "Depth") +
  # Apply minimalist theme
  theme_minimal() +
  # Remove legend
  guides(fill = "none") +  # Alternatively: theme(legend.position = "none")
  # Modify Y-axis labels to assign "baseline" and "24 hours"
  scale_x_discrete(labels = function(x) ifelse(x == "barcode30", "Baseline", 
                                               ifelse(x == "barcode32", "24 hours",
                                                      ifelse(x=="barcode33", "10 days",x)))) +
  # Customize title and labels
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank()
  )


#### BAR PLOT FOR ALL CATEGORIES


boxlef<-boxlef %>% rename(category=new_name) %>% filter(category!="hypothetical_protein")

# Load necessary libraries
library(dplyr)
library(ggplot2)

# 1. Calculate median depth by category and select the 8 categories with highest median
top_categories <- boxlef %>%
  group_by(category) %>%
  summarise(median_depth = median(depth, na.rm = TRUE)) %>%
  arrange(desc(median_depth)) %>%
  slice(1:9) %>%
  pull(category)


# 2. Filter data to include only selected categories
filtered_data <- boxlef %>%
  filter(category %in% top_categories)

# Dictionary with real names of COG categories
real_names <- c("C" = "Energy production and conversion",
                "E" = "Amino acid transport and metabolism",
                "H" = "Coenzyme transport and metabolism",
                "I" = "Lipid transport and metabolism",
                "J" = "Translation, ribosomal structure and biogenesis",
                "L" = "Replication, recombination and repair",
                "P" = "Inorganic ion transport and metabolism",
                "S" = "Function unknown",
                "M"="Cell wall/membrane/envelope biogenesis")

# Rename categories in 'category' column
filtered_data <- filtered_data %>%
  mutate(category = recode(category, !!!real_names))


p3 <- ggplot(filtered_data, aes(x = status, y = depth, fill=status)) +
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.6) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  facet_wrap(~ category, scales = "free_y") +
  theme_minimal() +
  labs(title = "Depth of each Category before and after treatment",
       x = NULL,  # Don't show x-axis label
       y = "Mean depth",
       color = "Sample") +
  scale_x_discrete(labels = function(x) ifelse(x == "barcode30", "Baseline", 
                                               ifelse(x == "barcode32", "24 hours",
                                                      ifelse(x=="barcode33", "10 days",x)))) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 9),  # Reduce x-axis label text size
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 9),  # Reduce facet text size
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside"
  ) +
  theme(axis.text.x = element_blank()) + # Hide x-axis text
  guides(fill = "none")
ggsave("boxplot_all_categories_faceted.png", plot = p3, width = 12, height = 10, dpi = 300)



data_all_no_barcode03 <- boxlef %>% filter(sample != 'barcode03')
boxlef2<-data_all_no_barcode03

# 1. Calculate median depth by category and select the 8 categories with highest median
top_categories <- boxlef2 %>%
  group_by(category) %>%
  summarise(median_depth = median(depth, na.rm = TRUE)) %>%
  arrange(desc(median_depth)) %>%
  slice(1:9) %>%
  pull(category)


# 2. Filter data to include only selected categories
filtered_data2 <- boxlef2 %>%
  filter(category %in% top_categories)

# Dictionary with real names of COG categories
real_names <- c("C" = "Energy production and conversion",
                "E" = "Amino acid transport and metabolism",
                "H" = "Coenzyme transport and metabolism",
                "I" = "Lipid transport and metabolism",
                "J" = "Translation, ribosomal structure and biogenesis",
                "L" = "Replication, recombination and repair",
                "P" = "Inorganic ion transport and metabolism",
                "S" = "Function unknown",
                "M"="Cell wall/membrane/envelope biogenesis")

# Rename categories in 'category' column
filtered_data2 <- filtered_data2 %>%
  mutate(category = recode(category, !!!real_names))


p4 <- ggplot(filtered_data2, aes(x = status, y = depth, fill=status)) +
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.6) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  facet_wrap(~ category, scales = "free_y") +
  theme_minimal() +
  labs(title = "Depth of each Category (without Barcode03) before and after treatment",
       x = NULL,  # Don't show x-axis label
       y = "Mean depth",
       color = "Sample") +
  scale_x_discrete(labels = function(x) ifelse(x == "barcode30", "Baseline", 
                                               ifelse(x == "barcode32", "24 hours",
                                                      ifelse(x=="barcode33", "10 days",x)))) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 9),  # Reduce x-axis label text size
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 9),  # Reduce facet text size
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside"
  ) +
  theme(axis.text.x = element_blank()) + # Hide x-axis text
  guides(fill = "none")
ggsave("boxplot_all_categories_faceted_no_barcode03.png", plot = p4, width = 12, height = 10, dpi = 300)



##### SIZE PLOT


file<-read_tsv("tamaño_barcode.tsv", col_names = c("Tamaño","Muestra"))
modified_data <- file %>%
  mutate(Muestra = str_extract(Muestra, "barcode\\d+")) %>%  
  filter(!(Muestra %in% c(paste0("barcode",30:33),"barcode07","barcode11","barcode13","barcode14","barcode18",paste0("barcode",20:22)))) %>% 
  mutate(status = case_when(
    Muestra %in% c('barcode01', 'barcode02', 'barcode03', 'barcode04', 'barcode05','barcode06') ~ 'barcode30',
    Muestra %in% c('barcode08', 'barcode09', 'barcode10', 'barcode12') ~ 'barcode32'  ,
    TRUE ~ 'barcode33'
  )) 



# Convert 'Tamaño' to numeric in Gigabytes
modified_data$Tamaño_GB <- as.numeric(sub("G", "", modified_data$Tamaño))
# Ensure 'status' is a factor with specified order
modified_data$status <- factor(modified_data$status, levels = c("barcode30", "barcode32", "barcode33"))

# Load necessary libraries
library(ggplot2)
library(ggpubr)

# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(Tamaño_GB ~ status, data = modified_data)
print(kruskal_result)

# Perform pairwise comparisons with Wilcoxon test
pairwise_wilcox <- pairwise.wilcox.test(modified_data$Tamaño_GB, modified_data$status, p.adjust.method = "BH")
print(pairwise_wilcox)

# Define custom colors for each 'status'
status_colors <- c("barcode30" = "#E69F00", "barcode32" = "#56B4E9", "barcode33" = "#009E73")

# Create plot
p5 <- ggplot(modified_data, aes(x = status, y = Tamaño_GB, fill = status)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.2, size = 3, alpha = 0.7, color = "black") +
  scale_fill_manual(values = status_colors) +
  labs(
    title = "Size Distribution by Status",
    x = NULL,
    y = "Size (GB)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  ) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(modified_data$Tamaño_GB) * 1.15,
    size = 5
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      c("barcode30", "barcode32"),
      c("barcode30", "barcode33"),
      c("barcode32", "barcode33")
    ),
    label = "p.signif",
    label.y = c(
      max(modified_data$Tamaño_GB) * 1.05,
      max(modified_data$Tamaño_GB) * 1.10,
      max(modified_data$Tamaño_GB) * 1.15
    ),
    size = 5
  )+
  scale_x_discrete(labels = function(x) ifelse(x == "barcode30", "Baseline", 
                                               ifelse(x == "barcode32", "24 hours",
                                                      ifelse(x=="barcode33", "10 days",x))))
ggsave("boxplot_size_distribution_by_status.png", plot = p5, width = 10, height = 8, dpi = 300)