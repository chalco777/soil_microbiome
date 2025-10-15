

# First we set the working directory and load required libraries
setwd("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/alen-belem/eggnog_abundance")
library(tidyverse)
library(gt)
library(stringr)
library(tidytext)

# Load plasmid annotation data with column names
general<-read_tsv("merged_results_all_barcodes.gff",
          col_names = c("Barcode","Contig","Feature","Start","End","Name", "Eggnog","Depth"))

# Alternative: Uncomment below to load a the virus dataset
#virus<-read_tsv("virus_depth_genes_all_barcodes_official.tsv",
#               col_names = c("Barcode","Contig","Feature","Start","End","Product","Depth"))

plasmid<-read_tsv("plasmid_merged_results_all_barcodes.gff",
          col_names = c("Barcode","Contig","Feature","Start","End","Name", "Eggnog","Depth"))

# ----------STATISTICS-------------------

# Count proteins without EggNOG annotation per barcode
not_asigned_egg <- plasmid %>%
  mutate(Eggnog = ifelse(is.na(Eggnog) | Eggnog == "None", "None", Eggnog)) %>%
  group_by(Barcode) %>%
  summarize(not_egg = sum(Eggnog == "None")) #missing to add those that are None

# Count hypothetical proteins per barcode
hypo<-plasmid %>%  
  filter(Name=="hypothetical protein") %>%
  group_by(Barcode) %>%
  summarise(num_hypo=length(Name)) %>%tail()

# Count total predicted proteins per barcode
total_protein<-plasmid %>% group_by(Barcode) %>% 
  summarise(num_total_protein=length(Name))


# Combine statistics and calculate proportions for each barcode
statistics<-inner_join(total_protein,hypo,by="Barcode") %>% 
  inner_join(not_asigned_egg, by="Barcode") %>% 
  mutate(prop_hypo=num_hypo/num_total_protein) %>% 
  mutate(prop_not_egg=not_egg/num_total_protein) %>%  tail(4)


# --- Here we create and save a summary table ---

tab<-statistics %>%  rename(
  "Total Predicted Proteins" = num_total_protein,
  "Hypothetical Proteins" = num_hypo,
  "Not EggNog Annotations" = not_egg,
  "Proportion Hypothetical" = prop_hypo,
  "Proportion Not EggNog Annotated" = prop_not_egg
) %>% gt() %>%
  tab_header(
  title = "plasmid: Statistics of Functional Gene Prediction Analysis by Barcode",
  ) %>%
  fmt_number(
  columns = c("Total Predicted Proteins", "Hypothetical Proteins", "Not EggNog Annotations"),
  decimals = 0
  ) %>%
  fmt_percent(
  columns = c("Proportion Hypothetical", "Proportion Not EggNog Annotated"),
  decimals = 1
  )
gtsave(tab, filename = "plasmid_statistics_table.pdf")

# Filter barcodes 30 to 33 and calculate mean depth per product for each barcode
top_products <- plasmid %>%
  filter(Barcode %in% c("barcode30", "barcode31", "barcode32", "barcode33")) %>%
  filter(!is.na(Eggnog) & Eggnog != "None") %>% 
  group_by(Barcode, Eggnog) %>%
  summarize(Relative_Abundance = sum(Depth, na.rm = TRUE)) %>%
  arrange(desc(Relative_Abundance)) %>%
  group_by(Barcode) %>% 
  top_n(10, Relative_Abundance) %>%  # Select the 10 most abundant products for each barcode
  mutate(Eggnog = substr(Eggnog, 1, 45)) %>% mutate(Eggnog = reorder_within(Eggnog, Relative_Abundance, Barcode))

# Create horizontal bar plot for each barcode with colors and centered title
ggplot(top_products, aes(x = Eggnog, y = Relative_Abundance, fill = Barcode)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Make the plot horizontal
  facet_wrap(~ Barcode, scales = "free") + # Create a panel for each barcode
  labs(title = "Top 10 classes with highest depth (relative abundance) in plasmid",
     x = "Product",
     y = "Total depth (Depth)") +
  scale_fill_brewer(palette = "Set3") + # Apply a color palette
  theme_minimal() +
  theme(
  plot.title = element_text(hjust = 0, size = 12, face = "bold", color = "#333333"), # Center and adjust the title
  axis.title.x = element_text(hjust=0,size = 10, face = "bold", color = "#666666"), # Style for X axis title
  axis.title.y = element_text(size = 10, face = "bold", color = "#666666"), # Style for Y axis title
  axis.text.x = element_text(size = 8, color = "#444444"), # Style for X axis text
  axis.text.y = element_text(size = 6, color = "#444444"), # Style for Y axis text
  strip.text = element_text(size = 8, face = "bold", color = "#555555"), # Style for facet text
  legend.position = "none", # Hide legend if not needed
  plot.background = element_rect(fill = "white", color = NA) # White background
  )+
  scale_x_discrete(labels = function(x) gsub("__.*$", "", x)) # Remove "__barcodeXX" from labels
