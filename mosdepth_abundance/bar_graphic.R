#Bar graphic
library(tidytext)
library(tidyverse)
plasmid<-read_tsv("plasmid_depth_genes_all_barcodes_official.tsv",
          col_names = c("Barcode","Contig","Feature","Start","End","Product","Depth"))
virus<-read_tsv("viurs_depth_genes_all_barcodes_official.tsv",
        col_names = c("Barcode","Contig","Feature","Start","End","Product","Depth"))


str(plasmid)
library(dplyr)
library(ggplot2)

# Filter barcodes 30 to 33 and calculate mean depth per product for each barcode
top_products <- plasmid %>%
  filter(Barcode %in% c("barcode30", "barcode31", "barcode32", "barcode33")) %>%
  filter(Product != "hypothetical protein") %>% 
  group_by(Barcode, Product) %>%
  summarize(Relative_Abundance = sum(Depth, na.rm = TRUE)) %>%
  arrange(desc(Relative_Abundance)) %>%
  group_by(Barcode) %>% 
  top_n(10, Relative_Abundance) %>%  # Select the 10 most abundant products for each barcode
  mutate(Product = reorder_within(Product, Relative_Abundance, Barcode))

  # Create horizontal bar plot for each barcode with colors and centered title
ggplot(top_products, aes(x = Product, y = Relative_Abundance, fill = Barcode)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Make the plot horizontal
  facet_wrap(~ Barcode, scales = "free") + # Create a panel for each barcode
  labs(title = "Top 10 products with highest depth (relative abundance) in plasmids",
     x = "Product",
     y = "Mean depth (Depth)") +
  scale_fill_brewer(palette = "Set3") + # Apply a color palette
  theme_minimal() +
  theme(
  plot.title = element_text(hjust = 0.5, size = 12, face = "bold", color = "#333333"), # Center and adjust the title
  axis.title.x = element_text(hjust=0,size = 10, face = "bold", color = "#666666"), # Style for X axis title
  axis.title.y = element_text(size = 10, face = "bold", color = "#666666"), # Style for Y axis title
  axis.text.x = element_text(size = 8, color = "#444444"), # Style for X axis text
  axis.text.y = element_text(size = 6, color = "#444444"), # Style for Y axis text
  strip.text = element_text(size = 8, face = "bold", color = "#555555"), # Style for facet text
  legend.position = "none", # Hide legend if not needed
  plot.background = element_rect(fill = "white", color = NA) # White background
  )+
  scale_x_discrete(labels = function(x) gsub("__.*$", "", x)) # Remove "__barcodeXX" from labels




# Filter barcodes 30 to 33 and calculate mean depth per product for each barcode
top_products <- virus %>%
  filter(Barcode %in% c("barcode30", "barcode31", "barcode32", "barcode33")) %>%
  filter(Product != "hypothetical protein") %>% 
  group_by(Barcode, Product) %>%
  summarize(Relative_Abundance = sum(Depth, na.rm = TRUE)) %>%
  arrange(desc(Relative_Abundance)) %>%
  group_by(Barcode) %>%
  top_n(10, Relative_Abundance) # Select the 10 most abundant products for each barcode

# Create horizontal bar plot for each barcode with colors and centered title
ggplot(top_products, aes(x = reorder(Product, Relative_Abundance), y = Relative_Abundance, fill = Barcode)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Make the plot horizontal
  facet_wrap(~ Barcode, scales = "free") + # Create a panel for each barcode
  labs(title = "Top 10 products with highest depth (relative abundance) in viruses",
     x = "Product",
     y = "Mean depth (Depth)") +
  scale_fill_brewer(palette = "Set3") + # Apply a color palette
  theme_minimal() +
  theme(
  plot.title = element_text(hjust = 0.5, size = 12, face = "bold", color = "#333333"), # Center and adjust the title
  axis.title.x = element_text(hjust=0,size = 10, face = "bold", color = "#666666"), # Style for X axis title
  axis.title.y = element_text(size = 10, face = "bold", color = "#666666"), # Style for Y axis title
  axis.text.x = element_text(size = 6, color = "#444444"), # Style for X axis text
  axis.text.y = element_text(size = 6, color = "#444444"), # Style for Y axis text
  strip.text = element_text(size = 8, face = "bold", color = "#555555"), # Style for facet text
  legend.position = "none", # Hide legend if not needed
  plot.background = element_rect(fill = "white", color = NA) # White background
  )
