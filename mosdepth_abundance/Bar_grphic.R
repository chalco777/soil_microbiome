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

# Filtrar los barcodes 30 a 33 y calcular la profundidad media por producto para cada barcode
top_products <- plasmid %>%
  filter(Barcode %in% c("barcode30", "barcode31", "barcode32", "barcode33")) %>%
  filter(Product != "hypothetical protein") %>% 
  group_by(Barcode, Product) %>%
  summarize(Relative_Abundance = sum(Depth, na.rm = TRUE)) %>%
  arrange(desc(Relative_Abundance)) %>%
  group_by(Barcode) %>% 
  top_n(10, Relative_Abundance) %>%  # Seleccionar los 10 productos más abundantes para cada barcode
  mutate(Product = reorder_within(Product, Relative_Abundance, Barcode))

  # Crear el gráfico de barras horizontal para cada barcode con colores y título centrado
ggplot(top_products, aes(x = Product, y = Relative_Abundance, fill = Barcode)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Hacer el gráfico horizontal
  facet_wrap(~ Barcode, scales = "free") + # Crear un panel para cada barcode
  labs(title = "Top 10 productos con mayor depth (relative abundance) en plásmidos",
       x = "Producto",
       y = "Profundidad media (Depth)") +
  scale_fill_brewer(palette = "Set3") + # Aplicar una paleta de colores
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold", color = "#333333"), # Centrar y ajustar el título
    axis.title.x = element_text(hjust=0,size = 10, face = "bold", color = "#666666"), # Estilo del título del eje X
    axis.title.y = element_text(size = 10, face = "bold", color = "#666666"), # Estilo del título del eje Y
    axis.text.x = element_text(size = 8, color = "#444444"), # Estilo de los textos del eje X
    axis.text.y = element_text(size = 6, color = "#444444"), # Estilo de los textos del eje Y
    strip.text = element_text(size = 8, face = "bold", color = "#555555"), # Estilo del texto de las facetas
    legend.position = "none", # Ocultar la leyenda si no es necesaria
    plot.background = element_rect(fill = "white", color = NA) # Fondo blanco
  )+
  scale_x_discrete(labels = function(x) gsub("__.*$", "", x)) # Remover "__barcodeXX" de las etiquetas




# Filtrar los barcodes 30 a 33 y calcular la profundidad media por producto para cada barcode
top_products <- virus %>%
  filter(Barcode %in% c("barcode30", "barcode31", "barcode32", "barcode33")) %>%
  filter(Product != "hypothetical protein") %>% 
  group_by(Barcode, Product) %>%
  summarize(Relative_Abundance = sum(Depth, na.rm = TRUE)) %>%
  arrange(desc(Relative_Abundance)) %>%
  group_by(Barcode) %>%
  top_n(10, Relative_Abundance) # Seleccionar los 10 productos más abundantes para cada barcode

# Crear el gráfico de barras horizontal para cada barcode con colores y título centrado
ggplot(top_products, aes(x = reorder(Product, Relative_Abundance), y = Relative_Abundance, fill = Barcode)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Hacer el gráfico horizontal
  facet_wrap(~ Barcode, scales = "free") + # Crear un panel para cada barcode
  labs(title = "Top 10 productos con mayor depth (relative abundance) en virus",
       x = "Producto",
       y = "Profundidad media (Depth)") +
  scale_fill_brewer(palette = "Set3") + # Aplicar una paleta de colores
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold", color = "#333333"), # Centrar y ajustar el título
    axis.title.x = element_text(hjust=0,size = 10, face = "bold", color = "#666666"), # Estilo del título del eje X
    axis.title.y = element_text(size = 10, face = "bold", color = "#666666"), # Estilo del título del eje Y
    axis.text.x = element_text(size = 6, color = "#444444"), # Estilo de los textos del eje X
    axis.text.y = element_text(size = 6, color = "#444444"), # Estilo de los textos del eje Y
    strip.text = element_text(size = 8, face = "bold", color = "#555555"), # Estilo del texto de las facetas
    legend.position = "none", # Ocultar la leyenda si no es necesaria
    plot.background = element_rect(fill = "white", color = NA) # Fondo blanco
  )





