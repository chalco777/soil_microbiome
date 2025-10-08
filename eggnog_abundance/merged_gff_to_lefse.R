setwd("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/alen-belem/eggnog_abundance")
library(tidyverse)
library(gt)
library(stringr)
library(tidytext)
plasmid<-read_tsv("merged_cog_results_plasmid.gff",
                  col_names = c("Barcode","Contig","Feature","Start","End","Name", "Eggnog","Depth"))
##Asignamos categorias definitivas
plasmid_new <- plasmid %>%
  mutate(new_name = ifelse(!is.na(Eggnog) & Eggnog != "None", Eggnog, 
                           ifelse(!is.na(Name), Name, Feature))) %>%
  select(Barcode, new_name, Depth)
##Sumamos el depth de los genes con el mismo nombre que están en un mismo barcode
plasmid_new<-plasmid_new %>% group_by(Barcode,new_name) %>% 
  summarize(sum_depth=sum(Depth))
##A formato ancho
lefse<-plasmid_new %>% 
  pivot_wider(names_from = Barcode, values_from = sum_depth) %>%
  replace(is.na(.), 0)

status = c("status",rep("barcode30", 6), "barcode31", rep("barcode32", 3), NA, "barcode32", NA,
           "barcode31", rep("barcode33", 3), NA, "barcode33", rep(NA,5))
#Añadimos el status y seleccionamos los coassemblies de interés
lefse_new <- rbind(status, lefse)  %>%
  mutate(new_name = gsub("\\s", "_", new_name)) %>%  # Reemplazar espacios por guiones bajos en la columna new_name
  select(new_name, where(~ .[1] %in% c("barcode30", "barcode32")))

write.table(lefse_new, file = "for_lefse_plasmid.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)



##PARA EL COG

plasmid_cog<-read_tsv("merged_cog_results_plasmid.gff",
                  col_names = c("Barcode","Contig","Feature","Start","End","Name", "Eggnog","cog","Depth"))
##Asignamos categorias definitivas
plasmid_new_cog <- plasmid_cog %>%
  mutate(new_name = ifelse(!is.na(Eggnog) & Eggnog != "None", cog, 
                           ifelse(!is.na(Name), Name, Feature))) %>%
  select(Barcode, new_name, Depth)
##Sumamos el depth de los genes con el mismo nombre que están en un mismo barcode
plasmid_new_cog<-plasmid_new_cog %>% group_by(Barcode,new_name) %>% 
  summarize(sum_depth=sum(Depth))
##A formato ancho
lefse_cog<-plasmid_new_cog %>% 
  pivot_wider(names_from = Barcode, values_from = sum_depth) %>%
  replace(is.na(.), 0)

status = c("status",rep("barcode30", 6), "barcode31", rep("barcode32", 3), NA, "barcode32", NA,
           "barcode31", rep("barcode33", 3), NA, "barcode33", rep(NA,5))
#Añadimos el status y seleccionamos los coassemblies de interés
lefse_new_cog <- rbind(status, lefse_cog)  %>%
  mutate(new_name = gsub("\\s", "_", new_name)) %>%  # Reemplazar espacios por guiones bajos en la columna new_name
  select(new_name, where(~ .[1] %in% c("barcode30", "barcode32","barcode33")))
##filas que estén presentes en todos
lefse_filtered_cog <- lefse_new_cog %>%
  filter(rowSums(select(., -new_name) != 0) > 0)

write.table(lefse_filtered_cog, file = "for_lefse_cog_plasmid.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)


###BOXPLOT

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

ggplot(boxlef, aes(x = status, y = depth, fill = status)) +
  # Añadir el boxplot con bordes negros y color según el status
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.6) +
  # Añadir los puntos de jitter sin colores personalizados
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  # Ajustar las etiquetas y el título
  labs(title = "Profundidad por Coassembly para 'Energy production and conversion'",
       x = "Coassembly",
       y = "Depth") +
  # Aplicar un tema minimalista
  theme_minimal() +
  # Eliminar la leyenda
  guides(fill = "none") +  # Alternativamente, puedes usar: theme(legend.position = "none")
  # Modificar las etiquetas del eje Y para asignar "baseline" y "24 horas"
  scale_x_discrete(labels = function(x) ifelse(x == "barcode30", "baseline", 
                                                 ifelse(x == "barcode32", "24 horas",
                                                        ifelse(x=="barcode33","10 días",x)))) +
  # Personalizar el título y las etiquetas
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank()
  )

# Paso 5: Crear el segundo gráfico boxplot (sin 'barcode03')
data_C_no_barcode03 <- boxlef %>% filter(sample != 'barcode03')
ggplot(data_C_no_barcode03, aes(x = status, y = depth, fill = status)) +
  # Añadir el boxplot con bordes negros y color según el status
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.6) +
  # Añadir los puntos de jitter sin colores personalizados
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  # Ajustar las etiquetas y el título
  labs(title = "Profundidad por Coassembly sin 'barcode03' para 'Energy production and conversion'",
       x = "Coassembly",
       y = "Depth") +
  # Aplicar un tema minimalista
  theme_minimal() +
  # Eliminar la leyenda
  guides(fill = "none") +  # Alternativamente, puedes usar: theme(legend.position = "none")
  # Modificar las etiquetas del eje Y para asignar "baseline" y "24 horas"
  scale_x_discrete(labels = function(x) ifelse(x == "barcode30", "Baseline", 
                                                 ifelse(x == "barcode32", "24 horas",
                                                        ifelse(x=="barcode33", "10 días",x)))) +
  # Personalizar el título y las etiquetas
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank()
  )


####GRAFICO DE BARRAS TODAS LAS CATEGORIAS


boxlef<-boxlef %>% rename(category=new_name) %>% filter(category!="hypothetical_protein")

# Cargar las librerías necesarias
library(dplyr)
library(ggplot2)

# 1. Calcular la profundidad mediana por categoría y seleccionar las 8 categorías con mayor mediana
top_categories <- boxlef %>%
  group_by(category) %>%
  summarise(median_depth = median(depth, na.rm = TRUE)) %>%
  arrange(desc(median_depth)) %>%
  slice(1:9) %>%
  pull(category)


# 2. Filtrar los datos para incluir solo las categorías seleccionadas
filtered_data <- boxlef %>%
  filter(category %in% top_categories)

# Diccionario con los nombres reales de las categorías COG
real_names <- c("C" = "Energy production and conversion",
                "E" = "Amino acid transport and metabolism",
                "H" = "Coenzyme transport and metabolism",
                "I" = "Lipid transport and metabolism",
                "J" = "Translation, ribosomal structure and biogenesis",
                "L" = "Replication, recombination and repair",
                "P" = "Inorganic ion transport and metabolism",
                "S" = "Function unknown",
                "M"="Cell wall/membrane/envelope biogenesis")

# Renombrar las categorías en la columna 'category'
filtered_data <- filtered_data %>%
  mutate(category = recode(category, !!!real_names))


ggplot(filtered_data, aes(x = status, y = depth, fill=status)) +
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.6) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  facet_wrap(~ category, scales = "free_y") +
  theme_minimal() +
  labs(title = "Profundidad de cada Categoría antes y después del tratamiento",
       x = NULL,  # No mostrar etiqueta del eje x
       y = "Profundidad media",
       color = "Muestra") +
  scale_x_discrete(labels = function(x) ifelse(x == "barcode30", "Baseline", 
                                               ifelse(x == "barcode32", "24 horas",
                                                      ifelse(x=="barcode33", "10 días",x)))) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 9),  # Reducir el tamaño del texto de las etiquetas del eje x
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 9),  # Reducir el tamaño del texto de los 'facet'
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside"
  ) +
  theme(axis.text.x = element_blank()) + # Esto oculta los textos del eje x
  guides(fill = "none")



data_all_no_barcode03 <- boxlef %>% filter(sample != 'barcode03')
boxlef2<-data_all_no_barcode03

# 1. Calcular la profundidad mediana por categoría y seleccionar las 8 categorías con mayor mediana
top_categories <- boxlef2 %>%
  group_by(category) %>%
  summarise(median_depth = median(depth, na.rm = TRUE)) %>%
  arrange(desc(median_depth)) %>%
  slice(1:9) %>%
  pull(category)


# 2. Filtrar los datos para incluir solo las categorías seleccionadas
filtered_data2 <- boxlef2 %>%
  filter(category %in% top_categories)

# Diccionario con los nombres reales de las categorías COG
real_names <- c("C" = "Energy production and conversion",
                "E" = "Amino acid transport and metabolism",
                "H" = "Coenzyme transport and metabolism",
                "I" = "Lipid transport and metabolism",
                "J" = "Translation, ribosomal structure and biogenesis",
                "L" = "Replication, recombination and repair",
                "P" = "Inorganic ion transport and metabolism",
                "S" = "Function unknown",
                "M"="Cell wall/membrane/envelope biogenesis")

# Renombrar las categorías en la columna 'category'
filtered_data2 <- filtered_data2 %>%
  mutate(category = recode(category, !!!real_names))


ggplot(filtered_data2, aes(x = status, y = depth, fill=status)) +
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.6) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  facet_wrap(~ category, scales = "free_y") +
  theme_minimal() +
  labs(title = "Profundidad de cada Categoría (sin Barcode03) antes y después del tratamiento",
       x = NULL,  # No mostrar etiqueta del eje x
       y = "Profundidad media",
       color = "Muestra") +
  scale_x_discrete(labels = function(x) ifelse(x == "barcode30", "Baseline", 
                                               ifelse(x == "barcode32", "24 horas",
                                                      ifelse(x=="barcode33", "10 días",x)))) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 9),  # Reducir el tamaño del texto de las etiquetas del eje x
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 9),  # Reducir el tamaño del texto de los 'facet'
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside"
  ) +
  theme(axis.text.x = element_blank()) + # Esto oculta los textos del eje x
  guides(fill = "none")



#####GRAFICA DE PESOS


file<-read_tsv("tamaño_barcode.tsv", col_names = c("Tamaño","Muestra"))
modified_data <- file %>%
  mutate(Muestra = str_extract(Muestra, "barcode\\d+")) %>%  
  filter(!(Muestra %in% c(paste0("barcode",30:33),"barcode07","barcode11","barcode13","barcode14","barcode18",paste0("barcode",20:22)))) %>% 
  mutate(status = case_when(
    Muestra %in% c('barcode01', 'barcode02', 'barcode03', 'barcode04', 'barcode05','barcode06') ~ 'barcode30',
    Muestra %in% c('barcode08', 'barcode09', 'barcode10', 'barcode12') ~ 'barcode32'  ,
    TRUE ~ 'barcode33'
  )) 



# Convertir 'Tamaño' a numérico en Gigabytes
modified_data$Tamaño_GB <- as.numeric(sub("G", "", modified_data$Tamaño))
# Asegurar que 'status' es un factor con el orden especificado
modified_data$status <- factor(modified_data$status, levels = c("barcode30", "barcode32", "barcode33"))

# Cargar las librerías necesarias
library(ggplot2)
library(ggpubr)

# Realizar la prueba de Kruskal-Wallis
kruskal_result <- kruskal.test(Tamaño_GB ~ status, data = modified_data)
print(kruskal_result)

# Realizar pruebas de comparaciones por pares con Wilcoxon
pairwise_wilcox <- pairwise.wilcox.test(modified_data$Tamaño_GB, modified_data$status, p.adjust.method = "BH")
print(pairwise_wilcox)

# Definir colores personalizados para cada 'status'
status_colors <- c("barcode30" = "#E69F00", "barcode32" = "#56B4E9", "barcode33" = "#009E73")

# Crear el gráfico
ggplot(modified_data, aes(x = status, y = Tamaño_GB, fill = status)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.2, size = 3, alpha = 0.7, color = "black") +
  scale_fill_manual(values = status_colors) +
  labs(
    title = "Distribución del Tamaño por Estado",
    x = NULL,
    y = "Tamaño (GB)"
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
                                               ifelse(x == "barcode32", "24 horas",
                                                      ifelse(x=="barcode33", "10 días",x))))


