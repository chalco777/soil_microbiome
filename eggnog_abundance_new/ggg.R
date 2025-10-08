plot_box <- gb$filtered_data
  #     dplyr::group_by(!!rlang::sym(group_by_string),status) %>%
  # dplyr::filter(sd(conteo_total) > 0) %>% 
  # dplyr::ungroup()%>%
library(ggplot2)
library(reshape2)
library(pheatmap)


res_sign <- gt$maaslin_results$results %>%
  filter(qval < 0.25)  # Filtramos por p-valor < 0.05
res_sign <- res_sign %>%
  mutate(
    qval_adj = if_else(qval == 0, 1e-300, qval),  # Evitar log10(0)
    heat_value = -log10(qval_adj) * sign(coef)
  )
heat_data <- res_sign %>%
  select(feature, value, heat_value) %>%
  pivot_wider(
    names_from = value,
    values_from = heat_value
  )

heat_mat <- as.data.frame(heat_data[,-1])  # Todas las columnas menos la 1
rownames(heat_mat) <- heat_data$feature
sig_features    <- rownames(heat_mat)  # features significativos

library(grid)

pheatmap(
  heat_mat,
  color = colorRampPalette(c("blue", "white", "red"))(100),  # Gradiente azul-blanco-rojo
  main = "Heatmap of Log-Scaled Significance of Features with Baseline as Reference",
  cluster_rows = FALSE,   # Sin agrupar filas
  cluster_cols = FALSE,   # Sin agrupar columnas
  fontsize_row = 9,       # Tamaño de letra para las filas
  fontsize_col = 11,       # Tamaño de letra para las columnas
  cellwidth = 25,         # Ancho de las celdas
  cellheight = 20,        # Altura de las celdas
    border_color = "black",       # Sin bordes entre celdas
  breaks = seq(-max(abs(heat_mat)), max(abs(heat_mat)), length.out = 101)  # Escala centrada en 0,

)


gtt

# 3A. Filtrar 'expanded_data' sólo con esos features significativos
expanded_filt <- gt$expanded_data %>%
  dplyr::filter(pathways%in% sig_features) %>%
  mutate(name=ifelse(status == "barcode30", "Baseline", 
                     ifelse(status == "barcode32", "24 horas",
                            ifelse(status=="barcode33", "10 días",status)))) %>% 
  dplyr::filter(name%in% c(res_sign$value,"Baseline"))  %>% group_by(pathways,name) %>%
  summarize(mean_value=mean(conteo_total)) %>% 
  pivot_wider(names_from=name,values_from = mean_value) %>%  ungroup() %>% 
  mutate(across(-c(pathways,Baseline),~ log2( (.x + 1) / (Baseline + 1) ), .names = "{.col}")) %>% 
  column_to_rownames(var = "pathways") %>% select(-Baseline)
 

print("ACA EXPANDED")
print("ACA FILT")
# 3B. Recodificar nombres en expanded_filt
old_feats_2      = unique(expanded_filt[[annotation_col]])
recoded_vec_2    = recode_feature_names(old_feats_2, annotation_type)

expanded_filt <- expanded_filt %>%
  dplyr::mutate(
    !!annotation_col := dplyr::recode(
      .data[[annotation_col]],
      !!!recoded_vec_2
    )
  )
range_val <- max(abs(expanded_filt), na.rm = TRUE)

# 3D. Crear y guardar el heatmap de log2(abund)
pheatmap::pheatmap(
  expanded_filt,
  color         = colorRampPalette(c("blue", "white", "red"))(100),
  main          = "Heatmap of log2(Abundance) of Features (Qval < 0.25) with Baseline as Reference",
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  fontsize_row  = 9,
  fontsize_col  = 11,
  cellwidth     = 25,
  cellheight    = 20,
  border_color  = "black",
  breaks = seq(-range_val, range_val, length.out = 101)  
  #filename      = paste0("heatmap_log2abund_", annotation_type, ".png")
)
library(grid)

# 1. Abres el dispositivo (PNG, PDF, etc.)
expanded_filt$feature <- rownames(expanded_filt)
df_long <- expanded_filt %>%
  pivot_longer(
    cols      = -feature,
    names_to  = "status",
    values_to = "value"
  )

p <- ggplot(df_long, aes(x = status, y = feature, fill = value)) +
  geom_tile(color = "black") +  # Borde negro en cada celda
  scale_fill_gradient2(
    midpoint = 0,   # Centra la escala en 0
    low = "blue", mid = "white", high = "red",
    name = "log2FC"
  ) +
  labs(
    title = "Heatmap of Log-Scaled Significance of Features with Baseline as Reference",
    x = "Condición",
    y = "Feature"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    # Centrar y dar formato al título
    plot.title = element_text(
      hjust  = 0.5,  # Centra el título horizontalmente
      face   = "bold",
      size   = 14
    ),
    # Ajustar el tamaño de las etiquetas del eje X (columnas) y eje Y (filas)
    axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x =element_text(size = 14),
    axis.title.y=element_text(size = 14),
    # Quitar las líneas de la cuadrícula
    panel.grid  = element_blank()
  ) +
  # Reducir el espacio vacío alrededor de las celdas
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  # Mantener la relación de aspecto fijo (hace que las celdas sean cuadradas).
  # Si lo quieres rectangular, ajusta 'ratio'
  coord_fixed()

# Ahora guardar con un tamaño más pequeño para que las celdas queden más chicas
ggsave(
  filename = "my_heatmap_multi_condition.png",
  plot     = p,
  width    = 10,   # Ajusta el ancho (en pulgadas si no se especifican units)
  height   = 8,   # Ajusta el alto
  dpi      = 300,
  bg       = "white"
)
