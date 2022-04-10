
# To Do List --------------------------------------------------------------

#' Agregar estadígrafos a los gráficos de diversidad alfa y beta
#' Mejorar pipeline de redes de co-ocurrencia



# Cargamos librerias ------------------------------------------------------

.packages <- c("microbiome",
               "microbiomeutilities",
               "dplyr",
               "tibble",
               "vegan",
               "RColorBrewer", 
               "viridis",
               "ggplot2",
               "ggpubr",
               "grid",
               "ranacapa", 
               "patchwork", 
               "DESeq2",
               "VennDiagram",
               "UpSetR",
               "cooccur",
               "visNetwork")

.to_install <-!.packages %in% installed.packages()


if(any(.to_install)){
        install.packages(.packages[.to_install])
}

sapply(.packages, require, character.only = T)




# Importamos los datos ----------------------------------------------------

data <- readRDS("data/ps.pool.ind.solo.bacte.rds")




# Limpieza de los datos ---------------------------------------------------


# 1) Cambiamos los nombres de las taxas: desde las secuencias a nombres en formatos "ASV+numero"
taxa_names(data) <- paste0("ASV", seq(ntaxa(data)))


# 2) Eliminamos las reads que:
#' no estan identificadas al nivel de phylum,
#' que son identificadas como Cloroplastos,
#' Mitocondrias o
#' Eukaryotas
data_clean <- subset_taxa(data, !Phylum %in% c("", "uncharacterized") &
                                !Order == "Chloroplast" &
                                !Family == "Mitochondria" &
                                !Kingdom == "Eukaryota")


# 3) Eliminamos las muestras agregadas (type == "Ind", en los metadatos)
data_clean <- subset_samples(data_clean, type == "Ind")


# 4) Nos quedamos con las muestras control, infectadas, y las infectadas + DH97 (probioticos)
my_data <- subset_samples(data_clean, condition %in% c("Ctrl", "Va", "DH97", "Va+DH97"))


cat("Número de taxas sin filtrar : ", ntaxa(data), 
    "\nNúmero de taxas después de filtrar : ", ntaxa(my_data))




# Rarefaccion -------------------------------------------------------------

rarefaction_curve <- ggrare(my_data, step = 100, color = "condition", label = "subject", plot = FALSE)

rarefaction_curve + 
        labs(title = "Curvas de rarefaccion para las muestras") +
        xlab(label = "Numero de secuencias en la muestra") +
        ylab(label = "Riqueza de especies") + 
        theme_bw() +
        scale_fill_viridis(discrete = TRUE) +
        scale_color_viridis(discrete = TRUE) 


my_data_rarefied <- rarefy_even_depth(my_data, 
                                      rngseed=01021197, 
                                      sample.size=0.9*min(sample_sums(my_data)), 
                                      replace=F)




# Graficos de abundancia a distintos niveles taxonomicos ------------------


# 1) Graficamos a nivel de phylum
data_phylum <- my_data %>% 
        aggregate_taxa(level = "Phylum") %>% 
        microbiome::transform("compositional")

plot_composition(data_phylum,
                 group_by = "condition") +
        labs(title = "Abundancias relativas a nivel de Phylum",
             fill = "Taxonomias") +
        ylab(label = "Abundancia relativa (%)") +
        xlab(label = "Muestras") +
        scale_y_continuous(labels = c("0%","25%","50%","75%","100%"), expand = c(0, 0)) +
        coord_cartesian(clip = "off") +
        theme_bw(base_line_size = 0, base_rect_size = 1) +
        theme(plot.title = element_text(face = "bold", size = 14),
              legend.title = element_text(face = "bold"),
              strip.background = element_rect(fill = "grey30"),
              strip.text = element_text(color = "white", face = "bold")) +
        scale_fill_brewer(type = "qual", palette = "Spectral")



# 2) Graficamos a nivel de clase
data_class <- my_data %>% 
        aggregate_taxa(level = "Class") %>% 
        microbiome::transform("compositional")

plot_composition(data_class,
                 group_by = "condition") +
    labs(title = "Abundancias relativas a nivel de Class",
         fill = "Taxonomias") +
    ylab(label = "Abundancia relativa (%)") +
    xlab(label = "Muestras") +
    scale_y_continuous(labels = c("0%","25%","50%","75%","100%"), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_bw(base_line_size = 0, base_rect_size = 1) +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.title = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey30"),
          strip.text = element_text(color = "white", face = "bold")) +
    scale_fill_brewer(type = "qual", palette = "Spectral")



# 3) Graficamos a nivel de orden     
data_order <- my_data %>% 
        aggregate_taxa(level = "Order") %>% 
        microbiome::transform("compositional")

my_cols <- c(brewer.pal(n=11,"PiYG")[-6], 
             brewer.pal(n=11,"BrBG")[-6], 
             brewer.pal(n=11, "PRGn")[c(1:5,7:8)],
             brewer.pal(n=11, "RdYlBu")[1:4],
             brewer.pal(n=8, "Dark2")[c(6,8)])

plot_composition(data_order,
                 group_by = "condition") +
    labs(title = "Abundancias relativas a nivel de Order",
         fill = "Taxonomias") +
    ylab(label = "Abundancia relativa (%)") +
    xlab(label = "Muestras") +
    scale_y_continuous(labels = c("0%","25%","50%","75%","100%"), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_bw(base_line_size = 0, base_rect_size = 1) +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.title = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey30"),
          strip.text = element_text(color = "white", face = "bold")) +
    scale_fill_manual(values = my_cols)



# 4) Graficamos a nivel de familia

data_family <- my_data %>% 
  microbiome::transform("compositional") %>% 
  aggregate_top_taxa(top = 14, level = "Family")


my_cols <- c(brewer.pal(11, "Spectral")[1:8],
             brewer.pal(8, "Dark2")[8],
             brewer.pal(11, "Spectral")[9:11],
             brewer.pal(8, "Dark2")[4:6])


plot_composition(data_family,
                 group_by = "condition") +
  labs(title = "Abundancias relativas a nivel de Order",
       fill = "Taxonomias") +
  ylab(label = "Abundancia relativa (%)") +
  xlab(label = "Muestras") +
  scale_y_continuous(labels = c("0%","25%","50%","75%","100%"), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_bw(base_line_size = 0, base_rect_size = 1) +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey30"),
        strip.text = element_text(color = "white", face = "bold")) +
  scale_fill_manual(values = my_cols)



# 4) Nota: Agrupar taxas de forma manual ----------------------------------

# En ese momento no sabia que se podia agrupar por top_taxa, asi que lo hice manual:

# # Identificamos las top 14 familias
# top_14_families <- my_data %>% aggregate_taxa(level = "Family") %>% top_taxa(n=14)
# 
# 
# # Tenemos todos los nombres de las familias en un vector. Si alguna de esas familias no es top, lo convertimos a "Others"
# all_families <- tax_table(my_data)[,"Family"] %>% as.data.frame() %>% pull("Family") %>% as.character()
# new_family_names <- ifelse(all_families %in% top_14_families, yes = all_families, no = "Others")
# 
# 
# # Agrupamos por familias y sustituimos el nombre de las taxas que no estan en el top por el nombre "Others"
# data_family <- my_data
# tax_table(data_family)[, "Family"] <- new_family_names
# 
# # Transformamos los datos a composicional y graficamos
# data_family <- data_family %>% microbiome::transform("compositional")
# 
# plot_composition(data_family)
# 
# my_cols <- c(brewer.pal(11, "Spectral")[1:8],
#              brewer.pal(8, "Dark2")[8],
#              brewer.pal(11, "Spectral")[9:11],
#              brewer.pal(8, "Dark2")[4:6])
# 
# ggplot(data = psmelt(data_family), 
#        aes(x = Sample,
#            y = Abundance,
#            fill = Family)) + 
#     geom_bar(position = "stack", stat = "identity") +
#     facet_wrap(facets = ~ condition, nrow = 1, scales = "free") +
#     labs(title = "Abundancias relativas a nivel de Family",
#          fill = "Taxonomias") +
#     ylab(label = "Abundancia relativa (%)") +
#     xlab(label = "Muestras") +
#     scale_y_continuous(labels = c("0%","25%","50%","75%","100%"), expand = c(0, 0)) +
#     coord_cartesian(clip = "off") +
#     theme_bw(base_line_size = 0, base_rect_size = 1) +
#     theme(plot.title = element_text(face = "bold", size = 14),
#           legend.title = element_text(face = "bold"),
#           strip.background = element_rect(fill = "grey30"),
#           strip.text = element_text(color = "white", face = "bold"),
#           strip.text.y = element_blank()) +
#     scale_fill_manual(values = my_cols) 




# Diversidad alfa ---------------------------------------------------------

# Ploteamos los graficos de diversidad y le agregamos las pruebas estadisticas
plot_richness(my_data_rarefied, x="condition", measures = c("Chao1", "Shannon"), color = "condition") +
    labs(title = "Indices de diversidad alfa",
         color = "Tratamientos") +
    ylab(label = "Indice de diversidad alfa") +
    xlab(label = "Tratamientos") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.title = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey30"),
          strip.text = element_text(color = "white", face = "bold")) +
    scale_color_viridis(discrete = TRUE)
    

# Hacemos pruebas estadisticas para inferir diferencias

#' Calculamos las dos medidas de diversidad que ploteamos anteriormente.
#' 
#' Se uso la funcion estimate_richness porque esa es la funcion que plot_richeness llama por detras
alpha_diversities <- estimate_richness(my_data_rarefied, measures = c("Chao1", "Shannon"))


#' A la tabla con diversidades le agregamos una columna con la informacion de los tratamientos 
alpha_diversities$condition <- meta(my_data_rarefied)[,"condition"]


#' Realizamos una prueba Kruskall walis para cada indice de diversidad alpha
kw_chao <- kruskal.test(data = alpha_diversities, Chao1 ~ condition)
kw_chao

kw_shannon <- kruskal.test(data = alpha_diversities, Shannon ~ condition)
kw_shannon




# Diversidad Beta ---------------------------------------------------------
pca_bray <- ordinate(my_data_rarefied, method = "MDS", distance = "bray")

p.bray <- plot_ordination(my_data, ordination = pca_bray, color = "condition") + 
    geom_point(size = 3) +
    labs(title = "Diversidad beta por Bray-Curtis",
         color = "Tratamientos") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          legend.title = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey30"),
          strip.text = element_text(color = "white", face = "bold")) +
    scale_color_viridis(discrete = TRUE) +
    stat_ellipse(type = "norm")



set.seed(01021997)
pca_wunifrac <- ordinate(my_data, method = "MDS", distance = "wunifrac")

p.wuni <- plot_ordination(my_data, ordination = pca_wunifrac, color = "condition") + 
    geom_point(size = 3) +
    labs(title = "Diversidad beta por Weighted Unifrac",
         color = "Tratamientos") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          legend.title = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey30"),
          strip.text = element_text(color = "white", face = "bold"),
          legend.position = "none") +
    scale_color_viridis(discrete = TRUE) +
    stat_ellipse(type = "norm")


ggarrange(p.bray, p.wuni, nrow = 1, labels = c("A", "B"), common.legend = TRUE)




# Permanovas de diversidad beta:
set.seed(1021997)
permanova_bray <- adonis(phyloseq::distance(my_data_rarefied,
                                            method="bray", 
                                            weighted=T) ~ condition, 
                         data = meta(my_data_rarefied))

summary(permanova_bray)
permanova_bray



set.seed(1021997)
permanova_wunif <- adonis(phyloseq::distance(my_data_rarefied, 
                                             method="unifrac", 
                                             weighted=T) ~ condition, 
                          data = meta(my_data_rarefied))

summary(permanova_wunif)
permanova_wunif


# Otra forma de visualizar la beta-diversidad: Heatmaps
data_phylum_filtered <- prune_taxa(taxa = taxa_sums(data_phylum) > 0, data_phylum)

plot_heatmap(data_phylum_filtered, 
             "PCoA", 
             "bray", 
             "condition", 
             "Phylum", 
             sample.order = "condition") +
  theme_biome_utils() +
  theme(axis.text.x = element_text(angle = 90))



# Core microbiota ---------------------------------------------------------

# Transformamos los datos a datos composicionales  y nos quedamos con las taxas que tienen más de 0 conteos
my_data_compositional <- my_data %>% microbiome::transform("compositional") 

my_data_compositional <-  prune_taxa(taxa = taxa_sums(my_data_compositional) > 0, x = my_data_compositional)


# Generamos un vector con las condiciones que fueron estudiadas
conditions <- my_data_compositional %>% psmelt() %>% select("condition") %>% unique() %>% pull() %>% as.character()


# Le cambiamos los nombres de las taxas desde el formato ASVX al formato ASVX:nivel_taxonomico
my_data_compositional_formated <- format_to_besthit(my_data_compositional)


list_core <- c() # Un vector vacio para guardar la infomacion


for (n in conditions){ # para cada condicion en el vector conditions
    
    # Nos quedamos con las muestras que son de la condicion n
    ps.sub <- subset_samples(my_data_compositional_formated, condition == n) 
    
    # Obtenemos los miembros de la core microbiota que estan en al menos el 75% de las muestras con un 0.1% de abundancia
    core_m <- core_members(ps.sub, 
                           detection = 0.001, 
                           prevalence = 0.75)
    
    # Output para saber cuantas taxas cumplen esos criterios
    print(paste0("No. de taxas de la microbiota core para la condicion ", n, " : ", length(core_m)))
    list_core[[n]] <- core_m 
    #print(list_core)
}


# Revisamos la lista que tiene todas las core microbiota de cada uno de los tratamientos
list_core


# Graficamos con Venn: No se imprime en pantalla, pues hay que verlo en la carpeta "plots"
venn.diagram(x= list_core,
             category.names = conditions,
             disable.logging = T,
             filename = "plots/venn.tiff",
             imagetype = "tiff",
             output = T,
             
             # Title
             main = "Core microbiota by condition",
             main.cex = 3,
             
             # Circles
             lty = "blank",
             fill = viridis(4),
             alpha = rep("0.4", 4),
             
             # Text within circles
             fontface = "italic",
             cex = 2,
             
             # Text of labels
             cat.fontface = "bold",
             cat.cex = 2)



# Graficamos con Upset Plot
upset(data = fromList(list_core),
      nintersects = NA,
      empty.intersections = "on",
      
      order.by = "freq",
      
      point.size = 3)


# Screening de la core microbiota (Se puede encontrar en https://microbiome.github.io/tutorials/Core.html)

#' NOTA: Se muestran la core microbiotas comun a todas las condiciones. 
#' Tambien se pueden hacer por separado (un phyloseq con la microbiota de cada condicion)

# 1) Core Line Plots
detections <- c(0, 0.1, 0.5, 2, 5, 20)/100 # Representa la abundancia de cada 
prevalences <- seq(from = 0.05, to = 1, by = 0.05)


plot_core(my_data_compositional,
          prevalences = prevalences,
          detections = detections,
          plot.type = "lineplot") + 
  xlab(label = "Relative abundances (%)") +
  ggtitle(label = "Core microbiota (Line plot)") +
  theme_biome_utils()


# 2) Core Heat Maps
prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-2), log10(.2), length = 10), 3)

plot_core(data_phylum,
          plot.type = "heatmap",
          colours = brewer.pal(5, "Spectral"),
          prevalences = prevalences,
          detections = detections, 
          min.prevalence = .5) +
  theme_biome_utils() +
  ggtitle(label = "Core microbiota (Heat map)") + 
  xlab("Detection Threshold (Relative Abundance (%))")



# Expresion diferencial ---------------------------------------------------

# 1) Agrupamos por generos
data_genus <- subset_taxa(my_data, !is.na(Genus))
data_genus %>% aggregate_taxa("Genus")


# 2) Convertimos el objeto phyloseq a un objeto DESeq2
my_deseq <- phyloseq_to_deseq2(data_genus, ~ condition)


# 3) Realizamos el analisis de expresion diferencial
dds <- DESeq(my_deseq, test = "Wald", fitType = "local", sfType = "poscounts")


# 4) Observamos los resultados
dh_contrast <- results(dds, contrast=c("condition", "Ctrl", "DH97"))
va_contrast <- results(dds, contrast=c("condition", "Ctrl", "Va"))
dh.va_contrast <- results(dds, contrast=c("condition", "Ctrl", "Va+DH97"))
va_dh.va_contrast <- results(dds, contrast = c("condition", "Va", "Va+DH97"))


# 5) Convertimos los resultados en un dataframe, eliminamos NA y valores no significativos
dh_df <- as.data.frame(dh_contrast) %>% rownames_to_column(var = "ID") %>% filter(padj < 0.05)
va_df <- as.data.frame(va_contrast) %>% rownames_to_column(var = "ID") %>% na.omit() %>% filter(padj < 0.05)
dh.va_df <- as.data.frame(dh.va_contrast) %>% rownames_to_column(var = "ID") %>% na.omit() %>% filter(padj < 0.05)
va_dh.va_df <- as.data.frame(va_dh.va_contrast) %>% rownames_to_column(var = "ID") %>%  na.omit() %>% filter(padj < 0.05)


# 6) Agregamos la clasificacion taxonomica a los ASV que fueorn identificados con abundacia diferencial

#' Si revisamos los objetos recien creados (dh_df, por ejemplo), vemos que las taxonomias son identificadas como ASV
#' Queremos que, en lugar de ASV, sean identificados por su taxonomia a nivel de Genus
#' Tenemos que hacer que la tax_table (que tiene las taxonomias) sea un dataframe y que tenga una columna con los ASV correspondientes
tax_table_genus <- as.data.frame(tax_table(data_genus)) %>% rownames_to_column(var = "ID")
tax_table_genus$ID <- paste0("ASV", seq(1:nrow(tax_table_genus)))


#' Ahora podemos unir las tablas con la informacion de la expresion diferencial y la tabla con taxonomias usando las ASV en comun de ambas tablas
dh_df <- left_join(dh_df, tax_table_genus, by = "ID")
va_df <- left_join(va_df, tax_table_genus, by = "ID")
dh.va_df <- left_join(dh.va_df, tax_table_genus, by = "ID")
va_dh.va_df <- left_join(va_dh.va_df, tax_table_genus, by = "ID")


# 6) Graficamos los taxas que presentan abundancias diferenciadas

my_qual_cols <- c(brewer.pal(n=12, "Set3")[-2],
                 brewer.pal(n=8, "Dark2")[6:8])


dh_da_p <- ggplot(dh_df, aes(x = Genus, y = log2FoldChange, color = Genus)) +
    geom_point(size = 4, alpha = 0.6) +
    labs(y = "\nLog2 Fold-Change for DH97 vs. Controls", x = "") +
    theme_bw() +
    theme(axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 12),
          axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.position = "none") +
    coord_flip() +
    geom_hline(yintercept = 0, linetype="dotted") +
    scale_color_manual(values = my_qual_cols)



va_da_p <- ggplot(va_df, aes(x = Genus, y = log2FoldChange, color = Genus)) +
    geom_point(size = 4, alpha = 0.6) +
    labs(y = "\nLog2 Fold-Change for Va vs. Controls", x = "") +
    theme_bw() +
    theme(axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 12),
          axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.position = "none") +
    coord_flip() +
    geom_hline(yintercept = 0, linetype="dotted") +
    scale_color_manual(values = my_qual_cols)



dh.va_da_p <- ggplot(dh.va_df, aes(x = Genus, y = log2FoldChange, color = Genus)) +
    geom_point(size = 4, alpha = 0.6) +
    labs(y = "\nLog2 Fold-Change for DH97+Va vs. Controls", x = "") +
    theme_bw() +
    theme(axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 12),
          axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.position = "none") +
    coord_flip() +
    geom_hline(yintercept = 0, linetype="dotted") +
    scale_color_manual(values = my_qual_cols)



va_dh.va_da_p <- ggplot(va_dh.va_df, aes(x = Genus, y = log2FoldChange, color = Genus)) +
    geom_point(size = 4, alpha = 0.6) +
    labs(y = "\nLog2 Fold-Change for Va vs. DH97+Va", x = "") +
    theme_bw() +
    theme(axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 12),
          axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.position = "none") +
    coord_flip() +
    geom_hline(yintercept = 0, linetype="dotted") +
    scale_color_manual(values = my_qual_cols)


(dh_da_p + va_da_p) / (dh.va_da_p + va_dh.va_da_p) + labs(caption = "Numeros positivos implican un aumento respecto al control")

ggarrange(dh_da_p, va_da_p, dh.va_da_p, va_dh.va_da_p, 
          nrow = 2, 
          ncol = 2, 
          common.legend = T,
          labels = c("A", "B", "C", "D"))



# Redes de coocurrencia ---------------------------------------------------

# 1) Agregamos una columna con las familias a la otu_table.

#' Creamos un vector con todas las familias
all_families <- my_data %>% subset_samples(condition == "Ctrl") %>%  tax_table() %>% as.data.frame() %>% pull("Family") %>% as.character()


#' Extraemos la otu_table de mis datos
otu_table <- my_data %>% subset_samples(condition == "DH97") %>% otu_table() %>% t() %>% as.data.frame()


#' Le agregamos una columna con las familias a la otu_table
otu_table$Family <- all_families


#' Colapsamos todas cada una de las observaciones de una misma familia en una sola row 
family_counts <- aggregate(. ~ Family, data = otu_table, FUN=sum) 


#' Convertimos la columna Family a rownames y luego eliminamos esa columna del df
family_counts <- column_to_rownames(family_counts, var="Family")


# 2) Calculamos las co-ocurrencias:

#' Con esta metodologia para determinar co-ocurrencias nos importa saber si una taxa esta presente
#' o ausente, no nos interesan los numeros de counts brutos.
#' 
#' Esta forma de determinar co-ocurrencia se basa en que si las taxas estan presentes de forma
#' conjunta, entonces se determina la probabilidad de que eso haya sido al azar. Si es baja, entonces
#' se determinan como taxas co-ocurrentes


#' Convertimos los counts a expresiones binarias: 1 si esta presenta; 0 si ausente
family_counts[family_counts > 0] = 1 


#' Calculamos la co-ocurrencia
co_ocurrencias <- cooccur(family_counts, spp_names = TRUE)
co_ocurrencias_df <- co_ocurrencias$results


#' Calculamos el tamaño de los nodos (como el numero de conexiones que tiene con otras strains)
nodes <- data.frame(id = 1:nrow(family_counts),
                    label = rownames(family_counts),
                    color = "#606482",
                    shadow = TRUE)

#' Creamos una tabla con las uniones entre nodos
edges <- data.frame(from = co_ocurrencias_df$sp1, to = co_ocurrencias_df$sp2,
                    color = ifelse(co_ocurrencias_df$p_gt <= 0.05, yes = "#B0B2C1", no = "3C3F51"),
                    dashes = ifelse(co_ocurrencias_df$p_gt <= 0.05, yes = TRUE, no = FALSE))

#' Ploteamos
network <- visNetwork(nodes = nodes, edges = edges) %>% 
    visIgraphLayout(layout = "layout_with_kk")

network

#



# Vias metabolicas --------------------------------------------------------

#' Ya esta hecho 
