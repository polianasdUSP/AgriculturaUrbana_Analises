# 2. Função para correlação + p-valor
calc_cor_p <- function(x, y) {
  valid <- complete.cases(x, y)
  if (sum(valid) > 2) {
    res <- suppressWarnings(cor.test(x[valid], y[valid]))
    return(c(correlation = as.numeric(res$estimate), p_value = res$p.value))
  } else {
    return(c(correlation = NA_real_, p_value = NA_real_))
  }
}

# 3. Criar lista para armazenar os resultados por variável
cor_list_diet <- list()

# 4. Loop pelas variáveis de dieta
for (var in diet_vars) {
  temp <- metadata_filtered %>%
    select(Sample.id, !!sym(var)) %>%
    left_join(SVs_long, by = "Sample.id") %>%
    group_by(Feature.ID) %>%
    summarise(
      correlation = calc_cor_p(Abundance, .data[[var]])["correlation"],
      p_value     = calc_cor_p(Abundance, .data[[var]])["p_value"]
    ) %>%
    ungroup() %>%
    mutate(variable = var)
  
  cor_list_diet[[var]] <- temp
}

# 5. Unir todos os resultados
cor_all_diet <- bind_rows(cor_list_diet)

# 6. Corrigir p-valor (FDR) e marcar significância
cor_all_diet <- cor_all_diet %>%
  group_by(variable) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  ungroup() %>%
  mutate(star = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01  ~ "**",
    p_adj < 0.05  ~ "*",
    TRUE ~ ""
  ))

# 7. Visualizar resultados ordenados
cor_all_diet %>%
  filter(p_adj < 0.05) %>%
  arrange(p_adj)

library(pheatmap)
library(tidyr)
library(dplyr)
library(tibble)

# 1. Filtrar ASVs com pelo menos uma correlação significativa com variáveis de dieta
asvs_signif_diet <- cor_all_diet %>%
  filter(p_adj < 0.05) %>%
  pull(Feature.ID) %>%
  unique()

# 2. Matriz de correlações
cor_matrix_diet <- cor_all_diet %>%
  filter(Feature.ID %in% asvs_signif_diet) %>%
  select(Feature.ID, variable, correlation) %>%
  pivot_wider(names_from = variable, values_from = correlation) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 3. Matriz de significância (*)
star_matrix_diet <- cor_all_diet %>%
  filter(Feature.ID %in% asvs_signif_diet) %>%
  select(Feature.ID, variable, star) %>%
  pivot_wider(names_from = variable, values_from = star) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 1. Remover NAs ou strings vazias das anotações antes de gerar as cores
annotation_row_diet_clean <- annotation_row_diet %>%
  mutate(across(everything(), ~ replace_na(., "Unknown"))) %>%
  mutate(across(everything(), ~ ifelse(. == "", "Unknown", .)))

# 2. Gerar as paletas de cores com número seguro de categorias (limitado pelas paletas)
library(RColorBrewer)

# Garantir que o número de categorias não ultrapasse o máximo suportado
phylum_colors <- brewer.pal(min(length(unique(annotation_row_diet$Phylum)), 12), "Set3")
family_colors <- brewer.pal(min(length(unique(annotation_row_diet$Family)), 12), "Paired")
genus_colors  <- brewer.pal(min(length(unique(annotation_row_diet$Genus)), 8), "Dark2")

ann_colors_diet <- list(
  Phylum = setNames(phylum_colors, unique(annotation_row_diet$Phylum)[1:length(phylum_colors)]),
  Family = setNames(family_colors, unique(annotation_row_diet$Family)[1:length(family_colors)]),
  Genus  = setNames(genus_colors,  unique(annotation_row_diet$Genus)[1:length(genus_colors)])
)

# 3. Repassar para o pheatmap
e <- pheatmap(cor_matrix_diet,
              display_numbers = star_matrix_diet,
              annotation_row = annotation_row_diet,
              annotation_colors = ann_colors_diet,
              clustering_method = "ward.D2",
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              color = color_palette,
              breaks = breaks,
              main = "Correlação entre ASVs e Variáveis da Dieta\n(* FDR < 0.05)",
              fontsize_row = 16,
              fontsize_col = 16,
              fontsize_number = 16,
              border_color = NA)



ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_Significant_correlation_ASVsvsDiet.png", e,  width = 35, height = 20, dpi = 300)


# Preparar nomes taxonômicos amigáveis para as ASVs
taxonomy_named <- taxonomy %>%
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right", extra = "drop") %>%
  mutate(across(Domain:Species, ~ gsub("^[a-z]__*", "", .))) %>%
  mutate(Taxon_Label = paste(Phylum, Family, Genus, sep = " | ")) %>%
  select(Feature.ID, Taxon_Label)

# Preparar dados com nome taxonômico
cor_all_health_named <- cor_all %>%
  left_join(taxonomy_named, by = "Feature.ID") %>%
  filter(!is.na(Taxon_Label))

# Selecionar ASVs significativas
taxa_signif_health <- cor_all_health_named %>%
  filter(p_adj < 0.05) %>%
  pull(Taxon_Label) %>%
  unique()

# Matriz de correlação
cor_matrix_health <- cor_all_health_named %>%
  filter(Taxon_Label %in% taxa_signif_health) %>%
  select(Taxon_Label, variable, correlation) %>%
  pivot_wider(names_from = variable, values_from = correlation, values_fn = mean) %>%
  column_to_rownames("Taxon_Label") %>%
  as.matrix()


# Matriz de asteriscos
star_matrix_health <- cor_all_health_named %>%
  filter(Taxon_Label %in% taxa_signif_health) %>%
  select(Taxon_Label, variable, star) %>%
  pivot_wider(names_from = variable, values_from = star, values_fn = ~ first(na.omit(.))) %>%
  column_to_rownames("Taxon_Label") %>%
  as.matrix()


# Heatmap
e2 <- pheatmap(cor_matrix_health,
               display_numbers = star_matrix_health,
               clustering_method = "ward.D2",
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               color = colorRampPalette(c("blue", "white", "red"))(100),
               breaks = seq(-0.5, 0.5, length.out = 101),
               main = "Taxons  vs Diet (n = 105) \n(FDR < 0.05)",
               fontsize_row = 7,
               fontsize_col = 9,
               fontsize_number = 10,
               border_color = NA)

ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_Significant_correlation_TaxonvsDiet.png", e2,  width = 35, height = 20, dpi = 300)

# Substituir IDs por taxonomia
cor_all_diet_named <- cor_all_diet %>%
  left_join(taxonomy_named, by = "Feature.ID") %>%
  filter(!is.na(Taxon_Label))



#P-Adjusted
# Selecionar apenas taxas com pelo menos uma correlação significativa
taxa_signif_diet <- cor_all_diet_named %>%
  filter(p_adj < 0.05) %>%
  pull(Taxon_Label) %>%
  unique()

# Matriz de correlação (média para duplicatas)
cor_matrix_diet <- cor_all_diet_named %>%
  filter(Taxon_Label %in% taxa_signif_diet) %>%
  select(Taxon_Label, variable, correlation) %>%
  pivot_wider(names_from = variable, values_from = correlation, values_fn = mean) %>%
  column_to_rownames("Taxon_Label") %>%
  as.matrix()

# Matriz de asteriscos (pegar o primeiro não nulo)
star_matrix_diet <- cor_all_diet_named %>%
  filter(Taxon_Label %in% taxa_signif_diet) %>%
  select(Taxon_Label, variable, star) %>%
  pivot_wider(names_from = variable, values_from = star, values_fn = ~ first(na.omit(.))) %>%
  column_to_rownames("Taxon_Label") %>%
  as.matrix()

# Gerar o heatmap
e2<- pheatmap(cor_matrix_diet,
              display_numbers = star_matrix_diet,
              clustering_method = "ward.D2",
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              color = colorRampPalette(c("blue", "white", "red"))(100),
              breaks = seq(-0.5, 0.5, length.out = 101),
              main = "Taxons  vs Diet (n = 105) \n(FDR < 0.05)",
              fontsize_row = 16,
              fontsize_col = 16,
              fontsize_number = 16,
              border_color = NA)


ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_Significant_correlation_TaxonvsDiet.png", e2,  width = 35, height = 20, dpi = 300)



### SEM FDR


library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# 1. Filtrar ASVs com pelo menos uma correlação significativa com p-valor bruto < 0.05
asvs_signif_rawp_df <- cor_all_diet %>%
  filter(p_value < 0.05) %>%
  select(Feature.ID, variable, correlation, p_value, p_adj, star)

# 2. Matriz de correlações (Feature.ID nas linhas, variáveis de dieta nas colunas)
cor_matrix_rawp <- asvs_signif_rawp_df %>%
  select(Feature.ID, variable, correlation) %>%
  pivot_wider(names_from = variable, values_from = correlation) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 3. Matriz de asteriscos com base nos p-VALUES não ajustados
star_matrix_rawp <- cor_all_diet %>%
  filter(Feature.ID %in% rownames(cor_matrix_rawp)) %>%
  select(Feature.ID, variable, p_value) %>%
  mutate(star = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ ""
  )) %>%
  select(Feature.ID, variable, star) %>%
  pivot_wider(names_from = variable, values_from = star) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 4. Anotação taxonômica (Phylum, Family, Genus)
annotation_row_rawp <- taxonomy %>%
  filter(Feature.ID %in% rownames(cor_matrix_rawp)) %>%
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right", extra = "drop") %>%
  mutate(across(Domain:Species, ~ gsub("^[a-z]__*", "", .))) %>%
  select(Feature.ID, Phylum, Family, Genus) %>%
  column_to_rownames("Feature.ID")

# Corrigir problema da paleta para mais categorias
generate_colors <- function(n, palette = "Set3") {
  colorRampPalette(brewer.pal(min(8, brewer.pal.info[palette, ]$maxcolors), palette))(n)
}

ann_colors_rawp <- list(
  Phylum = setNames(generate_colors(length(unique(annotation_row_rawp$Phylum)), "Set3"),
                    unique(annotation_row_rawp$Phylum)),
  Family = setNames(generate_colors(length(unique(annotation_row_rawp$Family)), "Paired"),
                    unique(annotation_row_rawp$Family)),
  Genus  = setNames(generate_colors(length(unique(annotation_row_rawp$Genus)), "Dark2"),
                    unique(annotation_row_rawp$Genus))
)

# Remover linhas com NA na matriz de correlação
na_rows <- rownames(cor_matrix_rawp)[apply(cor_matrix_rawp, 1, function(x) any(is.na(x)))]
cor_matrix_rawp_clean <- cor_matrix_rawp[!rownames(cor_matrix_rawp) %in% na_rows, ]
star_matrix_rawp_clean <- star_matrix_rawp[rownames(cor_matrix_rawp_clean), ]

# Atualizar anotação taxonômica
annotation_row_rawp_clean <- annotation_row_rawp[rownames(cor_matrix_rawp_clean), ]

# Heatmap final
pheatmap(cor_matrix_rawp_clean,
         display_numbers = star_matrix_rawp_clean,
         number_color = "white",
         clustering_method = "ward.D2",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = color_palette,
         breaks = breaks,
         annotation_row = annotation_row_rawp_clean,
         annotation_colors = ann_colors_rawp,
         main = "Correlação entre Taxons e Variáveis da Dieta \n(p-valor < 0.05, sem FDR)",
         fontsize = 16,
         fontsize_row = 7,
         fontsize_col = 9,
         fontsize_number = 12,
         border_color = NA)


#######



#quero salvar uma planilha com os pvalues:

library(dplyr)
library(tidyr)
library(readr)

# 1. Filtrar ASVs com pelo menos uma correlação significativa com variáveis de dieta (FDR < 0.05)
# Mantém todas as colunas importantes
asvs_signif_diet_df <- cor_all_diet %>%
  filter(p_adj < 0.05) %>%  # Apenas correlações significativas
  select(Feature.ID, variable, correlation, p_value, p_adj, star)  # Seleciona colunas de interesse

# 2. Juntar com a taxonomia baseada na coluna Feature.ID
asvs_signif_diet_df <- asvs_signif_diet_df %>%
  left_join(taxonomy, by = "Feature.ID")  # Adiciona coluna 'Taxon' e 'Confidence'

# 3. Separar a coluna Taxon em colunas taxonômicas distintas
asvs_signif_diet_df <- asvs_signif_diet_df %>%
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right", extra = "drop") %>%
  mutate(across(Domain:Species, ~ gsub("^[a-z]__*", "", .)))  # Remove prefixos como "g__", "f__" etc.

# 4. (Opcional) Reorganizar a ordem das colunas
asvs_signif_diet_df <- asvs_signif_diet_df %>%
  relocate(Feature.ID, Domain, Phylum, Class, Order, Family, Genus, Species,
           variable, correlation, p_value, p_adj, star)

# 5. (Opcional) Salvar em CSV
write_csv(asvs_signif_diet_df, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/asvs_significant_correlations_diet.csv")

