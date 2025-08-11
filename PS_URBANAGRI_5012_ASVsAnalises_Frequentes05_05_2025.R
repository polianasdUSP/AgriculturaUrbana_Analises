#_______________________________#
#     title: ASVs Análises"
#_______________________________#
  
  
#=== Carregar pacotes necessários ===#

library(tidyverse)
library(microbiomeMarker)
library(qiime2R)
library(phyloseq)
#install.packages("xfun")

library(tibble)

#====== Carregar Objetos =====#

# Carregando os dados
# Carregando os dados
metadados.all <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/metadados.all.corrigido - metadados.all.corrigido.csv") # Atualize com o caminho correto do arquivo

metadados.all.filtrado <- metadados.all[metadados.all$Sample.id != "S40142.F00", ]



# Exibir as primeiras linhas para verificar se os dados foram importados corretamente
head(metadados.all.filtrado)

# Verifique os nomes das colunas para garantir que foram carregados corretamente
colnames(metadados.all.filtrado)



SVs<-read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Qiime/Analises_Corrigidas/table_corrigida.qza")$data
taxonomy<-read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/taxonomy_silva.qza")$data



# Filtrar colunas de SVs com base nos IDs presentes em metadados.all.filtrado
SVs_filtrado <- SVs[, colnames(SVs) %in% metadados.all.filtrado$Sample.id]

#===== Filtragem, normalização e log das ASVs ==========#

# Calcular o número de voluntários em que cada ASV está presente
asv_presence <- as.data.frame(rowSums(SVs_filtrado > 0))

# Filtrar ASVs presentes em pelo menos 13 voluntários
SVs_filtered <- SVs_filtrado[asv_presence >= 10, ]


# Normalizar para porcentagem
SVs_normalized <- apply(SVs_filtered, 2, function(x) x / sum(x) * 100)  # Converte para porcentagem


# Transformar em data frame para manipulação
SVs_long <- SVs_normalized %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key = "Sample.id", value = "Abundance")


# Aplicar transformação logarítmica

SVs_long <- SVs_long %>%
  mutate(NormAbundance = log10(Abundance + 0.01))  # Adiciona 0.01 para evitar log(0)

# Exportar para uma planilha
write.csv(SVs_long, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/SVs_normalized_log_transformed.csv", row.names = FALSE)


#=========== Heatmap =============#


library(ggplot2)
library(pheatmap)

# Converter de volta para formato largo para o heatmap
SVs_matrix <- SVs_long %>%
  select(Feature.ID, Sample.id, NormAbundance) %>%
  spread(key = Sample.id, value = NormAbundance, fill = 0) %>%
  column_to_rownames("Feature.ID")



# Definir a paleta de cores com branco no centro
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)

# Gerar o heatmap com o ponto zero centralizado na cor branca
heatmap_ASVs <- pheatmap(SVs_matrix, 
                clustering_method = "ward.D2",  # Método de clusterização
                clustering_distance_rows = "euclidean",  # Distância para ASVs
                clustering_distance_cols = "euclidean",  # Distância para amostras
                scale = "row",  # Escala por linha para melhor visualização
                main = "Heatmap of the ASVs (Normalized and Log-Transformed)",
                color = color_palette)


ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_ASVs.png", heatmap_ASVs, width = 24, height = 16, dpi = 300)



library(pheatmap)

# Selecionar as 50 ASVs mais abundantes
top_ASVs <- rowSums(SVs_matrix) %>% sort(decreasing = TRUE) %>% head(50) %>% names()
SVs_matrix_top <- SVs_matrix[top_ASVs, ]

# Criar a paleta de cores com branco centralizado
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Gerar o heatmap com as top 50 ASVs
heatmap_ASVs_50_most <- pheatmap(SVs_matrix_top, 
                        clustering_method = "ward.D2",
                        clustering_distance_rows = "euclidean",
                        clustering_distance_cols = "euclidean",
                        scale = "row",
                        main = "Heatmap of the 50 Most Abundant ASVs (Normalized and Log-Transformed)",
         color = color_palette)


ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_ASVs_50_most.png", heatmap_ASVs_50_most, width = 24, height = 16, dpi = 300)

library(dplyr)



# Remover as colunas indesejadas do metadata
metadata_filtered <- metadados.all.filtrado %>%
  select(-c("Gravida", "Menopausa", "Raca", "Gravida", "Weight", "Height","Hip","WHR" ,  "energia2_kcal","sample.plate", "SoFAAS_Score", "Sat_Fat_Score", "Sodium_Score", "grain_roots_tubers_score", "milk_dairy_score",                 "total_fruit_score" ,  "total_vegetable_score" ,  "vegetable_oils_nuts_fishoil_score",  "dark_green_orange_veg_legume_score", "whole_fruit_score", "whole_grain_score"  , "ConsumoGrupo_NOVA_group_1" ,         "ConsumoGrupo_NOVA_group_2" ,  "ConsumoGrupo_NOVA_group_3" ,  "ConsumoCategoria"  , "BMI", "VAI" , "QUICKI", "METS_IR",  "TyG_BMI" ,  "TyG_WC" , "ACE" , "VAI"                                                                ))

str(metadata_filtered)


#==== juntar Metadada + ASVs =====#

# Juntar ASVs_long com metadata filtrado
# Assumindo que a coluna de junção seja "Sample.id" em ASVs_long e "SampleID" em metadata_filtered
metadata_ASVs <- SVs_long %>%
  inner_join(metadata_filtered, by = c("Sample.id" = "Sample.id"))


library(dplyr)

# Juntar o objeto taxonomy ao metadata_ASVs usando a coluna Feature.ID
metadata_ASVs <- metadata_ASVs %>%
  left_join(taxonomy %>% select(Feature.ID, Taxon, Confidence), by = "Feature.ID")

# Reorganizar para que Taxon seja a segunda coluna e Confidence a terceira
metadata_ASVs <- metadata_ASVs %>%
  select(Feature.ID, Taxon, Confidence, everything())

# Visualizar o resultado
metadata_ASVs

# Reorganizar para que Sample.id seja a primeira coluna
metadata_ASVs <- metadata_ASVs %>%
  select(Sample.id, everything())

# Visualizar o resultado
glimpse(metadata_ASVs)

#===================================#
#        Heatmap  
#===================================#

library(dplyr)
library(ggplot2)



# Plotar o heatmap
heatmap_normalized_abundance <- ggplot(metadata_ASVs, aes(x = Sample.id, y = Feature.ID, fill = NormAbundance)) +
  geom_tile() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_c(name = "NormAbundance") +
  labs(title = "Heatmap of the ASVs (Normalized and Log-Transformed)", x = "Sample ID", y = "ASV (Feature.ID e Taxonomia)")



ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_normalized_abundance.png", heatmap_normalized_abundance, width = 24, height = 16, dpi = 300)


#=================================================#
#-------- Inflammatory Markers vs ASVs -----------#
#=================================================#

#Criar um metadata+ ASVs inflamattory, tirando NAs de IL17A

metadata_filtered_inf <- metadata_filtered %>% filter(!is.na(IL17A))

# Juntar ASVs_long com metadata filtrado
# Assumindo que a coluna de junção seja "Sample.id" em ASVs_long e "Sample.id" em metadata_filtered
metadata_ASVs_inf <- SVs_long %>%
  inner_join(metadata_filtered_inf, by = c("Sample.id" = "Sample.id"))


library(dplyr)

# Juntar o objeto taxonomy ao metadata_ASVs usando a coluna Feature.ID
metadata_ASVs_inf <- metadata_ASVs_inf %>%
  left_join(taxonomy %>% select(Feature.ID, Taxon, Confidence), by = "Feature.ID")

# Reorganizar para que Taxon seja a segunda coluna e Confidence a terceira
metadata_ASVs_inf <- metadata_ASVs_inf %>%
  select(Feature.ID, Taxon, Confidence, everything())

# Visualizar o resultado
metadata_ASVs_inf

# Reorganizar para que Sample.id seja a primeira coluna
metadata_ASVs_inf <- metadata_ASVs_inf %>%
  select(Sample.id, everything())

# Visualizar o resultado
glimpse(metadata_ASVs_inf)
####

library(dplyr)
library(tidyr)
library(pheatmap)
library(tibble)

# 1. Selecionar apenas colunas relevantes de saúde
# Corrigido
inflammatory_markers_ASVs <- metadata_ASVs_inf %>%
  select(Sample.id,  Feature.ID, Taxon, Confidence,  Abundance, NormAbundance, IL17A, IFNGamma, TNF, IL10, IL6, IL4, IL2, INSULINA,  
         PCR) 


# 2. Calcular correlações para cada ASV com as variáveis de saúde
#cor_results <- inflammatory_markers_ASVs %>%
# group_by(Feature.ID) %>%
#  summarize(across(
#    .cols = c(IL17A, IFNGamma, TNF, IL10, IL6, IL4, IL2, INSULINA,  
#              PCR),
#    .fns = ~ cor(NormAbundance, .x, use = "complete.obs"),
#    .names = "cor_{.col}"
#  )) %>%
#  ungroup()



safe_cor <- function(x, y) {
  if (sum(complete.cases(x, y)) >= 2) {
    tryCatch(cor(x, y, use = "complete.obs"), error = function(e) NA_real_)
  } else {
    NA_real_
  }
}

cor_results <- inflammatory_markers_ASVs %>%
  group_by(Feature.ID) %>%
  summarise(
    cor_IL17A     = safe_cor(NormAbundance, IL17A),
    cor_IFNGamma  = safe_cor(NormAbundance, IFNGamma),
    cor_TNF       = safe_cor(NormAbundance, TNF),
    cor_IL10      = safe_cor(NormAbundance, IL10),
    cor_IL6       = safe_cor(NormAbundance, IL6),
    cor_IL4       = safe_cor(NormAbundance, IL4),
    cor_IL2       = safe_cor(NormAbundance, IL2),
    cor_INSULINA  = safe_cor(NormAbundance, INSULINA),
    cor_PCR       = safe_cor(NormAbundance, PCR)
  ) %>%
  ungroup()


#4. Transformar em formato longo
cor_long <- cor_results %>%
 pivot_longer(cols = starts_with("cor_"), names_to = "Health_Param", values_to = "Correlation") %>%
 mutate(Health_Param = gsub("cor_", "", Health_Param))

# 5. Criar matriz
cor_matrix <- cor_long %>%
  pivot_wider(names_from = Health_Param, values_from = Correlation) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 6. Heatmap
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)
breaks <- seq(-0.5, 0.5, length.out = 51)

p <- pheatmap(cor_matrix, 
         clustering_method = "ward.D2",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Correlation Heatmap between ASVs and Inflammatory Markers",
         color = color_palette, 
         breaks = breaks,
         border_color = NA)

ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_ASVvsInflamatMarkers.png", p,  width = 24, height = 16, dpi = 300)

library(dplyr)
library(tidyr)

# Função para calcular correlação e p-valor
calc_cor_p <- function(x, y) {
  valid <- complete.cases(x, y)
  if (sum(valid) > 2) {
    res <- suppressWarnings(cor.test(x[valid], y[valid]))
    return(c(correlation = as.numeric(res$estimate), p_value = res$p.value))
  } else {
    return(c(correlation = NA_real_, p_value = NA_real_))
  }
}


# Lista de variáveis de saúde
Inflammatory_vars <- c("IL17A", "IFNGamma", "TNF", "IL10", "IL6", "IL4", "IL2", "INSULINA",
                 "PCR")

# Criar uma lista vazia para guardar os resultados
cor_list <- list()

# Loop pelas variáveis de saúde
for (var in Inflammatory_vars) {
  
  temp <- inflammatory_markers_ASVs %>%
    group_by(Feature.ID) %>%
    summarise(
      correlation = calc_cor_p(NormAbundance, .data[[var]])["correlation"],
      p_value = calc_cor_p(NormAbundance, .data[[var]])["p_value"]
    ) %>%
    ungroup() %>%
    mutate(variable = var) %>%
    mutate(
      correlation = as.numeric(correlation),
      p_value = as.numeric(p_value)
    )
  
  # Armazenar na lista
  cor_list[[var]] <- temp
}


# Unir todos os resultados
cor_inflammatory <- bind_rows(cor_list)

# Ajustar p-valor com FDR
cor_inflammatory <- cor_inflammatory %>%
  group_by(variable) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  ungroup()

# Marcar significância
cor_inflammatory <- cor_inflammatory %>%
  mutate(star = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01 ~ "**",
    p_adj < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Visualizar resultados formatados
head(cor_inflammatory)


cor_signif_inflammatory <- cor_inflammatory %>%
  filter(p_adj < 0.05)

cor_inflammatory %>%
  filter(p_adj < 0.05) %>%
  arrange(desc(abs(correlation)))





library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)



# 1. Selecionar apenas ASVs com pelo menos uma correlação significativa (FDR < 0.05)
asvs_signif <- cor_inflammatory %>%
  filter(p_adj < 0.05) %>%
  pull(Feature.ID) %>%
  unique()

# 2. Criar matriz de correlação apenas com ASVs significativas
cor_matrix <- cor_inflammatory %>%
  filter(Feature.ID %in% asvs_signif) %>%
  select(Feature.ID, variable, correlation) %>%
  pivot_wider(names_from = variable, values_from = correlation) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 3. Criar matriz de asteriscos para significância (p_adj)
star_matrix <- cor_inflammatory %>%
  filter(Feature.ID %in% asvs_signif) %>%
  select(Feature.ID, variable, star) %>%
  pivot_wider(names_from = variable, values_from = star) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 4. Paleta de cores e quebras
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-0.5, 0.5, length.out = 101)

# 5. Gerar o heatmap com significância destacada
p <- pheatmap(cor_matrix,
         display_numbers = star_matrix,
         number_color = "white",  # cor dos asteriscos
         clustering_method = "ward.D2",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = color_palette,
         breaks = breaks,
         main = "Significant Correlations: ASVs vs Inflammatory Markers (n = 76) \n(FDR < 0.05)",
         fontsize = 18,
         fontsize_row = 18,
         fontsize_col = 18,
         fontsize_number = 26, # tamanho dos asteriscos
         border_color = NA)






ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_Significant_correlation_ASVvsInflamatMarkers.png", p,  width = 35, height = 20, dpi = 300)



# Filtrar apenas ASVs com significância
asvs_signif_df <- cor_inflammatory %>%
  filter(p_adj < 0.05) %>%
  distinct(Feature.ID)

# Juntar com a tabela de taxonomia
asvs_tax_signif <- asvs_signif_df %>%
  left_join(taxonomy, by = "Feature.ID")

# Visualizar
head(asvs_tax_signif)


library(dplyr)
library(tidyr)
library(pheatmap)

# 1. Criar coluna de nome taxonômico compacto
taxonomy_named <- taxonomy %>%
  filter(Feature.ID %in% rownames(cor_matrix)) %>%
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right", extra = "drop") %>%
  mutate(across(Domain:Species, ~ gsub("^[a-z]__*", "", .))) %>%
  mutate(Taxon_Label = paste(Phylum, Family, Genus, sep = " | ")) %>%
  select(Feature.ID, Taxon_Label)

# 2. Substituir rownames da matriz de correlação
rownames(cor_matrix) <- taxonomy_named$Taxon_Label[match(rownames(cor_matrix), taxonomy_named$Feature.ID)]
rownames(star_matrix) <- rownames(cor_matrix)  # para manter igual

# 3. Gerar o heatmap padrão com nomes taxonômicos
p <- pheatmap(cor_matrix,
              display_numbers = star_matrix,
              number_color = "white",
              clustering_method = "ward.D2",
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              color = color_palette,
              breaks = breaks,
              main = "Correlation Heatmap: Taxons vs Inflammatory Markers (n = 76)\n(FDR < 0.05)",
              fontsize = 18,
              fontsize_row = 18,
              fontsize_col = 18,
              fontsize_number = 26,
              border_color = NA)




ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_Significant_correlation_TaxonvsInflamatMarkers.png", p,  width = 35, height = 20, dpi = 300)



#=================================================#
#-------- Diet and BHEI vs ASVs -----------#
#=================================================#

#Criar um metadata+ ASVs Dieta, tirando NAs de BHEI_R_Score_Total

metadata_filtered_diet <- metadata_filtered %>% filter(!is.na(BHEI_R_Score_Total))

# Juntar ASVs_long com metadata filtrado
# Assumindo que a coluna de junção seja "Sample.id" em ASVs_long e "Sample.id" em metadata_filtered
metadata_ASVs_diet <- SVs_long %>%
  inner_join(metadata_filtered_diet, by = c("Sample.id" = "Sample.id"))


library(dplyr)

# Juntar o objeto taxonomy ao metadata_ASVs usando a coluna Feature.ID
metadata_ASVs_diet <- metadata_ASVs_diet %>%
  left_join(taxonomy %>% select(Feature.ID, Taxon, Confidence), by = "Feature.ID")

# Reorganizar para que Taxon seja a segunda coluna e Confidence a terceira
metadata_ASVs_diet <- metadata_ASVs_diet %>%
  select(Feature.ID, Taxon, Confidence, everything())

# Visualizar o resultado
metadata_ASVs_diet

# Reorganizar para que Sample.id seja a primeira coluna
metadata_ASVs_diet <- metadata_ASVs_diet %>%
  select(Sample.id, everything())

# Visualizar o resultado
glimpse(metadata_ASVs_diet)
####

library(dplyr)
library(tidyr)
library(pheatmap)
library(tibble)

# 1. Selecionar apenas colunas relevantes de saúde
# Corrigido
diet_ASVs <- metadata_ASVs_diet %>%
  select(Sample.id,  Feature.ID, Taxon, Confidence,  Abundance, NormAbundance, carboidrato_total_g, proteina_g, lipidios_g, fibra_alimentar_g, colesterol_mg,     acidos_graxos_saturados_g, acidos_graxos_monoinsaturados_g, acidos_graxos_poliinsaturados_g, acidos_graxos_trans_g, calcio_mg, ferro_mg,   sodio_mg, magnesio_mg, fosforo_mg, potassio_mg, manganes_mg, zinco_mg, cobre_mg, selenio_mcg, vitamina_A_RAE_mcg, vitamina_D_mcg,        vitamina_E_mg, tiamina_mg, riboflavina_mg, niacina_mg, vitamina_B6_mg, vitamina_B12_mcg, vitamina_C_mg, equivalente_de_folato_mcg, sal_de_adicao_g, acucar_de_adicao_g, BHEI_R_Score_Total) 





# 2. Calcular correlações para cada ASV com as variáveis de saúde
#cor_results <- diet_ASVs %>%
# group_by(Feature.ID) %>%
#  summarize(across(
#    .cols = c(IL17A, IFNGamma, TNF, IL10, IL6, IL4, IL2, INSULINA,  
#              PCR),
#    .fns = ~ cor(NormAbundance, .x, use = "complete.obs"),
#    .names = "cor_{.col}"
#  )) %>%
#  ungroup()



safe_cor <- function(x, y) {
  if (sum(complete.cases(x, y)) >= 2) {
    tryCatch(cor(x, y, use = "complete.obs"), error = function(e) NA_real_)
  } else {
    NA_real_
  }
}

cor_results <- diet_ASVs %>%
  group_by(Feature.ID) %>%
  summarise(cor_carboidrato_total_g     = safe_cor(NormAbundance, carboidrato_total_g),
            cor_proteina_g     = safe_cor(NormAbundance, proteina_g),
            cor_lipidios_g  = safe_cor(NormAbundance, lipidios_g),
            cor_fibra_alimentar_g       = safe_cor(NormAbundance, fibra_alimentar_g),
            cor_colesterol_mg      = safe_cor(NormAbundance, colesterol_mg),
            cor_acidos_graxos_saturados_g       = safe_cor(NormAbundance, acidos_graxos_saturados_g),
            cor_acidos_graxos_monoinsaturados_g       = safe_cor(NormAbundance, acidos_graxos_monoinsaturados_g),
            cor_acidos_graxos_poliinsaturados_g       = safe_cor(NormAbundance, acidos_graxos_poliinsaturados_g),
            cor_acidos_graxos_trans_g  = safe_cor(NormAbundance, acidos_graxos_trans_g),
            cor_calcio_mg       = safe_cor(NormAbundance, calcio_mg),
            cor_ferro_mg     = safe_cor(NormAbundance, ferro_mg),
            cor_sodio_mg  = safe_cor(NormAbundance, sodio_mg),
            cor_magnesio_mg       = safe_cor(NormAbundance, magnesio_mg),
            cor_colesterol_mg      = safe_cor(NormAbundance, colesterol_mg),
            cor_fosforo_mg       = safe_cor(NormAbundance, fosforo_mg),
            cor_potassio_mg       = safe_cor(NormAbundance, potassio_mg),
            cor_manganes_mg       = safe_cor(NormAbundance, manganes_mg),
            cor_zinco_mg  = safe_cor(NormAbundance, zinco_mg),
            cor_cobre_mg       = safe_cor(NormAbundance, cobre_mg),
            cor_selenio_mcg  = safe_cor(NormAbundance, selenio_mcg),
            cor_vitamina_A_RAE_mcg       = safe_cor(NormAbundance, vitamina_A_RAE_mcg),
            cor_vitamina_D_mcg      = safe_cor(NormAbundance, vitamina_D_mcg),
            cor_vitamina_E_mg       = safe_cor(NormAbundance, vitamina_E_mg),
            cor_tiamina_mg       = safe_cor(NormAbundance, tiamina_mg),
            cor_riboflavina_mg       = safe_cor(NormAbundance, riboflavina_mg),
            cor_niacina_mg  = safe_cor(NormAbundance, niacina_mg),
            cor_vitamina_B6_mg       = safe_cor(NormAbundance, vitamina_B6_mg),
            cor_vitamina_B12_mcg  = safe_cor(NormAbundance, vitamina_B12_mcg),
            cor_vitamina_C_mg       = safe_cor(NormAbundance, vitamina_C_mg),
            cor_equivalente_de_folato_mcg      = safe_cor(NormAbundance, equivalente_de_folato_mcg),
            cor_sal_de_adicao_g       = safe_cor(NormAbundance, sal_de_adicao_g),
            cor_acucar_de_adicao_g       = safe_cor(NormAbundance, acucar_de_adicao_g),
            cor_BHEI_R_Score_Total       = safe_cor(NormAbundance, BHEI_R_Score_Total)
            
            
            
  ) %>%
  ungroup()




#4. Transformar em formato longo
cor_long <- cor_results %>%
  pivot_longer(cols = starts_with("cor_"), names_to = "Health_Param", values_to = "Correlation") %>%
  mutate(Health_Param = gsub("cor_", "", Health_Param))

# 5. Criar matriz
cor_matrix <- cor_long %>%
  pivot_wider(names_from = Health_Param, values_from = Correlation) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 6. Heatmap
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)
breaks <- seq(-0.5, 0.5, length.out = 51)

p <- pheatmap(cor_matrix, 
              clustering_method = "ward.D2",
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              main = "Correlation Heatmap between ASVs and Diet",
              color = color_palette, 
              breaks = breaks,
              border_color = NA)

ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_ASVvsDiet.png", p,  width = 24, height = 16, dpi = 300)

library(dplyr)
library(tidyr)

# Função para calcular correlação e p-valor
calc_cor_p <- function(x, y) {
  valid <- complete.cases(x, y)
  if (sum(valid) > 2) {
    res <- suppressWarnings(cor.test(x[valid], y[valid]))
    return(c(correlation = as.numeric(res$estimate), p_value = res$p.value))
  } else {
    return(c(correlation = NA_real_, p_value = NA_real_))
  }
}


# Lista de variáveis de saúde
diet_vars <- c("carboidrato_total_g", "proteina_g", "lipidios_g", "fibra_alimentar_g", "colesterol_mg", "acidos_graxos_saturados_g", "acidos_graxos_monoinsaturados_g", "acidos_graxos_poliinsaturados_g", "acidos_graxos_trans_g", "calcio_mg", "ferro_mg", "sodio_mg", "magnesio_mg", "fosforo_mg", "potassio_mg", "manganes_mg", "zinco_mg", "cobre_mg", "selenio_mcg", "vitamina_A_RAE_mcg", "vitamina_D_mcg", "vitamina_E_mg", "tiamina_mg", "riboflavina_mg", "niacina_mg", "vitamina_B6_mg", "vitamina_B12_mcg", "vitamina_C_mg", "equivalente_de_folato_mcg", "sal_de_adicao_g", "acucar_de_adicao_g", "BHEI_R_Score_Total")

# Criar uma lista vazia para guardar os resultados
cor_list <- list()

# Diagnóstico rápido
diet_vars[!sapply(diet_ASVs[diet_vars], is.numeric)]


# Loop pelas variáveis de saúde
for (var in diet_vars) {
  
  temp <- diet_ASVs %>%
    group_by(Feature.ID) %>%
    summarise(
      correlation = calc_cor_p(NormAbundance, .data[[var]])["correlation"],
      p_value = calc_cor_p(NormAbundance, .data[[var]])["p_value"]
    ) %>%
    ungroup() %>%
    mutate(variable = var) %>%
    mutate(
      correlation = as.numeric(correlation),
      p_value = as.numeric(p_value)
    )
  
  # Armazenar na lista
  cor_list[[var]] <- temp
}


# Unir todos os resultados
cor_diet <- bind_rows(cor_list)

# Ajustar p-valor com FDR
cor_diet <- cor_diet %>%
  group_by(variable) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  ungroup()

# Marcar significância
cor_diet <- cor_diet %>%
  mutate(star = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01 ~ "**",
    p_adj < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Visualizar resultados formatados
head(cor_diet)


cor_signif_diet <- cor_diet %>%
  filter(p_adj < 0.05)

cor_diet %>%
  filter(p_adj < 0.05) %>%
  arrange(desc(abs(correlation)))





library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)



# 1. Selecionar apenas ASVs com pelo menos uma correlação significativa (FDR < 0.05)
asvs_signif <- cor_diet %>%
  filter(p_adj < 0.05) %>%
  pull(Feature.ID) %>%
  unique()

# 2. Criar matriz de correlação apenas com ASVs significativas
cor_matrix <- cor_diet %>%
  filter(Feature.ID %in% asvs_signif) %>%
  select(Feature.ID, variable, correlation) %>%
  pivot_wider(names_from = variable, values_from = correlation) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 3. Criar matriz de asteriscos para significância (p_adj)
star_matrix <- cor_diet %>%
  filter(Feature.ID %in% asvs_signif) %>%
  select(Feature.ID, variable, star) %>%
  pivot_wider(names_from = variable, values_from = star) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 4. Paleta de cores e quebras
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-0.5, 0.5, length.out = 101)

# 5. Gerar o heatmap com significância destacada
p <- pheatmap(cor_matrix,
              display_numbers = star_matrix,
              number_color = "white",  # cor dos asteriscos
              clustering_method = "ward.D2",
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              color = color_palette,
              breaks = breaks,
              main = "Significant Correlations: ASVs vs Diet Markers (n = 105) \n(FDR < 0.05)",
              fontsize = 18,
              fontsize_row = 18,
              fontsize_col = 18,
              fontsize_number = 20, # tamanho dos asteriscos
              border_color = NA)






ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_Significant_correlation_ASVvsDiet.png", p,  width = 35, height = 20, dpi = 300)



# Filtrar apenas ASVs com significância
asvs_signif_df <- cor_diet %>%
  filter(p_adj < 0.05) %>%
  distinct(Feature.ID)

# Juntar com a tabela de taxonomia
asvs_tax_signif <- asvs_signif_df %>%
  left_join(taxonomy, by = "Feature.ID")

# Visualizar
head(asvs_tax_signif)


library(dplyr)
library(tidyr)
library(pheatmap)

# 1. Criar coluna de nome taxonômico compacto
taxonomy_named <- taxonomy %>%
  filter(Feature.ID %in% rownames(cor_matrix)) %>%
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right", extra = "drop") %>%
  mutate(across(Domain:Species, ~ gsub("^[a-z]__*", "", .))) %>%
  mutate(Taxon_Label = paste(Phylum, Family, Genus, sep = " | ")) %>%
  select(Feature.ID, Taxon_Label)

# 2. Substituir rownames da matriz de correlação
rownames(cor_matrix) <- taxonomy_named$Taxon_Label[match(rownames(cor_matrix), taxonomy_named$Feature.ID)]
rownames(star_matrix) <- rownames(cor_matrix)  # para manter igual

# 3. Gerar o heatmap padrão com nomes taxonômicos
p <- pheatmap(cor_matrix,
              display_numbers = star_matrix,
              number_color = "white",
              clustering_method = "ward.D2",
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              color = color_palette,
              breaks = breaks,
              main = "Correlation Heatmap: Taxons vs Diet (n = 105)\n(FDR < 0.05)",
              fontsize = 18,
              fontsize_row = 18,
              fontsize_col = 18,
              fontsize_number = 16,
              border_color = NA)




ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_Significant_correlation_TaxonvsDiet.png", p,  width = 35, height = 20, dpi = 300)


#=================================================#
#-------- NOVA Markers vs ASVs -----------#
#=================================================#

#Criar um metadata+ ASVs NOVA, tirando NAs de IL17A

metadata_filtered_NOVA <- metadata_filtered %>% filter(!is.na(Percentual_NOVA_group_1))

# Juntar ASVs_long com metadata filtrado
# Assumindo que a coluna de junção seja "Sample.id" em ASVs_long e "Sample.id" em metadata_filtered
metadata_ASVs_NOVA <- SVs_long %>%
  inner_join(metadata_filtered_NOVA, by = c("Sample.id" = "Sample.id"))


library(dplyr)

# Juntar o objeto taxonomy ao metadata_ASVs usando a coluna Feature.ID
metadata_ASVs_NOVA <- metadata_ASVs_NOVA %>%
  left_join(taxonomy %>% select(Feature.ID, Taxon, Confidence), by = "Feature.ID")

# Reorganizar para que Taxon seja a segunda coluna e Confidence a terceira
metadata_ASVs_NOVA <- metadata_ASVs_NOVA %>%
  select(Feature.ID, Taxon, Confidence, everything())

# Visualizar o resultado
metadata_ASVs_NOVA

# Reorganizar para que Sample.id seja a primeira coluna
metadata_ASVs_NOVA <- metadata_ASVs_NOVA %>%
  select(Sample.id, everything())

# Visualizar o resultado
glimpse(metadata_ASVs_NOVA)
####

library(dplyr)
library(tidyr)
library(pheatmap)
library(tibble)

# 1. Selecionar apenas colunas relevantes de saúde
# Corrigido
NOVA_ASVs <- metadata_ASVs_NOVA %>%
  select(Sample.id,  Feature.ID, Taxon, Confidence,  Abundance, NormAbundance, Percentual_NOVA_group_1, Percentual_NOVA_group_2,          Percentual_NOVA_group_3) 


# 2. Calcular correlações para cada ASV com as variáveis de saúde
#cor_results <- NOVA_ASVs %>%
# group_by(Feature.ID) %>%
#  summarize(across(
#    .cols = c(IL17A, IFNGamma, TNF, IL10, IL6, IL4, IL2, INSULINA,  
#              PCR),
#    .fns = ~ cor(NormAbundance, .x, use = "complete.obs"),
#    .names = "cor_{.col}"
#  )) %>%
#  ungroup()



safe_cor <- function(x, y) {
  if (sum(complete.cases(x, y)) >= 2) {
    tryCatch(cor(x, y, use = "complete.obs"), error = function(e) NA_real_)
  } else {
    NA_real_
  }
}

cor_results <- NOVA_ASVs %>%
  group_by(Feature.ID) %>%
  summarise(
    cor_Percentual_NOVA_group_1     = safe_cor(NormAbundance, Percentual_NOVA_group_1),
    cor_Percentual_NOVA_group_2  = safe_cor(NormAbundance, Percentual_NOVA_group_2),
    cor_Percentual_NOVA_group_3       = safe_cor(NormAbundance, Percentual_NOVA_group_3)
    
  ) %>%
  ungroup()


#4. Transformar em formato longo
cor_long <- cor_results %>%
  pivot_longer(cols = starts_with("cor_"), names_to = "Health_Param", values_to = "Correlation") %>%
  mutate(Health_Param = gsub("cor_", "", Health_Param))

# 5. Criar matriz
cor_matrix <- cor_long %>%
  pivot_wider(names_from = Health_Param, values_from = Correlation) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 6. Heatmap
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)
breaks <- seq(-0.5, 0.5, length.out = 51)

p <- pheatmap(cor_matrix, 
              clustering_method = "ward.D2",
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              main = "Correlation Heatmap between ASVs and NOVA Markers (n=96)",
              color = color_palette, 
              breaks = breaks,
              border_color = NA)

ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_ASVvsNOVA.png", p,  width = 24, height = 16, dpi = 300)

library(dplyr)
library(tidyr)

# Função para calcular correlação e p-valor
calc_cor_p <- function(x, y) {
  valid <- complete.cases(x, y)
  if (sum(valid) > 2) {
    res <- suppressWarnings(cor.test(x[valid], y[valid]))
    return(c(correlation = as.numeric(res$estimate), p_value = res$p.value))
  } else {
    return(c(correlation = NA_real_, p_value = NA_real_))
  }
}


# Lista de variáveis de saúde
NOVA_vars <- c("Percentual_NOVA_group_1", "Percentual_NOVA_group_2", "Percentual_NOVA_group_3")

# Criar uma lista vazia para guardar os resultados
cor_list <- list()

# Loop pelas variáveis de saúde
for (var in NOVA_vars) {
  
  temp <- NOVA_ASVs %>%
    group_by(Feature.ID) %>%
    summarise(
      correlation = calc_cor_p(NormAbundance, .data[[var]])["correlation"],
      p_value = calc_cor_p(NormAbundance, .data[[var]])["p_value"]
    ) %>%
    ungroup() %>%
    mutate(variable = var) %>%
    mutate(
      correlation = as.numeric(correlation),
      p_value = as.numeric(p_value)
    )
  
  # Armazenar na lista
  cor_list[[var]] <- temp
}


# Unir todos os resultados
cor_NOVA <- bind_rows(cor_list)


# 1. Marcar significância com base no p-valor cru (antes de FDR)
cor_NOVA <- cor_NOVA %>%
  mutate(star_raw = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  ))


# 2. Matriz de correlação (ASVs com pelo menos uma correlação significativa pelo p-valor cru)
asvs_signif_raw <- cor_NOVA %>%
  filter(p_value < 0.01) %>%
  pull(Feature.ID) %>%
  unique()

cor_matrix_raw <- cor_NOVA %>%
  filter(Feature.ID %in% asvs_signif_raw) %>%
  select(Feature.ID, variable, correlation) %>%
  pivot_wider(names_from = variable, values_from = correlation) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 3. Matriz de asteriscos com base no p-valor cru
star_matrix_raw <- cor_NOVA %>%
  filter(Feature.ID %in% asvs_signif_raw) %>%
  select(Feature.ID, variable, star_raw) %>%
  pivot_wider(names_from = variable, values_from = star_raw) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 4. Substituir rownames por nomes taxonômicos
taxonomy_named <- taxonomy %>%
  filter(Feature.ID %in% rownames(cor_matrix_raw)) %>%
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right", extra = "drop") %>%
  mutate(across(Domain:Species, ~ gsub("^[a-z]__*", "", .))) %>%
  mutate(Taxon_Label = paste(Phylum, Family, Genus, sep = " | ")) %>%
  select(Feature.ID, Taxon_Label)

# Aplicar os nomes aos rownames
rownames(cor_matrix_raw) <- taxonomy_named$Taxon_Label[match(rownames(cor_matrix_raw), taxonomy_named$Feature.ID)]
rownames(star_matrix_raw) <- rownames(cor_matrix_raw)

# 5. Criar o heatmap com p-valor cru
library(pheatmap)

n <- pheatmap(cor_matrix_raw,
              display_numbers = star_matrix_raw,
              number_color = "white",
              clustering_method = "ward.D2",
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              color = colorRampPalette(c("blue", "white", "red"))(100),
              breaks = seq(-0.5, 0.5, length.out = 101),
              main = "Correlation Heatmap: Taxons vs NOVA Index (n = 96)\n(* = p < 0.05 before FDR)",
              fontsize = 18,
              fontsize_row = 18,
              fontsize_col = 18,
              fontsize_number = 24,
              border_color = NA)


ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_Significant_correlation_ASVvsNOVA_semfdr.png", n,  width = 35, height = 20, dpi = 300)

# Ajustar p-valor com FDR
cor_NOVA <- cor_NOVA %>%
  group_by(variable) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  ungroup()

# Marcar significância
cor_NOVA <- cor_NOVA %>%
  mutate(star = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01 ~ "**",
    p_adj < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Visualizar resultados formatados
head(cor_NOVA)


cor_signif_NOVA <- cor_NOVA %>%
  filter(p_adj < 0.05)

cor_NOVA %>%
  filter(p_adj < 0.05) %>%
  arrange(desc(abs(correlation)))





library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)



# 1. Selecionar apenas ASVs com pelo menos uma correlação significativa (FDR < 0.05)
asvs_signif <- cor_NOVA %>%
  filter(p_adj < 0.05) %>%
  pull(Feature.ID) %>%
  unique()

# 2. Criar matriz de correlação apenas com ASVs significativas
cor_matrix <- cor_NOVA %>%
  filter(Feature.ID %in% asvs_signif) %>%
  select(Feature.ID, variable, correlation) %>%
  pivot_wider(names_from = variable, values_from = correlation) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 3. Criar matriz de asteriscos para significância (p_adj)
star_matrix <- cor_NOVA %>%
  filter(Feature.ID %in% asvs_signif) %>%
  select(Feature.ID, variable, star) %>%
  pivot_wider(names_from = variable, values_from = star) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 4. Paleta de cores e quebras
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-0.5, 0.5, length.out = 101)

# 5. Gerar o heatmap com significância destacada
p <- pheatmap(cor_matrix,
              display_numbers = star_matrix,
              number_color = "white",  # cor dos asteriscos
              clustering_method = "ward.D2",
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              color = color_palette,
              breaks = breaks,
              main = "Significant Correlations: ASVs vs NOVA Index (n = 76) \n(FDR < 0.05)",
              fontsize = 18,
              fontsize_row = 18,
              fontsize_col = 18,
              fontsize_number = 26, # tamanho dos asteriscos
              border_color = NA)






ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_Significant_correlation_ASVvsNOVA.png", p,  width = 35, height = 20, dpi = 300)



# Filtrar apenas ASVs com significância
asvs_signif_df <- cor_NOVA %>%
  filter(p_adj < 0.05) %>%
  distinct(Feature.ID)

# Juntar com a tabela de taxonomia
asvs_tax_signif <- asvs_signif_df %>%
  left_join(taxonomy, by = "Feature.ID")

# Visualizar
head(asvs_tax_signif)


library(dplyr)
library(tidyr)
library(pheatmap)

# 1. Criar coluna de nome taxonômico compacto
taxonomy_named <- taxonomy %>%
  filter(Feature.ID %in% rownames(cor_matrix)) %>%
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right", extra = "drop") %>%
  mutate(across(Domain:Species, ~ gsub("^[a-z]__*", "", .))) %>%
  mutate(Taxon_Label = paste(Phylum, Family, Genus, sep = " | ")) %>%
  select(Feature.ID, Taxon_Label)

# 2. Substituir rownames da matriz de correlação
rownames(cor_matrix) <- taxonomy_named$Taxon_Label[match(rownames(cor_matrix), taxonomy_named$Feature.ID)]
rownames(star_matrix) <- rownames(cor_matrix)  # para manter igual

# 3. Gerar o heatmap padrão com nomes taxonômicos
p <- pheatmap(cor_matrix,
              display_numbers = star_matrix,
              number_color = "white",
              clustering_method = "ward.D2",
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              color = color_palette,
              breaks = breaks,
              main = "Correlation Heatmap: Taxons vs NOVA INDEX (n = 76)\n(FDR < 0.05)",
              fontsize = 18,
              fontsize_row = 18,
              fontsize_col = 18,
              fontsize_number = 26,
              border_color = NA)




ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_Significant_correlation_TaxonvsNOVA.png", p,  width = 35, height = 20, dpi = 300)


#=================================================#
#-------- biochemical Markers vs ASVs -----------#
#=================================================#

#Criar um metadata+ ASVs inflamattory, tirando NAs de HbA1c

metadata_filtered_bio <- metadata_filtered_bio %>% filter(!is.na(Systolic))

# Juntar ASVs_long com metadata filtrado
# Assumindo que a coluna de junção seja "Sample.id" em ASVs_long e "Sample.id" em metadata_filtered
metadata_ASVs_bio <- SVs_long %>%
  inner_join(metadata_filtered_bio, by = c("Sample.id" = "Sample.id"))


library(dplyr)

# Juntar o objeto taxonomy ao metadata_ASVs usando a coluna Feature.ID
metadata_ASVs_bio <- metadata_ASVs_bio %>%
  left_join(taxonomy %>% select(Feature.ID, Taxon, Confidence), by = "Feature.ID")

# Reorganizar para que Taxon seja a segunda coluna e Confidence a terceira
metadata_ASVs_bio <- metadata_ASVs_bio %>%
  select(Feature.ID, Taxon, Confidence, everything())

# Visualizar o resultado
metadata_ASVs_bio

# Reorganizar para que Sample.id seja a primeira coluna
metadata_ASVs_bio <- metadata_ASVs_bio %>%
  select(Sample.id, everything())

# Visualizar o resultado
glimpse(metadata_ASVs_bio)
####

library(dplyr)
library(tidyr)
library(pheatmap)
library(tibble)

# 1. Selecionar apenas colunas relevantes de saúde
# Corrigido
biochemical_markers_ASVs <- metadata_ASVs_bio %>%
  select(Sample.id,  Feature.ID, Taxon, Confidence,  Abundance, NormAbundance, Systolic, Diastolic, UREIA, CREATININA, HbA1c, COLESTEROL, LDL, HDL, VLDL, TRIGLICERIDES, TGO, TGP, GGT, GLICOSE) 



# 2. Calcular correlações para cada ASV com as variáveis de saúde
#cor_results <- biochemical_markers_ASVs %>%
# group_by(Feature.ID) %>%
#  summarize(across(
#    .cols = c(IL17A, IFNGamma, TNF, IL10, IL6, IL4, IL2, INSULINA,  
#              PCR),
#    .fns = ~ cor(NormAbundance, .x, use = "complete.obs"),
#    .names = "cor_{.col}"
#  )) %>%
#  ungroup()



safe_cor <- function(x, y) {
  if (sum(complete.cases(x, y)) >= 2) {
    tryCatch(cor(x, y, use = "complete.obs"), error = function(e) NA_real_)
  } else {
    NA_real_
  }
}

cor_results <- biochemical_markers_ASVs %>%
  group_by(Feature.ID) %>%
  summarise(
    cor_Systolic     = safe_cor(NormAbundance, Systolic),
    cor_Diastolic  = safe_cor(NormAbundance, Diastolic),
    cor_UREIA       = safe_cor(NormAbundance, UREIA),
    cor_CREATININA      = safe_cor(NormAbundance, CREATININA),
    cor_HbA1c       = safe_cor(NormAbundance, HbA1c),
    cor_COLESTEROL       = safe_cor(NormAbundance, COLESTEROL),
    cor_LDL       = safe_cor(NormAbundance, LDL),
    cor_HDL  = safe_cor(NormAbundance, HDL),
    cor_VLDL       = safe_cor(NormAbundance, VLDL),
    cor_TRIGLICERIDES       = safe_cor(NormAbundance, TRIGLICERIDES),
    cor_TGO       = safe_cor(NormAbundance, TGO),
    cor_TGP       = safe_cor(NormAbundance, TGP),
    cor_GGT       = safe_cor(NormAbundance, GGT),
    cor_GLICOSE       = safe_cor(NormAbundance, GLICOSE)
  ) %>%
  ungroup()




#4. Transformar em formato longo
cor_long <- cor_results %>%
  pivot_longer(cols = starts_with("cor_"), names_to = "Health_Param", values_to = "Correlation") %>%
  mutate(Health_Param = gsub("cor_", "", Health_Param))

# 5. Criar matriz
cor_matrix <- cor_long %>%
  pivot_wider(names_from = Health_Param, values_from = Correlation) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 6. Heatmap
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)
breaks <- seq(-0.5, 0.5, length.out = 51)

p <- pheatmap(cor_matrix, 
              clustering_method = "ward.D2",
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              main = "Correlation Heatmap between ASVs and Biochemical Markers",
              color = color_palette, 
              breaks = breaks,
              border_color = NA)

ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_ASVvsBiochemicalMarkers.png", p,  width = 24, height = 16, dpi = 300)

library(dplyr)
library(tidyr)

# Função para calcular correlação e p-valor
calc_cor_p <- function(x, y) {
  valid <- complete.cases(x, y)
  if (sum(valid) > 2) {
    res <- suppressWarnings(cor.test(x[valid], y[valid]))
    return(c(correlation = as.numeric(res$estimate), p_value = res$p.value))
  } else {
    return(c(correlation = NA_real_, p_value = NA_real_))
  }
}


# Lista de variáveis de saúde
biochemical_vars <- c("Systolic", "Diastolic", "UREIA", "CREATININA", "HbA1c", "COLESTEROL", "LDL", "HDL", "VLDL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "GLICOSE")

# Criar uma lista vazia para guardar os resultados
cor_list <- list()

# Loop pelas variáveis de saúde
for (var in biochemical_vars) {
  
  temp <- biochemical_markers_ASVs %>%
    group_by(Feature.ID) %>%
    summarise(
      correlation = calc_cor_p(NormAbundance, .data[[var]])["correlation"],
      p_value = calc_cor_p(NormAbundance, .data[[var]])["p_value"]
    ) %>%
    ungroup() %>%
    mutate(variable = var) %>%
    mutate(
      correlation = as.numeric(correlation),
      p_value = as.numeric(p_value)
    )
  
  # Armazenar na lista
  cor_list[[var]] <- temp
}


# Unir todos os resultados
cor_biochemical <- bind_rows(cor_list)

# Ajustar p-valor com FDR
cor_biochemical <- cor_biochemical %>%
  group_by(variable) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  ungroup()

# Marcar significância
cor_biochemical <- cor_biochemical %>%
  mutate(star = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01 ~ "**",
    p_adj < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Visualizar resultados formatados
head(cor_biochemical)


cor_signif_biochemical <- cor_biochemical %>%
  filter(p_adj < 0.05)

cor_biochemical %>%
  filter(p_adj < 0.05) %>%
  arrange(desc(abs(correlation)))





library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)



# 1. Selecionar apenas ASVs com pelo menos uma correlação significativa (FDR < 0.05)
asvs_signif <- cor_biochemical %>%
  filter(p_adj < 0.05) %>%
  pull(Feature.ID) %>%
  unique()

# 2. Criar matriz de correlação apenas com ASVs significativas
cor_matrix <- cor_biochemical %>%
  filter(Feature.ID %in% asvs_signif) %>%
  select(Feature.ID, variable, correlation) %>%
  pivot_wider(names_from = variable, values_from = correlation) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 3. Criar matriz de asteriscos para significância (p_adj)
star_matrix <- cor_biochemical %>%
  filter(Feature.ID %in% asvs_signif) %>%
  select(Feature.ID, variable, star) %>%
  pivot_wider(names_from = variable, values_from = star) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 4. Paleta de cores e quebras
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-0.5, 0.5, length.out = 101)

# 5. Gerar o heatmap com significância destacada
p <- pheatmap(cor_matrix,
              display_numbers = star_matrix,
              number_color = "white",  # cor dos asteriscos
              clustering_method = "ward.D2",
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              color = color_palette,
              breaks = breaks,
              main = "Significant Correlations: ASVs vs Biochemical Markers (n = 91) \n(FDR < 0.05)",
              fontsize = 18,
              fontsize_row = 18,
              fontsize_col = 18,
              fontsize_number = 26, # tamanho dos asteriscos
              border_color = NA)






ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_Significant_correlation_ASVvsBiochemicalMarkers.png", p,  width = 35, height = 20, dpi = 300)



# Filtrar apenas ASVs com significância
asvs_signif_df <- cor_biochemical %>%
  filter(p_adj < 0.05) %>%
  distinct(Feature.ID)

# Juntar com a tabela de taxonomia
asvs_tax_signif <- asvs_signif_df %>%
  left_join(taxonomy, by = "Feature.ID")

# Visualizar
head(asvs_tax_signif)


library(dplyr)
library(tidyr)
library(pheatmap)

# 1. Criar coluna de nome taxonômico compacto
taxonomy_named <- taxonomy %>%
  filter(Feature.ID %in% rownames(cor_matrix)) %>%
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right", extra = "drop") %>%
  mutate(across(Domain:Species, ~ gsub("^[a-z]__*", "", .))) %>%
  mutate(Taxon_Label = paste(Phylum, Family, Genus, sep = " | ")) %>%
  select(Feature.ID, Taxon_Label)

# 2. Substituir rownames da matriz de correlação
rownames(cor_matrix) <- taxonomy_named$Taxon_Label[match(rownames(cor_matrix), taxonomy_named$Feature.ID)]
rownames(star_matrix) <- rownames(cor_matrix)  # para manter igual

# 3. Gerar o heatmap padrão com nomes taxonômicos
p <- pheatmap(cor_matrix,
              display_numbers = star_matrix,
              number_color = "white",
              clustering_method = "ward.D2",
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              color = color_palette,
              breaks = breaks,
              main = "Correlation Heatmap: Taxons vs Biochemical Markers (n = 91)\n(FDR < 0.05)",
              fontsize = 18,
              fontsize_row = 18,
              fontsize_col = 18,
              fontsize_number = 26,
              border_color = NA)




ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_Significant_correlation_TaxonvsBiochemicalMarkers.png", p,  width = 35, height = 20, dpi = 300)




#============================================#
#           Core Microbiota
#============================================#


#Encontrar quais ASVs (linhas da matriz) estão presentes em pelo menos X% das amostras de um determinado grupo — esse conjunto é chamado de "core microbiota" para aquele grupo.
# Função para calcular core microbiota por grupo

get_core_microbiota <- function(SV_matrix, group_vector, group_name, threshold = 0.5) {
  
  #threshold de 0.5 = 50% das amostras
  
  # Subset da matriz de ASVs por grupo
  group_samples <- colnames(SV_matrix)[group_vector == group_name]
  SV_group <- SV_matrix[, group_samples]
  
  # Calcular a presença relativa (em % das amostras)
  presença <- rowSums(SV_group > 0) / length(group_samples)
  
  # Selecionar ASVs com presença acima do limiar
  core_asvs <- presença[presença >= threshold]
  
  return(core_asvs)
}



#===========================================#
#     TyG e Core Microbiota
#==========================================#


# Carregar pacotes necessários
library(tidyverse)
library(ggplot2)
library(FSA)
library(dplyr)

# Criar grupo binário com base na mediana de TyG
metadados.bioquimicos.alpha <- metadados.bioquimicos.alpha %>%
  mutate(TyG_group = ifelse(TyG <= median(TyG, na.rm = TRUE), "Low", "High"))

#Criar a coluna TyG_clinical com 3 categorias
metadados.bioquimicos.alpha <- metadados.bioquimicos.alpha %>%
  mutate(TyG_clinical = case_when(
    TyG < 8.0 ~ "Low risk",
    TyG >= 8.0 & TyG < 8.5 ~ "Intermediate risk",
    TyG >= 8.5 ~ "High risk",
    TRUE ~ NA_character_
  ))

#Verificar a distribuição:
table(metadados.bioquimicos.alpha$TyG_clinical, useNA = "ifany")

#Novo agrupamento: "At Risk" vs "Low risk"
metadados.bioquimicos.alpha <- metadados.bioquimicos.alpha %>%
  mutate(TyG_risk_group = case_when(
    TyG_clinical == "Low risk" ~ "Low",
    TyG_clinical %in% c("Intermediate risk", "High risk") ~ "At risk",
    TRUE ~ NA_character_
  ))

#E transformar em fator com ordem:
metadados.bioquimicos.alpha$TyG_risk_group <- factor(
  metadados.bioquimicos.alpha$TyG_risk_group,
  levels = c("Low", "At risk")
)


#Verificar a nova distribuição:
table(metadados.bioquimicos.alpha$TyG_risk_group, useNA = "ifany")



# 1. Filtrar apenas as amostras com TyG_group definido
dados_tyg <- metadados.bioquimicos.alpha %>% 
  filter(!is.na(TyG_group))

# 2. Descobrir quais IDs estão duplicados
#duplicados <- dados_tyg$Sample.id[duplicated(dados_tyg$Sample.id)]
#print(unique(duplicados))
# [1] "S40061.F00" "S40141.F00"

# Visualizar as duplicatas para inspecionar (opcional)
#dados_tyg %>%
#  filter(Sample.id %in% c("S40061.F00", "S40141.F00")) %>%
#  arrange(Sample.id)

# Remover duplicatas mantendo a primeira ocorrência
#dados_tyg <- dados_tyg %>% distinct(Sample.id, .keep_all = TRUE)

# 3. Alinhar os metadados com as colunas da matriz de SVs
# Filtrar para manter só amostras presentes na matriz
dados_tyg <- dados_tyg %>% 
  filter(Sample.id %in% colnames(SVs_filtrado)) %>%
  arrange(match(Sample.id, colnames(SVs_filtrado)))

# Reordenar a matriz SVs de acordo com os metadados
SVs_tyg <- SVs_filtrado[, dados_tyg$Sample.id]

# 4. Calcular abundância relativa (se ainda não estiver feita)
SVs_tyg_rel_filtrado <- sweep(SVs_tyg, 2, colSums(SVs_tyg), FUN = "/")

# Criar vetor de grupos na ordem correta
grupo_tyg <- dados_tyg$TyG_group

# Verificar se o tamanho bate: deve ser 57 (mesmo número de colunas)
if(length(grupo_tyg) != ncol(SVs_tyg_rel_filtrado)) {
  stop("O comprimento de grupo_tyg não bate com o número de amostras na matriz SVs_tyg_rel_filtrado")
}

# 5. Separar as amostras dos grupos "Low" e "High"
SVs_tyg_grupo_low <- SVs_tyg_rel_filtrado[, grupo_tyg == "Low"]
SVs_tyg_grupo_high <- SVs_tyg_rel_filtrado[, grupo_tyg == "High"]

# 6. Calcular médias para cada ASV e log2 Fold Change
media_low <- rowMeans(SVs_tyg_grupo_low)
media_high <- rowMeans(SVs_tyg_grupo_high)

# Adiciona um pseudocount para evitar divisão por zero e log de zero
log2fc <- log2((media_high + 1e-6) / (media_low + 1e-6))

# Calcular p-values para cada ASV com o teste de Wilcoxon
pvalues <- apply(SVs_tyg_rel_filtrado, 1, function(x) {
  grupo1 <- x[grupo_tyg == "Low"]
  grupo2 <- x[grupo_tyg == "High"]
  # Realiza o teste somente se houver variação em ambos os grupos
  if (length(unique(grupo1)) > 1 && length(unique(grupo2)) > 1) {
    wilcox.test(grupo1, grupo2)$p.value
  } else {
    NA  # Caso não haja variação, retorna NA
  }
})

# Criar o dataframe com os resultados
resultado_volcano_tyg <- data.frame(
  ASV = rownames(SVs_tyg_rel_filtrado),
  mean_low = media_low,
  mean_high = media_high,
  log2FoldChange = log2fc,
  pvalue = pvalues
)

# Classificar a significância com base em pvalue bruto
resultado_volcano_tyg <- resultado_volcano_tyg %>%
  mutate(Significativo = ifelse(!is.na(pvalue) & pvalue < 0.05 & abs(log2FoldChange) > 1, "Sim", "Não"))

# Remover linhas com NA em pvalue ou log2FoldChange para plotar
resultado_volcano_tyg_filtrado <- resultado_volcano_tyg %>%
  filter(!is.na(pvalue), !is.na(log2FoldChange))

# 7. Gerar o Volcano Plot sem correção FDR

  v <- ggplot(resultado_volcano_tyg_filtrado, aes(x = log2FoldChange, y = -log10(pvalue), color = Significativo)) +
  geom_point(alpha = 0.8, size = 6) +  # aumenta o tamanho dos pontos
  scale_color_manual(values = c("Sim" = "red", "Não" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot - Triglyceride-Glucose (TyG) Index Group (pvalue)",
    x = "log2(Fold Change) (High vs Low)",
    y = "-log10(p-valor)",
    color = "Significativo"
  ) +
  theme_minimal(base_size = 18) +  # aumenta base geral
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  )


ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/volcanoPlot_Tyg_pvalue.png", v,  width = 35, height = 20, dpi = 300)


# Supondo que sua tabela de taxonomia se chame 'taxonomy' com coluna 'Feature.ID' e 'Taxon'
resultado_volcano_tyg_tax <- resultado_volcano_tyg %>%
  left_join(taxonomy, by = c("ASV" = "Feature.ID"))

library(tidyr)

resultado_volcano_tyg_tax_sep <- resultado_volcano_tyg_tax %>%
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right", remove = FALSE)


resultado_volcano_tyg_tax_sep <- resultado_volcano_tyg_tax_sep %>%
  mutate(Significativo_pval = ifelse(!is.na(pvalue) & pvalue < 0.05 & abs(log2FoldChange) > 1, "Sim", "Não"))


resultado_volcano_tyg_tax_sep %>%
  filter(Significativo_pval == "Sim") %>%
  select(ASV, log2FoldChange, pvalue, Phylum, Family, Genus) %>%
  arrange(pvalue)


library(ggplot2)
library(ggrepel)  # para rótulos que não se sobrepõem


library(ggplot2)
library(ggrepel)
library(dplyr)

# 1. Corrigir prefixos e criar Label mais informativo
resultado_volcano_tyg_tax_sep <- resultado_volcano_tyg_tax_sep %>%
  mutate(
    Genus_clean  = gsub("^g__", "", Genus),
    Family_clean = gsub("^f__", "", Family),
    Order_clean  = gsub("^o__", "", Order),
    Label = case_when(
      !is.na(Genus_clean) & Genus_clean != "" ~ paste0("", Genus_clean),
      !is.na(Family_clean) & Family_clean != "" ~ paste0("", Family_clean),
      !is.na(Order_clean) & Order_clean != "" ~ paste0("", Order_clean),
      TRUE ~ "unclassified"
    )
  )

# 2. Volcano plot
v2 <- ggplot(resultado_volcano_tyg_tax_sep, aes(x = log2FoldChange, y = -log10(pvalue), color = Significativo_pval)) +
  geom_point(alpha = 0.85, size = 4) +
  scale_color_manual(values = c("Sim" = "red", "Não" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(
    data = subset(resultado_volcano_tyg_tax_sep, Significativo_pval == "Sim"),
    aes(label = Label),
    size = 6,
    color = "navy",
    max.overlaps = 50,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50"
  ) +
  labs(
    title = "Volcano Plot - Triglyceride-Glucose (TyG) Index Group",
    x = "log2(Fold Change) (High vs Low)",
    y = "-log10(p-valor)",
    color = "Significativo"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  )

# 3. Visualizar
print(v2)




ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/volcanoPlot_Tyg_taxon.png", v2,  width = 35, height = 20, dpi = 300)



#======> Volcano plot com TyG Clinical 

# Filtrar apenas os IDs com grupo TyG definido
amostras_validas <- metadados.bioquimicos.alpha %>%
  filter(!is.na(TyG_risk_group)) %>%
  select(Sample.id, TyG_risk_group)

# Filtrar a matriz de abundância para essas amostras
SVs_filtrado <- SVs[, amostras_validas$Sample.id]

# Garantir mesma ordem
grupo <- amostras_validas$TyG_risk_group
names(grupo) <- amostras_validas$Sample.id


resultados <- apply(SVs_filtrado, 1, function(x) {
  grupo1 <- x[grupo == "Low"]
  grupo2 <- x[grupo == "At risk"]
  
  # Evitar erros com poucos dados
  if (length(grupo1) > 2 & length(grupo2) > 2) {
    p <- tryCatch(wilcox.test(grupo1, grupo2)$p.value, error = function(e) NA)
    log2FC <- log2(mean(grupo2 + 1e-6) / mean(grupo1 + 1e-6))
    c(log2FC = log2FC, pvalue = p)
  } else {
    c(log2FC = NA, pvalue = NA)
  }
})

resultados_df <- as.data.frame(t(resultados))
resultados_df$padj <- p.adjust(resultados_df$pvalue, method = "fdr")
resultados_df$ASV <- rownames(resultados_df)

library(ggplot2)

ggplot(resultados_df, aes(x = log2FC, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05)) +
  scale_color_manual(values = c("gray", "blue")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    title = "Volcano Plot: TyG Risk (Low vs At Risk)",
    x = "log2 Fold Change",
    y = "-log10 Adjusted p-value"
  ) +
  theme_minimal()

#Identificar quais são esses ASVs
sig_asvs_tax <- sig_asvs %>%
  left_join(taxonomy, by = c("ASV" = "Feature.ID"))

print(sig_asvs)



library(tidyr)

#Separar a taxonomia em colunas
taxonomy_sep <- taxonomy %>%
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = "; ", fill = "right", remove = FALSE)

#Juntar com sig_asvs
sig_asvs_tax_sep <- sig_asvs %>%
  left_join(taxonomy_sep, by = c("ASV" = "Feature.ID"))

#Criar coluna de rótulo mais limpa
sig_asvs_tax_sep <- sig_asvs_tax_sep %>%
  mutate(Taxon_label = case_when(
    !is.na(Genus) & Genus != "g__uncultured" ~ Genus,
    !is.na(Family) & Family != "f__uncultured" ~ Family,
    !is.na(Order) ~ Order,
    TRUE ~ ASV
  ))



library(ggrepel)

ggplot(resultados_df, aes(x = log2FC, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05)) +
  geom_text_repel(
    data = sig_asvs_tax_sep,
    aes(x = log2FC, y = -log10(padj), label = Taxon_label),
    size = 3,
    max.overlaps = 10
  ) +
  scale_color_manual(values = c("gray", "blue")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    title = "Volcano Plot: TyG Risk (Low vs At Risk)",
    x = "log2 Fold Change",
    y = "-log10 Adjusted p-value"
  ) +
  theme_minimal()


#g__Clostridia_vadinBB60_group (mais abundante em At Risk)
#Esse grupo pertence à classe Clostridia, filo Firmicutes, e costuma aparecer em contextos:
  
#  Associação com dietas ricas em gordura saturada

#Ambientes de disbiose intestinal

#Já foi ligado a resistência à insulina e estresse metabólico

#Exemplo de artigo:
  
#  Li et al., 2020 (PMID: 32888540) – mostraram aumento de Clostridia_vadinBB60 em camundongos com dieta obesogênica.

#Também associado em estudos com pacientes com alterações glicêmicas e perfil lipídico alterado.

#📌 Interpretação: a associação com TyG alto faz sentido, pois esse índice é um marcador de resistência à insulina e disfunção metabólica.


#o__Rhodospirillales (também mais abundante em At Risk)
#Ordem dentro da classe Alphaproteobacteria, filo Proteobacteria

#O filo Proteobacteria como um todo é frequentemente descrito como marcador de disbiose

#Rhodospirillales pode estar envolvido em:
  
#  Resposta inflamatória

#Ambientes intestinais alterados

#Já foi descrito como aumentado em indivíduos com síndrome metabólica ou doença hepática gordurosa não alcoólica (NAFLD).

#📌 Interpretação: a presença aumentada em indivíduos com TyG elevado pode refletir um padrão inflamatório ou disbiótico, consistente com resistência à insulin


library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Filtrar ASVs significativas
asvs_significativas <- resultado_volcano_tyg_tax_sep %>%
  filter(Significativo_pval == "Sim") %>%
  pull(ASV)

# 2. Filtrar matriz de abundância relativa para essas ASVs
SVs_signif <- SVs_tyg_rel_filtrado[asvs_significativas, ]

# 3. Pegar taxonomia (Genus, Family, Order) das ASVs significativas
taxas_signif <- resultado_volcano_tyg_tax_sep %>%
  filter(ASV %in% asvs_significativas) %>%
  select(ASV, Genus, Family, Order, Class, Phylum, Domain)

# 4. Substituir nomes das linhas pela ASV (temporário, para manter referência)
rownames(SVs_signif) <- taxas_signif$ASV

# 5. Transpor matriz e converter para long format
SVs_long <- as.data.frame(t(SVs_signif))
SVs_long$Sample.id <- rownames(SVs_long)

SVs_long <- SVs_long %>%
  pivot_longer(-Sample.id, names_to = "ASV", values_to = "Abundancia")

grupo_tyg_df <- data.frame(
  Sample.id = colnames(SVs_tyg_rel_filtrado),
  TyG_group = grupo_tyg
)

# 6. Adicionar grupo TyG
SVs_long <- SVs_long %>%
  left_join(grupo_tyg_df, by = "Sample.id")  # grupo_tyg_df deve ter Sample.id e TyG_group

# 7. Juntar taxonomia
SVs_long <- SVs_long %>%
  left_join(taxas_signif, by = "ASV") %>%
  mutate(
    Genus_clean  = gsub("^g__", "", Genus),
    Family_clean = gsub("^f__", "", Family),
    Order_clean  = gsub("^o__", "", Order),
    Label = case_when(
      !is.na(Genus_clean) & Genus_clean != "" ~ paste0("", Genus_clean),
      !is.na(Family_clean) & Family_clean != "" ~ paste0("", Family_clean),
      !is.na(Order_clean) & Order_clean != "" ~ paste0("o_", Order_clean),
      TRUE ~ ASV
    )
  )

# 8. Calcular média por taxa (Label) e grupo TyG
medias <- SVs_long %>%
  group_by(Label, TyG_group) %>%
  summarise(Abundancia_media = mean(Abundancia), .groups = "drop")

# 9. Gráfico
h <- ggplot(medias, aes(x = reorder(Label, -Abundancia_media), y = Abundancia_media, fill = TyG_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(
    title = "Relative Abundance of Significant Taxa by Triglyceride-Glucose (TyG) Index Group",
    x = "Taxon",
    y = "Mean Relative Abundance"
  ) +
  scale_fill_manual(values = c("Low" = "#1f77b4", "High" = "#ff7f0e")) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(size = 26, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  )

# 10. Salvar
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/histograma_Tyg_taxon.png",
       h, width = 35, height = 20, dpi = 300)


#============================================================#
#         Tercis de Hb1Ac e Abundancia relativa    
#============================================================#

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# 1. Criar tercis de HbA1c
metadados_com_tercil <- metadados.bioquimicos.alpha %>%
  filter(Sample.id %in% colnames(SVs_filtrado)) %>%
  mutate(HbA1c_tercil = ntile(HbA1c, 3))

# 2. Calcular abundância relativa (caso ainda não tenha)
SVs_rel <- sweep(SVs_filtrado, 2, colSums(SVs_filtrado), FUN = "/") %>% as.data.frame()

# 3. Converter para formato longo
SVs_long <- SVs_rel %>%
  tibble::rownames_to_column("Feature.ID") %>%
  pivot_longer(-Feature.ID, names_to = "Sample.id", values_to = "Abundancia")

# 4. Juntar taxonomia e tercis
SVs_long <- SVs_long %>%
  left_join(taxonomy, by = "Feature.ID") %>%
  left_join(metadados_com_tercil[, c("Sample.id", "HbA1c_tercil")], by = "Sample.id") %>%
  mutate(Genus = stringr::str_extract(Taxon, "g__[^;]*"),
         Genus = gsub("g__", "", Genus)) %>%
  filter(!is.na(HbA1c_tercil), !is.na(Genus))

# 5. Calcular abundância média por tercil
abund_tercil_hba1c <- SVs_long %>%
  group_by(Genus, HbA1c_tercil) %>%
  summarise(Abundancia_Relativa = mean(Abundancia), .groups = "drop")

# 6. Gráfico de barras
ggplot(abund_tercil_hba1c, aes(x = fct_reorder(Genus, -Abundancia_Relativa), 
                               y = Abundancia_Relativa, fill = factor(HbA1c_tercil))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Relative Abundance of Genera by HbA1c Terciles",
       x = "Genus", y = "Mean Relative Abundance", fill = "HbA1c Tercile") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3"))


#====> Só as significativas


# 1. Selecionar ASVs significativas do volcano
asvs_significativas <- resultado_volcano_tyg_tax_sep %>%
  filter(Significativo_pval == "Sim") %>%  # ou use Significativo_fdr se tiver FDR
  pull(ASV)

# 2. Filtrar matriz SVs
SVs_sig <- SVs_filtrado[asvs_significativas, ]

# 3. Calcular abundância relativa
SVs_rel <- sweep(SVs_sig, 2, colSums(SVs_sig), FUN = "/") %>% as.data.frame()

# 4. Converter para formato longo
SVs_long <- SVs_rel %>%
  tibble::rownames_to_column("Feature.ID") %>%
  pivot_longer(-Feature.ID, names_to = "Sample.id", values_to = "Abundancia")

# 5. Juntar taxonomia e metadados com tercis de HbA1c
metadados_com_tercil <- metadados.bioquimicos.alpha %>%
  filter(Sample.id %in% colnames(SVs_filtrado)) %>%
  mutate(HbA1c_tercil = ntile(HbA1c, 3))

SVs_long <- SVs_long %>%
  left_join(taxonomy, by = "Feature.ID") %>%
  left_join(metadados_com_tercil[, c("Sample.id", "HbA1c_tercil")], by = "Sample.id") %>%
  mutate(
    Genus = stringr::str_extract(Taxon, "g__[^;]*"),
    Genus = gsub("g__", "", Genus)
  ) %>%
  filter(!is.na(HbA1c_tercil), !is.na(Genus))

# 6. Calcular abundância média por gênero e tercil
abund_tercil_sig <- SVs_long %>%
  group_by(Genus, HbA1c_tercil) %>%
  summarise(Abundancia_Relativa = mean(Abundancia), .groups = "drop")

h2 <- ggplot(abund_tercil_sig, aes(x = fct_reorder(Genus, -Abundancia_Relativa), 
                                   y = Abundancia_Relativa, fill = factor(HbA1c_tercil))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Relative Abundance of Significant Genera by HbA1c Terciles",
       x = "Genus", y = "Mean Relative Abundance", fill = "HbA1c Tercile") +
  scale_fill_manual(values = c("1" = "#1b9e77", "2" = "#d95f02", "3" = "#7570b3")) +
  theme_minimal(base_size = 18) +  # base geral
  theme(
    plot.title = element_text(size = 26, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text.x = element_text(size = 22, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  )


# 10. Salvar
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/histograma_HbA1ctercis_taxon.png",
       h2, width = 40, height = 20, dpi = 300)





#============> Pontos de corte clinicos


library(tidyverse)
library(forcats)

# 1. Classificar indivíduos conforme pontos de corte clínicos da HbA1c
metadados_com_grupo <- metadados.bioquimicos.alpha %>%
  filter(Sample.id %in% colnames(SVs_filtrado)) %>%
  mutate(HbA1c_cat = case_when(
    HbA1c < 5.7 ~ "Normal",
    HbA1c >= 5.7 & HbA1c < 6.5 ~ "Prediabetes",
    HbA1c >= 6.5 ~ "Diabetes",
    TRUE ~ NA_character_
  ))

# 2. Selecionar ASVs significativas
asvs_significativas <- resultado_volcano_tyg_tax_sep %>%
  filter(Significativo_pval == "Sim") %>%
  pull(ASV)

# 3. Filtrar matriz SVs para essas ASVs
SVs_sig <- SVs_filtrado[asvs_significativas, ]

# 4. Calcular abundância relativa
SVs_rel <- sweep(SVs_sig, 2, colSums(SVs_sig), FUN = "/") %>% as.data.frame()

# 5. Converter para formato longo
SVs_long <- SVs_rel %>%
  rownames_to_column("Feature.ID") %>%
  pivot_longer(-Feature.ID, names_to = "Sample.id", values_to = "Abundancia")

# 6. Juntar taxonomia e classificação da HbA1c
SVs_long <- SVs_long %>%
  left_join(taxonomy, by = "Feature.ID") %>%
  left_join(metadados_com_grupo[, c("Sample.id", "HbA1c_cat")], by = "Sample.id") %>%
  mutate(
    Genus = stringr::str_extract(Taxon, "g__[^;]*"),
    Genus = gsub("g__", "", Genus)
  ) %>%
  filter(!is.na(HbA1c_cat), !is.na(Genus))

# 7. Calcular média por grupo clínico
abund_grupo_clinico <- SVs_long %>%
  group_by(Genus, HbA1c_cat) %>%
  summarise(Abundancia_Relativa = mean(Abundancia, na.rm = TRUE), .groups = "drop")

# 8. Gráfico
h3 <- ggplot(abund_grupo_clinico, aes(x = fct_reorder(Genus, -Abundancia_Relativa), 
                                      y = Abundancia_Relativa, fill = HbA1c_cat)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Relative Abundance of Significant Genera by HbA1c Clinical Categories",
       x = "Genus", y = "Mean Relative Abundance", fill = "HbA1c Category") +
  theme_minimal(base_size = 16) +
  scale_fill_manual(values = c("Normal" = "#1b9e77", "Prediabetes" = "#d95f02", "Diabetes" = "#7570b3")) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(size = 22, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

# 9. Salvar (opcional)
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/histograma_HbA1c_clinico.png",
       h3, width = 35, height = 20, dpi = 300)



#==============Volcano Plot com Health Groups =========#

#1. Preparar os dados

# Garantir que SVs esteja com ASVs nas linhas e Sample.id nas colunas
# Garantir que metadados_clusterizado tem Sample.id e HealthGroup_named

# Filtrar metadados para as amostras presentes nos SVs
ids_comuns <- intersect(colnames(SVs), metadados_clusterizado$Sample.id)

SVs_filtrado <- SVs[, ids_comuns]
metadados_filtrado <- metadados_clusterizado %>%
  filter(Sample.id %in% ids_comuns)


#2. Função para gerar dados do volcano

library(dplyr)

gerar_volcano_data_sem_fdr <- function(grupo1, grupo2, metadados, matriz_SVs) {
  ids1 <- metadados %>% filter(HealthGroup_named == grupo1) %>% pull(Sample.id)
  ids2 <- metadados %>% filter(HealthGroup_named == grupo2) %>% pull(Sample.id)
  
  resultado <- apply(matriz_SVs, 1, function(x) {
    grupo1_vals <- as.numeric(x[ids1])
    grupo2_vals <- as.numeric(x[ids2])
    
    if (length(grupo1_vals) > 2 & length(grupo2_vals) > 2) {
      p <- tryCatch(wilcox.test(grupo1_vals, grupo2_vals, exact = FALSE)$p.value, error = function(e) NA)
      log2FC <- log2(mean(grupo1_vals + 1e-6) / mean(grupo2_vals + 1e-6))
      c(log2FC = log2FC, pvalue = p)
    } else {
      c(log2FC = NA, pvalue = NA)
    }
  })
  
  df <- as.data.frame(t(resultado))
  df$Feature.ID <- rownames(df)
  return(df)
}


#3. Gerar os dois datasets

volcano_risk_vs_healthy_raw <- gerar_volcano_data_sem_fdr("High Risk Metabolic Profile", "Optimal Health", metadados_filtrado, SVs_filtrado)

volcano_transition_vs_healthy <- gerar_volcano_data_sem_fdr("Optimal Health","Transition", metadados_filtrado, SVs_filtrado)


#4. Plotar os volcano plots

ggplot(volcano_risk_vs_healthy_raw, aes(x = log2FC, y = -log10(pvalue))) +
  geom_point(aes(color = pvalue < 0.05), alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("gray", "blue")) +
  labs(
    title = "Volcano Plot: High Risk vs Optimal Health (Raw p-values)",
    x = "log2 Fold Change",
    y = "-log10 p-value"
  ) +
  theme_minimal()


volcano_tax <- volcano_risk_vs_healthy_raw %>%
  left_join(taxonomy_sep, by = c("Feature.ID" = "Feature.ID"))

volcano_tax <- volcano_tax %>%
  mutate(taxon_label = case_when(
    !is.na(Genus) & Genus != "g__uncultured" ~ Genus,
    !is.na(Family) & Family != "f__uncultured" ~ Family,
    !is.na(Order) ~ Order,
    TRUE ~ Feature.ID
  ))




library(ggplot2)
library(ggrepel)



library(ggplot2)
library(ggrepel)

p2 <- ggplot(volcano_tax, aes(x = log2FC, y = -log10(pvalue))) +
  geom_point(aes(color = ifelse(pvalue < 0.05 & abs(log2FC) > 1.2, "Significant", "Not significant")), alpha = 0.7, size = 7) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = c(-1.2, 1.2), linetype = "dotted", color = "black", linewidth = 0.8) +
  
  geom_text_repel(
    data = subset(volcano_tax, pvalue < 0.05 & (log2FC < -1.2 | log2FC > 1.2)),
    aes(label = taxon_label),
    size = 12,              # AUMENTA o tamanho dos rótulos
    max.overlaps = 25,
    box.padding = 3,
    segment.size = 2
  ) +
  
  scale_color_manual(
    values = c("Not significant" = "gray", "Significant" = "blue"),
    name = "Significance"
  ) +
  
  labs(
    title = "Volcano Plot: High Risk vs Optimal Health",
    x = "log2 Fold Change",
    y = "-log10 p-value"
  ) +
  
  theme_minimal(base_size = 22) +        # AUMENTA TODA A BASE DE FONTE
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 28, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20)
  )


ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/volcano_risk_vs_healthy.png",
       p2, width = 35, height = 20, dpi = 300)


#====> volcano_transition_vs_healthy

# Juntar com a taxonomia
volcano_tax_transition <- volcano_transition_vs_healthy %>%
  left_join(taxonomy_sep, by = c("Feature.ID" = "Feature.ID"))

# Criar coluna com o melhor nível taxonômico disponível
volcano_tax_transition <- volcano_tax_transition %>%
  mutate(taxon_label = case_when(
    !is.na(Genus) & Genus != "g__uncultured" ~ Genus,
    !is.na(Family) & Family != "f__uncultured" ~ Family,
    !is.na(Order) ~ Order,
    TRUE ~ Feature.ID
  ))

# Criar o gráfico
p3 <- ggplot(volcano_tax_transition, aes(x = log2FC, y = -log10(pvalue))) +
  geom_point(aes(color = ifelse(pvalue < 0.05 & abs(log2FC) > 1.2, "Significant", "Not significant")),
             alpha = 0.7, size = 7) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = c(-1.2, 1.2), linetype = "dotted", color = "black", linewidth = 0.8) +
  
  geom_text_repel(
    data = subset(volcano_tax_transition, pvalue < 0.05 & (log2FC < -1.2 | log2FC > 1.2)),
    aes(label = taxon_label),
    size = 12,
    max.overlaps = 25,
    box.padding = 3,
    segment.size = 2
  ) +
  
  scale_color_manual(
    values = c("Not significant" = "gray", "Significant" = "blue"),
    name = "Significance"
  ) +
  
  labs(
    title = "Volcano Plot: Transition vs Optimal Health",
    x = "log2 Fold Change",
    y = "-log10 p-value"
  ) +
  
  theme_minimal(base_size = 22) +
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 28, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20)
  )

# Salvar o gráfico em alta resolução
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/volcano_transition_vs_healthy.png",
       p3, width = 35, height = 20, dpi = 300)

#=====>  Histograma

library(tidyverse)
library(forcats)

# 1. Juntar metadados com os grupos
metadados_com_grupo <- metadados_clusterizado %>%
  filter(Sample.id %in% colnames(SVs_filtrado)) %>%
  select(Sample.id, HealthGroup_named)

# 2. Selecionar ASVs significativas do volcano plot
asvs_significativas <- volcano_tax %>%
  filter(pvalue < 0.05 & abs(log2FC) > 1.2) %>%
  pull(Feature.ID)

# 3. Filtrar matriz SVs para essas ASVs
SVs_sig <- SVs_filtrado[asvs_significativas, ]

# 4. Calcular abundância relativa
SVs_rel <- sweep(SVs_sig, 2, colSums(SVs_sig), FUN = "/") %>% as.data.frame()

# 5. Converter para formato longo
SVs_long <- SVs_rel %>%
  rownames_to_column("Feature.ID") %>%
  pivot_longer(-Feature.ID, names_to = "Sample.id", values_to = "Abundancia")

# 6. Juntar taxonomia e classificação dos grupos
SVs_long <- SVs_long %>%
  left_join(taxonomy, by = "Feature.ID") %>%
  left_join(metadados_com_grupo, by = "Sample.id") %>%
  mutate(
    Genus = stringr::str_extract(Taxon, "g__[^;]*"),
    Genus = gsub("g__", "", Genus)
  ) %>%
  filter(!is.na(HealthGroup_named), !is.na(Genus))

# 7. Calcular média por grupo
abund_grupo <- SVs_long %>%
  group_by(Genus, HealthGroup_named) %>%
  summarise(Abundancia_Relativa = mean(Abundancia, na.rm = TRUE), .groups = "drop")


SVs_long <- SVs_long %>%
  mutate(
    Family = stringr::str_extract(Taxon, "f__[^;]*") %>% gsub("f__", "", .),
    Order = stringr::str_extract(Taxon, "o__[^;]*") %>% gsub("o__", "", .),
    Genero_label = case_when(
      !is.na(Genus) & Genus != "uncultured" & !grepl("^[A-Z0-9\\-_]{5,}$", Genus) ~ Genus,
      !is.na(Family) & Family != "uncultured" ~ paste0("f_", Family),
      !is.na(Order) ~ paste0("o_", Order),
      TRUE ~ "Unknown"
    )
  ) %>%
  filter(!is.na(HealthGroup_named), !is.na(Genero_label))


SVs_long <- SVs_rel %>%
  rownames_to_column("Feature.ID") %>%
  pivot_longer(-Feature.ID, names_to = "Sample.id", values_to = "Abundancia") %>%
  left_join(taxonomy %>% select(Feature.ID, Taxon), by = "Feature.ID") %>%
  left_join(metadados_com_grupo %>% select(Sample.id, HealthGroup_named), by = "Sample.id") %>%
  mutate(
    Genus = str_extract(Taxon, "g__[^;]*") %>% gsub("g__", "", .),
    Family = str_extract(Taxon, "f__[^;]*") %>% gsub("f__", "", .),
    Order = str_extract(Taxon, "o__[^;]*") %>% gsub("o__", "", .),
    Class = str_extract(Taxon, "c__[^;]*") %>% gsub("c__", "", .),
    Phylum = str_extract(Taxon, "p__[^;]*") %>% gsub("p__", "", .),
    
    Genero_label = case_when(
      !is.na(Genus) & Genus != "uncultured" & !grepl("^[A-Z0-9\\-_]{4,}$", Genus) ~ Genus,
      !is.na(Family) & Family != "uncultured" & !grepl("^[A-Z0-9\\-_]{4,}$", Family) ~ paste0("f_", Family),
      !is.na(Order) & Order != "uncultured" ~ paste0("o_", Order),
      !is.na(Class) & Class != "uncultured" ~ paste0("c_", Class),
      !is.na(Phylum) & Phylum != "uncultured" ~ paste0("p_", Phylum),
      TRUE ~ "Unknown"
    )
  ) %>%
  filter(!is.na(HealthGroup_named), !is.na(Genero_label))

#Ordem dos grupos no gráfico (Optimal → Transition → High Risk)
#Você pode ajustar isso com factor() + levels= no mutate():
SVs_long$HealthGroup_named <- factor(
  SVs_long$HealthGroup_named,
  levels = c("Optimal Health", "In Transition", "High Risk Metabolic Profile")
)


#entender pq treponema parece alto
SVs_long %>%
  filter(Genero_label == "Treponema") %>%
  group_by(HealthGroup_named) %>%
  summarise(mean_abund = mean(Abundancia), sd_abund = sd(Abundancia), .groups = "drop")

#quantas pessoas tem treponema?!
SVs_long %>%
  filter(Genero_label == "Treponema", Abundancia > 0) %>%
  group_by(HealthGroup_named) %>%
  summarise(n = n(), .groups = "drop")



library(ggplot2)
library(forcats)

# Calcular abundância média por grupo de saúde
# Ordenar pelo grupo mais abundante
abund_grupo_saude <- abund_grupo_saude %>%
  group_by(Genero_label) %>%
  mutate(max_abund = max(Abundancia_Relativa)) %>%
  ungroup() %>%
  arrange(desc(max_abund))


# Gráfico de barras
h4 <- ggplot(abund_grupo_saude, aes(x = fct_reorder(Genero_label, -Abundancia_Relativa), 
                                    y = Abundancia_Relativa, fill = HealthGroup_named)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Relative Abundance of Significant Taxa by Health Profile",
       x = "Genus / Family / Order",
       y = "Mean Relative Abundance",
       fill = "Health Group") +
  theme_minimal(base_size = 16) +
  scale_fill_manual(values = c(
    "Optimal Health" = "#1b9e77",
    "In Transition" = "#d95f02",
    "High Risk Metabolic Profile" = "#7570b3"
  )) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(size = 40, face = "bold"),
    axis.title = element_text(size = 36),
    axis.text = element_text(size = 26),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 20)
  )

# Exibir o gráfico
print(h4)



# 9. Salvar (opcional)
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/histograma_core_microbiome_clusters.png",
       h4, width = 35, height = 20, dpi = 300)


h5 <- ggplot(abund_grupo_saude, aes(x = fct_reorder(Genero_label, -Abundancia_Relativa), 
                                    y = Abundancia_Relativa, fill = HealthGroup_named)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Relative Abundance of Significant Taxa by Health Profile",
       x = "Genus / Family / Order",
       y = "Mean Relative Abundance",
       fill = "Health Group") +
  theme_minimal(base_size = 16) +
  scale_fill_manual(
    values = c(
      "Optimal Health" = "#1b9e77",
      "In Transition" = "#d95f02",
      "High Risk Metabolic Profile" = "#7570b3"
    )
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(size = 40, face = "bold"),
    axis.title = element_text(size = 36),
    axis.text = element_text(size = 26),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 20)
  )



