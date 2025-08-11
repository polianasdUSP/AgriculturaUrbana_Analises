Análises ASVs com Metabolic Score


# Carregar pacotes necessários
library(tidyverse)


# Instalar microbiomeMarker via Bioconductor
#BiocManager::install("microbiomeMarker")

library(microbiomeMarker)

library(qiime2R)
library(phyloseq)
#install.packages("xfun")

library(tibble)

#metadata <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/metadados_raw_V01_local_nascimento2.csv", 
#sep = ";", 
#header = TRUE, 
#stringsAsFactors = FALSE, 
#fileEncoding = "UTF-8")

#metadata <-metadados.raw

# Manter apenas as primeiras 130 linhas
#metadata <- metadata[1:130, ]

# Remover a primeira coluna
metadata <- metadados.metabolic.score


# Remover quaisquer linhas completamente vazias (se necessário)
metadata <- metadata[rowSums(is.na(metadata)) != ncol(metadata), ]



# Exibir as primeiras linhas para verificar se os dados foram importados corretamente
head(metadata)

# Verifique os nomes das colunas para garantir que foram carregados corretamente
colnames(metadata)



SVs<-read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/table.qza")$data
taxonomy<-read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/taxonomy_silva.qza")$data

SVs <- SVs_filtered

# Calcular o número de voluntários em que cada ASV está presente
asv_presence <- rowSums(SVs > 0)

# Filtrar ASVs presentes em pelo menos 13 voluntários
SVs_filtered <- SVs[asv_presence >= 10, ]




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
#write.csv(SVs_long, "SVs_normalized_log_transformed.csv", row.names = FALSE)




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
pheatmap(SVs_matrix, 
         clustering_method = "ward.D2",  # Método de clusterização
         clustering_distance_rows = "euclidean",  # Distância para ASVs
         clustering_distance_cols = "euclidean",  # Distância para amostras
         scale = "row",  # Escala por linha para melhor visualização
         main = "Heatmap de Abundância Normalizada e Log-transformada das ASVs",
         color = color_palette)

library(pheatmap)

# Selecionar as 50 ASVs mais abundantes
top_ASVs <- rowSums(SVs_matrix) %>% sort(decreasing = TRUE) %>% head(50) %>% names()
SVs_matrix_top <- SVs_matrix[top_ASVs, ]

# Criar a paleta de cores com branco centralizado
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Gerar o heatmap com as top 50 ASVs
pheatmap(SVs_matrix_top, 
         clustering_method = "ward.D2",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         scale = "row",
         main = "Heatmap das 50 ASVs Mais Abundantes (Normalizado e Log-transformado)",
         color = color_palette)


library(dplyr)





# Remover as colunas indesejadas do metadata
metadata_filtered <- metadados.all %>%
  select(-c("Gravida", "Menopausa", "Raca", "Gravida", "Weight", "Height",  "Hip","WHR"                ))

str(metadata_filtered)


# Juntar ASVs_long com metadata filtrado
# Assumindo que a coluna de junção seja "Sample.id" em ASVs_long e "SampleID" em metadata_filtered
metadata_ASVs <- SVs_long %>%
  inner_join(metadados.metabolic.score, by = c("Sample.id" = "Sample.id"))

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


library(dplyr)
library(ggplot2)

# Filtrar apenas colunas numéricas de metadata e manter Sample.id
numeric_metadata <- metadados.metabolic.score %>%
  select(Sample.id, where(is.numeric))

# Preparar dados para o heatmap usando SVs_long e combinando com colunas numéricas de metadata
heatmap_data <- SVs_long %>%
  left_join(numeric_metadata, by = "Sample.id") %>%
  left_join(taxonomy, by = "Feature.ID") %>%
  mutate(Feature = paste(Feature.ID, Taxon)) %>%
  mutate(Feature = gsub("[kpcofgs]__", "", Feature))  # Remover prefixos na taxonomia

# Plotar o heatmap
ggplot(heatmap_data, aes(x = Sample.id, y = Feature, fill = NormAbundance)) +
  geom_tile() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_c(name = "NormAbundance") +
  labs(title = "Heatmap de Abundância Normalizada das ASVs", x = "Sample ID", y = "ASV (Feature.ID e Taxonomia)")

# Salvar o gráfico como PDF
#ggsave("heatmap.pdf", height = 4, width = 11, device = "pdf")

library(dplyr)
library(ggplot2)


# Preparar dados para o heatmap usando SVs_long e combinando com colunas numéricas de metadata
heatmap_data <- SVs_long %>%
  left_join(numeric_metadata, by = "Sample.id") %>%
  left_join(taxonomy, by = "Feature.ID") 

# Plotar o heatmap
ggplot(heatmap_data, aes(x = Sample.id, y = Feature.ID, fill = NormAbundance)) +
  geom_tile() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_c(name = "NormAbundance") +
  labs(title = "Heatmap de Abundância Normalizada das ASVs", x = "Sample ID", y = "Feature.ID")

# Salvar o gráfico como PDF
#ggsave("heatmap_larger.pdf", height = 10, width = 20, device = "pdf")




library(dplyr)
library(tidyr)
library(pheatmap)
library(tibble)

# 1. Selecionar apenas colunas relevantes de saúde
# Corrigido
health_data <- metadata_ASVs %>%
  select(Sample.id, IL17A, IFNGamma, TNF, IL10, IL6, IL4, IL2, Age,
         Systolic, Diastolic, HbA1c, COLESTEROL, LDL,
         HDL, VLDL, TRIGLICERIDES, TGO, TGP, GGT, GLICOSE, INSULINA, HOMA.IR, 
         PCR, IMC, WHR, Health_Status) %>%
  distinct(Sample.id, .keep_all = TRUE)


# 2. Juntar com abundância
combined_data <- SVs_long %>%
  left_join(health_data, by = "Sample.id")

#Transformar Health_Status em numeros
combined_data <- combined_data %>%
  mutate(Health_Status_num = case_when(
    Health_Status == "Low Metabolic Risk" ~ 0,
    Health_Status == "Intermediate Metabolic Risk" ~ 1,
    Health_Status == "High Metabolic Risk" ~ 2
  ))

cor_results <- combined_data %>%
  group_by(Feature.ID) %>%
  summarize(across(
    .cols = c(IL17A, IFNGamma, TNF, IL10, IL6, IL4, IL2, Age,
              Systolic, Diastolic, HbA1c, COLESTEROL, LDL,
              HDL, VLDL, TRIGLICERIDES, TGO, TGP, GGT, GLICOSE,
              INSULINA, HOMA.IR, PCR, IMC, WHR, Health_Status_num),
    .fns = ~ suppressWarnings(cor(NormAbundance, .x, method = "spearman", use = "complete.obs")),
    .names = "cor_{.col}"
  )) %>%
  ungroup()


# 4. Transformar em formato longo
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

pheatmap(cor_matrix, 
         clustering_method = "ward.D2",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Heatmap de Correlação entre ASVs e Parâmetros de Saúde",
         color = color_palette, 
         breaks = breaks,
         border_color = NA)

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
health_vars <- c("IL17A", "IFNGamma", "TNF", "IL10", "IL6", "IL4", "IL2", "Age",
                 "Systolic", "Diastolic", "HbA1c", "COLESTEROL", "LDL",
                 "HDL", "VLDL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "GLICOSE", "INSULINA", "HOMA.IR", 
                 "PCR", "IMC", "WHR", "Health_Status_num")

# Criar uma lista vazia para guardar os resultados
cor_list <- list()

# Loop pelas variáveis de saúde
for (var in health_vars) {
  
  temp <- combined_data %>%
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
cor_all <- bind_rows(cor_list)

# Ajustar p-valor com FDR
cor_all <- cor_all %>%
  group_by(variable) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  ungroup()

# Marcar significância
cor_all <- cor_all %>%
  mutate(star = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01 ~ "**",
    p_adj < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Visualizar resultados formatados
head(cor_all)

cor_signif <- cor_all %>%
  filter(p_adj < 0.05)

cor_all %>%
  filter(p_adj < 0.05) %>%
  arrange(desc(abs(correlation)))


library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)



# 1. Selecionar apenas ASVs com pelo menos uma correlação significativa (FDR < 0.05)
asvs_signif <- cor_all %>%
  filter(p_adj < 0.05) %>%
  pull(Feature.ID) %>%
  unique()

# 2. Criar matriz de correlação apenas com ASVs significativas
cor_matrix <- cor_all %>%
  filter(Feature.ID %in% asvs_signif) %>%
  select(Feature.ID, variable, correlation) %>%
  pivot_wider(names_from = variable, values_from = correlation) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 3. Criar matriz de asteriscos para significância (p_adj)
star_matrix <- cor_all %>%
  filter(Feature.ID %in% asvs_signif) %>%
  select(Feature.ID, variable, star) %>%
  pivot_wider(names_from = variable, values_from = star) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# 4. Paleta de cores e quebras
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-0.5, 0.5, length.out = 101)

# 5. Gerar o heatmap com significância destacada
pheatmap(cor_matrix,
         display_numbers = star_matrix,
         clustering_method = "ward.D2",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = color_palette,
         breaks = breaks,
         main = "Correlação entre ASVs e Parâmetros de Saúde\n(* FDR < 0.05)",
         fontsize_row = 6,
         fontsize_col = 9,
         fontsize_number = 10,
         border_color = NA)


# Filtrar apenas ASVs com significância
asvs_signif_df <- cor_all %>%
  filter(p_adj < 0.05) %>%
  distinct(Feature.ID)

# Juntar com a tabela de taxonomia
asvs_tax_signif <- asvs_signif_df %>%
  left_join(taxonomy, by = "Feature.ID")

# Visualizar
head(asvs_tax_signif)

library(tidyr)

asvs_tax_signif <- asvs_tax_signif %>%
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right", remove = FALSE, extra = "drop") %>%
  mutate(across(Domain:Species, ~ gsub("^[a-z]__*", "", .)))  # Remove prefixos tipo "g__"

# Garantir que os Feature.IDs do heatmap estão no objeto com taxonomia
tax_annot <- taxonomy %>%
  filter(Feature.ID %in% rownames(cor_matrix)) %>%
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right", extra = "drop", remove = FALSE) %>%
  mutate(across(Domain:Species, ~ gsub("^[a-z]__*", "", .))) %>%
  select(Feature.ID, Phylum, Family, Genus)  # escolha os níveis que quer mostrar

# Criar data frame com os Feature.IDs como rownames
annotation_row <- tax_annot %>%
  column_to_rownames("Feature.ID")  # Importante: rownames devem ser iguais ao cor_matrix

pheatmap(cor_matrix,
         display_numbers = star_matrix,
         clustering_method = "ward.D2",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = color_palette,
         breaks = breaks,
         main = "Correlação entre ASVs e Parâmetros de Saúde\n(* FDR < 0.05)",
         fontsize_row = 6,
         fontsize_col = 9,
         fontsize_number = 10,
         border_color = NA,
         annotation_row = annotation_row)

#====================================#
#Dieta
#====================================#



library(dplyr)
library(tidyr)

# 1. Vetor com os nomes das variáveis de dieta
diet_vars <- c("carboidrato_total_g", "proteina_g", "lipidios_g",  "fibra_alimentar_g", "colesterol_mg",                   "acidos_graxos_saturados_g", "acidos_graxos_monoinsaturados_g", "acidos_graxos_poliinsaturados_g", "calcio_mg",                    "ferro_mg", "sodio_mg", "magnesio_mg", "fosforo_mg", "potassio_mg", "manganes_mg", "zinco_mg", "cobre_mg", "selenio_mcg",                 "vitamina_A_RAE_mcg", "vitamina_D_mcg", "vitamina_E_mg", "tiamina_mg", "riboflavina_mg", "niacina_mg", "vitamina_C_mg", "equivalente_de_folato_mcg", "sal_de_adicao_g", "acucar_de_adicao_g", "BHEI_R_Score_Total", "Percentual_NOVA_group_1",     "Percentual_NOVA_group_2",  "Percentual_NOVA_group_3", "Health_Status_num"    
)

#Transformar Health_Status em numeros
#metadados.metabolic.score_diet <- metadados.metabolic.score_diet %>%
 # mutate(Health_Status_num = case_when(
  #  Health_Status == "Low Metabolic Risk" ~ 0,
   # Health_Status == "Intermediate Metabolic Risk" ~ 1,
    #Health_Status == "High Metabolic Risk" ~ 2
#  ))

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
  temp <- metadados.metabolic.score_diet %>%
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

# 4. Anotação taxonômica
annotation_row_diet <- taxonomy %>%
  filter(Feature.ID %in% rownames(cor_matrix_diet)) %>%
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right", extra = "drop") %>%
  mutate(across(Domain:Species, ~ gsub("^[a-z]__*", "", .))) %>%
  select(Feature.ID, Phylum, Family, Genus) %>%
  column_to_rownames("Feature.ID")

# 5. Paleta de cores para anotação
library(RColorBrewer)

ann_colors_diet <- list(
  Phylum = setNames(brewer.pal(length(unique(annotation_row_diet$Phylum)), "Set3"),
                    unique(annotation_row_diet$Phylum)),
  Family = setNames(brewer.pal(length(unique(annotation_row_diet$Family)), "Paired"),
                    unique(annotation_row_diet$Family)),
  Genus = setNames(brewer.pal(length(unique(annotation_row_diet$Genus)), "Dark2"),
                   unique(annotation_row_diet$Genus))
)

# 6. Paleta e breaks para correlação
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-0.5, 0.5, length.out = 101)

# 7. Gerar heatmap
pheatmap(cor_matrix_diet,
         display_numbers = star_matrix_diet,
         annotation_row = annotation_row_diet,
         #annotation_colors = ann_colors_diet,
         clustering_method = "ward.D2",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = color_palette,
         breaks = breaks,
         main = "Correlação entre ASVs e Variáveis da Dieta\n(* FDR < 0.05)",
         fontsize_row = 6,
         fontsize_col = 9,
         fontsize_number = 10,
         border_color = NA)

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
pheatmap(cor_matrix_health,
         display_numbers = star_matrix_health,
         clustering_method = "ward.D2",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-0.5, 0.5, length.out = 101),
         main = "Correlação entre Taxa e Parâmetros de Saúde\n(* FDR < 0.05)",
         fontsize_row = 7,
         fontsize_col = 9,
         fontsize_number = 10,
         border_color = NA)

# Substituir IDs por taxonomia
cor_all_diet_named <- cor_all_diet %>%
  left_join(taxonomy_named, by = "Feature.ID") %>%
  filter(!is.na(Taxon_Label))

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
pheatmap(cor_matrix_diet,
         display_numbers = star_matrix_diet,
         clustering_method = "ward.D2",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-0.5, 0.5, length.out = 101),
         main = "Correlação entre Taxa e Variáveis da Dieta\n(* FDR < 0.05)",
         fontsize_row = 7,
         fontsize_col = 9,
         fontsize_number = 10,
         border_color = NA)


#============================================#
#Core Microbiota
#============================================#

# Pacotes
library(dplyr)
library(ggplot2)

# 1. Importar metadados
#metadados.metabolic.score <- read.csv(
 # "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados.metabolic.score.csv",
  #header = TRUE,
  #sep = ",",
  #stringsAsFactors = FALSE
#)

# 2. Função para core microbiota
get_core_microbiota <- function(SV_matrix, group_vector, group_name, threshold = 0.5) {
  group_samples <- colnames(SV_matrix)[group_vector == group_name]
  SV_group <- SV_matrix[, group_samples]
  presença <- rowSums(SV_group > 0) / length(group_samples)
  core_asvs <- presença[presença >= threshold]
  return(core_asvs)
}

# 3. Extrair vetor de grupo
Metabolic_group <- metadados.metabolic.score$Health_Status

# 4. Rodar core microbiota (ajuste o nome do objeto SVs_filtered se necessário)
core_low         <- get_core_microbiota(SVs_filtered, Metabolic_group, group_name = "Low Metabolic Risk", threshold = 0.5)
core_intermediate <- get_core_microbiota(SVs_filtered, Metabolic_group, group_name = "Intermediate Metabolic Risk", threshold = 0.5)
core_high        <- get_core_microbiota(SVs_filtered, Metabolic_group, group_name = "High Metabolic Risk", threshold = 0.5)

# 5. Transpor a matriz e juntar com metadados
SVs_t <- as.data.frame(t(SVs_filtered))
SVs_t$Sample.id <- rownames(SVs_t)

SVs_metabolic <- merge(
  SVs_t,
  metadados.metabolic.score[, c("Sample.id", "Health_Status")],
  by = "Sample.id"
)

# 6. Remover NAs
SVs_metabolic <- SVs_metabolic %>% filter(!is.na(Health_Status))

# 7. Preparar dados para análise
grupo_metabolic <- SVs_metabolic$Health_Status
SVs_test <- SVs_metabolic[, -c(1, ncol(SVs_metabolic))]  # Remove Sample.id e Health_Status

# 8. Kruskal-Wallis para p-valor
pvals_kw <- apply(SVs_test, 2, function(x) {
  kruskal.test(x ~ grupo_metabolic)$p.value
})

# 9. Calcular log2 Fold Change (High vs Low)
log2FC_metabolic <- apply(SVs_test, 2, function(x) {
  log2(mean(x[grupo_metabolic == "High Metabolic Risk"] + 1) /
         mean(x[grupo_metabolic == "Low Metabolic Risk"] + 1))
})

# 10. Dataframe de resultados
df_metabolic <- data.frame(
  ASV = names(pvals_kw),
  pvalue = pvals_kw,
  log2FC = log2FC_metabolic,
  FDR = p.adjust(pvals_kw, method = "fdr")
)

# 11. Marcar significância por FDR
df_metabolic$significant <- ifelse(df_metabolic$FDR < 0.05 & abs(df_metabolic$log2FC) > 1, "Yes", "No")

# 12. Plot Volcano (FDR)
df_metabolic_filtrado <- df_metabolic %>%
  filter(!is.na(log2FC), !is.na(FDR), is.finite(log2FC), is.finite(FDR))

ggplot(df_metabolic_filtrado, aes(x = log2FC, y = -log10(FDR), color = significant)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  scale_color_manual(values = c("No" = "gray", "Yes" = "red")) +
  labs(
    title = "Volcano Plot - Health Status (FDR)",
    x = "log2 Fold Change (High vs Low)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()

# 13. Marcar significância por p-valor cru
df_metabolic$pval_signif <- ifelse(!is.na(df_metabolic$pvalue) & df_metabolic$pvalue < 0.05 & abs(df_metabolic$log2FC) > 1, "Yes", "No")

# 14. Plot Volcano (p-valor)
ggplot(df_metabolic, aes(x = log2FC, y = -log10(pvalue), color = pval_signif)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  scale_color_manual(values = c("No" = "gray", "Yes" = "red")) +
  labs(
    title = "Volcano Plot - Health Status (p-valor cru)",
    x = "log2 Fold Change (High vs Low)",
    y = "-log10(p-value)"
  ) +
  theme_minimal()


#Separar taxonomias em colunas

library(tidyr)

taxonomy_sep <- taxonomy %>%
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right", remove = FALSE)



#Juntar com df_metabolic

df_metabolic_tax <- df_metabolic %>%
  left_join(taxonomy_sep, by = c("ASV" = "Feature.ID"))


#Adicionar os nomes em Volcano Plot

df_metabolic_tax <- df_metabolic_tax %>%
  mutate(
    Label = ifelse(!is.na(Genus) & Genus != "", paste0(Genus, " (", Family, ")"), ASV)
  )



#filtrar apenas os significativos e adicionar rótulos ao gráfico

df_metabolic_sig <- df_metabolic_tax %>%
  filter(pval_signif == "Yes")

ggplot(df_metabolic_tax, aes(x = log2FC, y = -log10(pvalue), color = pval_signif)) +
  geom_point(size = 2) +
  geom_text(
    data = df_metabolic_tax %>% filter(pval_signif == "Yes"),
    aes(label = Label),
    color = "black",           # <- Aqui é o ajuste importante
    vjust = -1, size = 3
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  scale_color_manual(values = c("No" = "gray", "Yes" = "red")) +
  labs(
    title = "Volcano Plot - Health Status (Raw p-value)",
    x = "log2 Fold Change (High vs Low)",
    y = "-log10(p-value)",
    color = "Significativo"
  ) +
  theme_minimal()


df_metabolic_tax <- df_metabolic_tax %>%
  mutate(
    Label = case_when(
      !is.na(Genus) & Genus != "" ~ paste0("", Genus, " (", Family, ")"),
      !is.na(Family) & Family != "" ~ paste0("", Family, " (", Order, ")"),
      !is.na(Order) & Order != "" ~ paste0("", Order, " (", Class, ")"),
      !is.na(Class) & Class != "" ~ paste0("", Class, " (", Phylum, ")"),
      !is.na(Phylum) & Phylum != "" ~ paste0("", Phylum),
      TRUE ~ ASV
    )
  )



df_metabolic_sig <- df_metabolic_tax %>% filter(pval_signif == "Yes")

library(ggrepel)

ggplot(df_metabolic_tax, aes(x = log2FC, y = -log10(pvalue), color = pval_signif)) +
  geom_point(size = 2) +
  geom_text_repel(
    data = df_metabolic_tax %>% filter(pval_signif == "Yes"),
    aes(label = Label),
    color = "black",
    size = 3,
    max.overlaps = Inf,
    force = 1.2,           # força de repulsão (ajustável)
    min.segment.length = 0 # mantém as setas curtas
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  scale_color_manual(values = c("No" = "gray", "Yes" = "red")) +
  labs(
    title = "Volcano Plot - Health Status (Raw p-value)",
    x = "log2 Fold Change (High vs Low)",
    y = "-log10(p-value)",
    color = "Significativo"
  ) +
  theme_minimal()



# Salvar em PNG com 300 dpi
ggsave(
  filename = "volcano_health_status_pvalor.png",
  plot = last_plot(),      # ou use o nome do seu objeto ggplot, se tiver salvo
  path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/",  # ajuste o caminho se quiser salvar em outra pasta
  dpi = 300,
  width = 12,     # em polegadas
  height = 6,
  units = "in"
)













#Esse volcano plot com p-valor cru mostra que algumas ASVs estão significativamente diferentes entre os grupos de Health_Status, com base nos critérios:

# p-valor < 0.05

# |log2FC| > 1 (ou seja, pelo menos 2x de diferença na média da abundância relativa entre os grupos)

# Essas ASVs estão destacadas em vermelho no gráfico.


asvs_signif_metabolic <- df_metabolic %>%
  filter(pvalue < 0.05, abs(log2FC) > 1) %>%
  arrange(pvalue)

# Visualizar
View(asvs_signif_metabolic)

# Supondo que o objeto com a taxonomia se chame 'taxonomy'
asvs_signif_metabolic_tax <- left_join(asvs_signif_metabolic, taxonomy, by = c("ASV" = "Feature.ID"))

library(ggrepel)

# Filtrar significativas com p < 0.05 e |log2FC| > 1
df_metabolic_tax <- df_metabolic %>%
  filter(pvalue < 0.05, abs(log2FC) > 1) %>%
  left_join(taxonomy, by = c("ASV" = "Feature.ID")) %>%
  mutate(label = stringr::str_extract(Taxon, "g__[^;]*"))  # extrai o gênero

# Volcano plot com rótulos
ggplot(df_metabolic, aes(x = log2FC, y = -log10(pvalue))) +
  geom_point(aes(color = pval_signif)) +
  geom_text_repel(data = df_metabolic_tax, aes(label = label), size = 3, max.overlaps = 10) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  scale_color_manual(values = c("gray", "red")) +
  labs(title = "Volcano Plot - metabolic (ASVs) with Taxonomy",
       x = "log2 Fold Change", y = "-log10(p-valor)") +
  theme_minimal()


# Salvar em PNG com 300 dpi
ggsave(
  filename = "volcano_health_status_pvalor_genero.png",
  plot = last_plot(),      # ou use o nome do seu objeto ggplot, se tiver salvo
  path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/",  # ajuste o caminho se quiser salvar em outra pasta
  dpi = 300,
  width = 12,     # em polegadas
  height = 6,
  units = "in"
)


#=======> Histograma

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# 1. Calcular abundância relativa por amostra (linhas = amostras)
SVs_rel <- sweep(SVs_filtered, 1, rowSums(SVs_filtered), FUN = "/")

# 2. Converter para data frame e adicionar Sample.id
SVs_long <- as.data.frame(t(SVs_rel))
SVs_long$Sample.id <- rownames(SVs_long)

# 3. Transformar em formato longo
dados_long <- SVs_long %>%
  pivot_longer(-Sample.id, names_to = "ASV", values_to = "Abundancia")

# 4. Juntar com metadados
dados_long <- dados_long %>%
  left_join(metadados.metabolic.score_diet[, c("Sample.id", "Health_Status")], by = "Sample.id") %>%
  filter(!is.na(Health_Status))

# 5. Juntar com taxonomia e extrair Genus
dados_long <- dados_long %>%
  left_join(taxonomy[, c("Feature.ID", "Taxon")], by = c("ASV" = "Feature.ID")) %>%
  mutate(
    Genus = str_extract(Taxon, "g__[^;]+"),
    Genus = gsub("g__", "", Genus),
    Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus)
  )

# 6. Calcular abundância média por grupo e gênero
dados_genus <- dados_long %>%
  group_by(Health_Status, Genus) %>%
  summarise(Media_Abundancia = mean(Abundancia), .groups = "drop")

# 7. Selecionar os 30 gêneros mais abundantes
top_genera <- dados_genus %>%
  group_by(Genus) %>%
  summarise(Total = sum(Media_Abundancia), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = 30) %>%
  pull(Genus)

# Filtrar para esses 30 gêneros
dados_top30 <- dados_genus %>%
  filter(Genus %in% top_genera)

# Garantir ordem dos grupos
dados_top30$Health_Status <- factor(dados_top30$Health_Status,
                                    levels = c("Low Metabolic Risk", "Intermediate Metabolic Risk", "High Metabolic Risk"))

# 8. Gráfico final
ggplot(dados_top30, aes(x = reorder(Genus, -Media_Abundancia), y = Media_Abundancia, fill = Health_Status)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Gênero",
    y = "Abundância relativa média",
    title = "Abundância relativa dos 30 gêneros mais abundantes por grupo de Health Status"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = c(
    "Low Metabolic Risk" = "#66c2a5",
    "Intermediate Metabolic Risk" = "#fc8d62",
    "High Metabolic Risk" = "#8da0cb"
  ))



ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/grafico_abundancia_generos_HealthStatus.png",
       width = 12, height = 8, dpi = 300)

