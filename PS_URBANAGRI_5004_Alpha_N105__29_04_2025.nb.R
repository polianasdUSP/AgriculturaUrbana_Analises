#======================================#
#   title: "Alpha Diversity Analisys"
#     N=105
#======================================#

#Por causa de contaminação das amostras, terei que calcular novamente as métricas de alpha e beta diversidade, returando as amostras contaminadas
# Novos ids removidos: "S10041.F00", "S10232.F00", "S30092.F00", "S30191.F00"
 
#=========== Carregar as bibliotecas necessárias=================#

library(ggplot2)
library(reshape2)
library(dplyr)
library(ggcorrplot)
library(pheatmap)
library(grid)
library(patchwork)
library(ggpubr)
library(ggpmisc)
library(phyloseq)
library(microbiomeMarker)
library(tidyverse)
library(vegan)
library(DESeq2)
library(data.table)
library(qiime2R)
library(viridis)
library(ape)




#============ Importar o arquivo CSV ==============================#



metadados.all.filtrado <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/AgriculturaUrbana_Analises/metadados_all_filtrado.csv", header = TRUE, sep = ",")


#metadata <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Qiime/Analises_Corrigidas/metadados_completo_N109 - metadados_completo_N109.csv", 
   #                  sep = ",", 
  #                   header = TRUE, 
    #                 stringsAsFactors = FALSE, 
     #                fileEncoding = "UTF-8")

#metadados.all <- read.csv("metadados_all - metadados_all.csv", 
#                          sep = ",", 
#                          header = TRUE, 
#                          stringsAsFactors = FALSE, 
#                          fileEncoding = "UTF-8")
#




#metadata_diet <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Qiime/Analises_Corrigidas/metadata_residual_diet - diet_data.csv", 
#                          sep = ",", 
#                          header = TRUE, 
#                          stringsAsFactors = FALSE, 
#                          fileEncoding = "UTF-8")

#metadados.all <- merge(metadata, metadata_diet, by.x = "Sample.id", by.y = "Sample.id")

#bhei_score_n109 <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/Dieta/BHEI_R_scores_N109.xlsx - Sheet1.csv", 
#                            sep = ",", 
#                            header = TRUE, 
#                            stringsAsFactors = FALSE, 
#                            fileEncoding = "UTF-8")



#metadados.all <- merge(metadados.all, bhei_score_n109, by.x = "Sample.id", by.y = "IdVoluntario")

# Juntar colunas selecionadas de metadados.alpha.all a metadados.all
#metadados.all <- metadados.all %>%
#  left_join(metadados.alpha.all %>%
#              select(Sample.id,
#                     ConsumoGrupo_NOVA_group_1, ConsumoGrupo_NOVA_group_2,
#                     ConsumoGrupo_NOVA_group_3, Percentual_NOVA_group_1,
#                     Percentual_NOVA_group_2, Percentual_NOVA_group_3,
#                     ConsumoCategoria, BMI, TyG, VAI,
#                     QUICKI, METS_IR, TyG_BMI, TyG_WC, WHR),
#            by = "Sample.id")


# salvar .csv
#write.csv(metadados_grupos_saude_binario, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/metadados_grupos_saude_binario_105.csv", row.names = FALSE
)

#=========== Heatmap Alpha Diversidade com Dieta ==============#
#ORGANIZANDO OS DADOS EM UM OBJETO PARA DEPOIS RODAR A CORRELAÇÃO



# Selecionar as colunas numéricas de interesse e renomeá-las para metadados_shannon_selected
metadados.dieta.alpha <- metadados.all.filtrado %>%
  select(shannon_entropy, simpson, pielou_evenness, observed_features, chao1, faith_pd, Region, Region_type, carboidrato_total_g, proteina_g, lipidios_g, fibra_alimentar_g, colesterol_mg, acidos_graxos_saturados_g, acidos_graxos_monoinsaturados_g, acidos_graxos_poliinsaturados_g, acidos_graxos_trans_g, calcio_mg, ferro_mg, sodio_mg, magnesio_mg, fosforo_mg, potassio_mg, manganes_mg, zinco_mg, cobre_mg, selenio_mcg, vitamina_A_RAE_mcg, vitamina_D_mcg, vitamina_E_mg, tiamina_mg, riboflavina_mg, niacina_mg, vitamina_B6_mg, vitamina_B12_mcg, vitamina_C_mg, equivalente_de_folato_mcg, sal_de_adicao_g, acucar_de_adicao_g, BHEI_R_Score_Total)



colnames(metadados.dieta.alpha)

# Renomeando as colunas para remover a palavra "_g"
#colnames(metadados.dieta.alpha) <- gsub("_g$", "", colnames(metadados.dieta.alpha))

# Renomeando as colunas para remover a palavra "_mg"
#colnames(metadados.dieta.alpha) <- gsub("_mg$", "", colnames(metadados.dieta.alpha))

# Renomeando as colunas para remover a palavra "_mcg"
#colnames(metadados.dieta.alpha) <- gsub("_mcg$", "", colnames(metadados.dieta.alpha))


#colnames(metadados.dieta.alpha)


# Transformar rownames de alpha.shannon em coluna Sample.id
#alpha.shannon <- alpha.shannon %>%
#  tibble::rownames_to_column(var = "Sample.id")

# Juntar a coluna shannon_entropy a metadados.all
#metadados.all <- metadados.all %>%
#  left_join(alpha.shannon %>% select(Sample.id, shannon_entropy),
#            by = "Sample.id")

# Converter rownames em coluna Sample.id
#alpha.evenness <- alpha.evenness %>%
#  tibble::rownames_to_column(var = "Sample.id")

# Agora pode fazer o merge ou usar left_join (mais seguro e legível)
#metadados.all <- metadados.all %>%
#  left_join(alpha.evenness %>% select(Sample.id, pielou_evenness),
#            by = "Sample.id")


# Converter rownames em coluna Sample.id
#alpha.chao1 <- alpha.chao1 %>%
#  tibble::rownames_to_column(var = "Sample.id")

# Converter rownames em coluna Sample.id
#alpha.observed.features <- alpha.observed.features %>%
#  tibble::rownames_to_column(var = "Sample.id")


# Converter rownames em coluna Sample.id
#alpha.simpson <- alpha.simpson %>%
#  tibble::rownames_to_column(var = "Sample.id")

# Agora pode fazer o merge ou usar left_join (mais seguro e legível)
#metadados.all <- metadados.all %>%
#  left_join(alpha.evenness %>% select(Sample.id, pielou_evenness),
#            by = "Sample.id")

# Exportar para uma planilha
#write.csv(metadados.all, "metadados_all.csv", row.names = FALSE)





# Juntar todas as tabelas de diversidade
#metadados.all <- merge(metadados.all, alpha.shannon[, c("Sample.id", "shannon_entropy")], by.x = #"Sample.id", by.y = "Sample.id", all.x = TRUE)


#metadados.all <- merge(metadados.all, alpha.evenness[, c("Sample.id", "pielou_evenness")], by.x = #"Sample.id", by.y = "Sample.id", all.x = TRUE)


#metadados.all <- merge(metadados.all, alpha.pd[, c("Sample.id", "faith_pd")], by.x = "Sample.id", by.y = "Sample.id", all.x = TRUE)


#metadados.all <- merge(metadados.all, alpha.observed.features[, c("Sample.id", "observed_features")], by.x = "Sample.id", by.y = "Sample.id", all.x = TRUE)

#metadados.all <- merge(metadados.all, alpha.simpson[, c("Sample.id", "simpson")], by.x = "Sample.id", by.y = "Sample.id", all.x = TRUE)


#metadados.all <- merge(metadados.all, alpha.chao1[, c("Sample.id", "chao1")], by.x = "Sample.id", by.y = "Sample.id", all.x = TRUE)




#metadados.all <- metadados.all %>%
#  left_join(ace_df %>% select(SampleID, ACE), 
#            by = c("Sample.id" = "SampleID"))


#============ Limpando os objetos ASV e criando Phyloseq ===================#
#Para criar phyloseq

# Carregar o objeto!

#asv_table <- fread("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Qiime/exported-asv-table/feature-table.tsv",
#                   sep = "\t", header = TRUE)


#print(colnames(asv_table))



library(data.table)

# Certificar que asv_table está no formato data.table
setDT(asv_table, keep.rownames = TRUE) # Garante que os nomes das linhas sejam preservados como coluna

# Renomear a primeira coluna para garantir que ela não seja removida
colnames(asv_table)[1] <- "Feature.ID"  # Ajuste o nome conforme necessário

# Substituir "-" por "." nos nomes das colunas (caso haja hífens)
setnames(asv_table, old = colnames(asv_table), new = gsub("-", ".", colnames(asv_table)))

# Identificar a primeira coluna (a que contém os IDs das ASVs)
id_column <- "Feature.ID"  # Nome da coluna que contém os IDs

# Selecionar apenas as colunas que terminam com ".F00", mantendo a primeira coluna
cols_to_keep <- c(id_column, grep("\\.F00$", colnames(asv_table), value = TRUE))

# Criar um novo dataframe com a primeira coluna + colunas selecionadas
asv_table_filtered <- asv_table[, ..cols_to_keep]

# Verificar se funcionou
head(asv_table_filtered)



# Encontra a posição da coluna "S40401.F00"
pos_limite <- which(colnames(asv_table_filtered) == "S40371.F00")

# Mantém apenas as colunas até essa posição
asv_table_filtered <- asv_table_filtered[, 1:pos_limite, drop = FALSE]

# Verifica os nomes das colunas restantes
print(colnames(asv_table_filtered))

asv_table <- asv_table_filtered

rm(asv_table_filtered)

asv_table_LIMPA <-asv_table


# Converter para matrix e definir ASV como rownames
# Remover a coluna Feature.ID antes da conversão para matriz
asv_matrix <- as.matrix(asv_table_LIMPA[, -1, with = FALSE]) 

# Definir os nomes das linhas com Feature.ID
rownames(asv_matrix) <- asv_table_LIMPA$Feature.ID 

# Confirmar que os valores dentro da matriz são numéricos
mode(asv_matrix)  # Deve retornar "numeric"


#Agora, confira se os valores são numéricos:
sapply(asv_matrix, class)  # Deve retornar todas as colunas como "numeric"

#Criar OTU table
otu_table <- otu_table(asv_matrix, taxa_are_rows = TRUE)



#criar objeto phyloseq
physeq_obj <- phyloseq(otu_table)


otu_table(physeq_obj)

phy_tree(physeq_obj)  # Deve mostrar a árvore


sample_names(physeq_obj)  # Deve retornar os nomes das amostras

taxa_names(physeq_obj)  # Deve retornar os nomes dos ASVs

rank_names(physeq_obj)  # Verifica se há tabela taxonômica associada



# IDs das ASVs na OTU Table
otu_ids <- taxa_names(physeq_obj)

# IDs na tabela de taxonomia
taxa_ids <- taxonomy$Feature.ID

# IDs das amostras na OTU Table e Metadata
sample_ids_otu <- sample_names(otu_table)
sample_ids_meta <- row.names(metadata)

# Verificações
sum(taxa_ids %in% otu_ids)      # Deve ser igual ao número de linhas da taxonomia
sum(sample_ids_meta %in% sample_ids_otu)  # Deve ser igual ao número de amostras



rownames(taxonomy) <- taxonomy$Feature.ID  # Define os nomes das linhas como IDs das ASVs
taxonomy <- taxonomy[, -1]  # Remove a coluna original "Feature.ID"




# Criar os componentes do phyloseq
otu_table_ps <- otu_table(otu_table, taxa_are_rows = TRUE)
tax_table_ps <- tax_table(as.matrix(taxonomy))

sample_data_ps <- sample_data(metadata)
phy_tree_ps <- phy_tree(tree)

physeq_obj <- phyloseq(otu_table_ps, 
                        tax_table_ps, 
                        phy_tree_ps)


# Verifique os IDs das OTUs na matriz de abundância
otu_ids <- taxa_names(otu_table_ps)

# Verifique os IDs na tabela taxonômica
taxa_ids <- taxa_names(tax_table_ps)

# Verifique os IDs na árvore filogenética
tree_ids <- phy_tree(phy_tree_ps)$tip.label

# Veja quantas OTUs da tabela estão na taxonomia
sum(otu_ids %in% taxa_ids)  # Deve ser igual ao número de OTUs

# Veja quantas OTUs da árvore estão na matriz OTU
sum(otu_ids %in% tree_ids)  # Deve ser igual ao número de OTUs

# Veja quantas OTUs da árvore estão na taxonomia
sum(tree_ids %in% taxa_ids) # Deve ser igual ao número de OTUs



#====================#
#Alpha Diversidade
#===================#


library(vegan)



# Lê a primeira linha para ver os nomes reais das amostras:
first_line <- readLines("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Qiime/Analises_Corrigidas/feature-table - feature-table.tsv", n = 2)[2]
cat(first_line)


abund_table <- read.delim("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Qiime/Analises_Corrigidas/feature-table - feature-table.tsv", 
                          skip = 1, row.names = 1, check.names = FALSE)

# Conferir os nomes das colunas:
head(colnames(abund_table), 10)




# Conferir se importou certo:
dim(abund_table)
head(abund_table[, 1:5])  # primeiras colunas


abund_table_t <- t(abund_table)




#==========================#
#ALPHA DIVERSITY
#=========================#



# Remover as colunas correspondentes
abund_table <- abund_table[, !(colnames(abund_table) %in% c("S10041.F00", "S10232.F00", "S30092.F00", "S30191.F00"))]


# Selecionar as colunas numéricas de interesse e renomeá-las para metadados_shannon_selected

#N=76
cor.inflamatory.alpha <- metadados.all.filtrado %>%
  select("shannon_entropy", "simpson", "pielou_evenness", "observed_features", "chao1", "faith_pd", "IFNGamma" , "IL2", "IL4", "IL6", "IL10", "IL17A", "IFNGamma", "TNF" )

cor.inflamatory.alpha <- na.omit(cor.inflamatory.alpha)

# N=73
cor.HU.alpha <- metadados.all.filtrado %>%
  select("shannon_entropy", "simpson", "pielou_evenness", "observed_features", "chao1", "faith_pd", "INSULINA",  "HOMA.IR","PCR")

cor.HU.alpha <- na.omit(cor.HU.alpha)

# N=97
cor.bioquimicos.alpha <- metadados.all.filtrado %>%
select( "shannon_entropy", "simpson", "pielou_evenness", "observed_features", "chao1", "faith_pd", "HbA1c",  "GLICOSE", "Systolic",  "Diastolic", "COLESTEROL", "LDL", "HDL",  "VLDL" , "TRIGLICERIDES" , "TGO", "TGP", "GGT", "IMC","TyG")

cor.bioquimicos.alpha <- cor.bioquimicos.alpha[!is.na(cor.bioquimicos.alpha$GLICOSE), ]


# N=105
cor.residuais.alpha <- metadados.all.filtrado %>%
  select("shannon_entropy", "simpson", "pielou_evenness", "observed_features", "chao1", "faith_pd",  "carboidrato_total_g", "proteina_g", "lipidios_g", "fibra_alimentar_g", "colesterol_mg", "acidos_graxos_saturados_g", "acidos_graxos_monoinsaturados_g", "acidos_graxos_poliinsaturados_g", "acidos_graxos_trans_g", "calcio_mg", "ferro_mg", "sodio_mg", "magnesio_mg", "fosforo_mg", "potassio_mg", "manganes_mg", "zinco_mg", "cobre_mg", "selenio_mcg", "vitamina_A_RAE_mcg", "vitamina_D_mcg", "vitamina_E_mg", "tiamina_mg", "riboflavina_mg", "niacina_mg", "vitamina_B6_mg", "vitamina_B12_mcg", "vitamina_C_mg", "equivalente_de_folato_mcg", "sal_de_adicao_g", "acucar_de_adicao_g")

cor.residuais.alpha <- na.omit(cor.residuais.alpha)

# N= 95
metadados.NOVA.index.alpha <- metadados.all.filtrado %>%
  select("shannon_entropy", "simpson", "pielou_evenness", "observed_features", "chao1", "faith_pd", "Percentual_NOVA_group_1", "Percentual_NOVA_group_2", "Percentual_NOVA_group_3")

cor.diet.index.alpha <- na.omit(cor.diet.index.alpha)

summary(cor.diet.index.alpha)



#=====================================#
#===== correlacao para citocinas =====#
#=====================================#

library(ggplot2)
library(reshape2)

#=====================Juntar citocinas com HU (Insulina e PCR n73) =====================#

library(dplyr)

metadados.inflamatory.alpha <- metadados.inflamatory.alpha %>%
  left_join(
    metadados.HU.alpha %>% select(Sample.id, INSULINA, PCR),
    by = "Sample.id"
  )


#=== P-value Raw ======#

# Calcular correlações (se ainda não tiver)
n <- ncol(metadados.inflamatory.alpha)
cor_matrix <- matrix(NA, nrow = n, ncol = n)
p_matrix <- matrix(NA, nrow = n, ncol = n)
colnames(cor_matrix) <- rownames(cor_matrix) <- colnames(metadados.inflamatory.alpha)
colnames(p_matrix) <- rownames(p_matrix) <- colnames(metadados.inflamatory.alpha)

for (i in 2:n) {
  for (j in 1:(i - 1)) {
    x <- as.numeric(as.character(metadados.inflamatory.alpha[[i]]))
    y <- as.numeric(as.character(metadados.inflamatory.alpha[[j]]))
    cor_val <- suppressWarnings(cor(x, y, method = "spearman", use = "pairwise.complete.obs"))
    idx <- complete.cases(x, y)
    x_clean <- x[idx]
    y_clean <- y[idx]
    if (length(x_clean) > 2) {
      test <- cor.test(x_clean, y_clean, method = "spearman")
      cor_matrix[i, j] <- cor_val
      p_matrix[i, j] <- test$p.value
    }
  }
}

# Asteriscos com p-valor bruto
asterisks_raw <- ifelse(p_matrix < 0.001, "***",
                        ifelse(p_matrix < 0.01, "**",
                               ifelse(p_matrix < 0.05, "*", "")))

# Derreter para data.frame
library(reshape2)
df_plot_raw <- reshape2::melt(cor_matrix, na.rm = TRUE)
colnames(df_plot_raw) <- c("Var1", "Var2", "cor")
df_plot_raw$p <- reshape2::melt(p_matrix, na.rm = TRUE)[, 3]
df_plot_raw$asterisks <- reshape2::melt(asterisks_raw, na.rm = TRUE)[, 3]



library(ggplot2)
library(reshape2)

# Organizar fatores para garantir a ordem correta dos eixos
df_plot_raw$Var1 <- factor(df_plot_raw$Var1, levels = colnames(cor_matrix))
df_plot_raw$Var2 <- factor(df_plot_raw$Var2, levels = colnames(cor_matrix))

# Manter apenas triângulo inferior sem diagonal
df_plot_raw <- df_plot_raw[as.numeric(df_plot_raw$Var1) > as.numeric(df_plot_raw$Var2), ]

# Plot
p <- ggplot(df_plot_raw, aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = asterisks), size = 3) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, limits = c(-1, 1), name = "Spearman"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Alpha Diversity and Inflamatory Markers (raw p-values) n=76",
    x = "", y = ""
  )


# Salvar em alta resolução
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_cytokines_raw_pvalues.png",
       plot = p, width = 20, height = 15, units = "cm", dpi = 300)






#====== FDR adjusted p-value =======#


#Agora vamos adaptar o seu código para calcular a correlação de Spearman com p-valores ajustados por FDR (False Discovery Rate) — método de Benjamini-Hochberg.


# Número de variáveis
n <- ncol(metadados.inflamatory.alpha)
cor_matrix <- matrix(NA, nrow = n, ncol = n)
p_matrix <- matrix(NA, nrow = n, ncol = n)
colnames(cor_matrix) <- rownames(cor_matrix) <- colnames(metadados.inflamatory.alpha)
colnames(p_matrix) <- rownames(p_matrix) <- colnames(metadados.inflamatory.alpha)

# Loop para calcular correlações e p-valores
for (i in 2:n) {
  for (j in 1:(i - 1)) {
    x <- as.numeric(as.character(metadados.inflamatory.alpha[[i]]))
    y <- as.numeric(as.character(metadados.inflamatory.alpha[[j]]))
    cor_val <- suppressWarnings(cor(x, y, method = "spearman", use = "pairwise.complete.obs"))
    idx <- complete.cases(x, y)
    x_clean <- x[idx]
    y_clean <- y[idx]
    if (length(x_clean) > 2) {
      test <- cor.test(x_clean, y_clean, method = "spearman")
      cor_matrix[i, j] <- cor_val
      p_matrix[i, j] <- test$p.value
    }
  }
}

# Ajuste FDR (Benjamini-Hochberg)
p_adjusted <- matrix(NA, nrow = n, ncol = n)
p_adjusted[lower.tri(p_matrix)] <- p.adjust(p_matrix[lower.tri(p_matrix)], method = "fdr")

# Gera asteriscos
asterisks_fdr <- ifelse(p_adjusted < 0.001, "***",
                        ifelse(p_adjusted < 0.01, "**",
                               ifelse(p_adjusted < 0.05, "*", "")))

# Derreter para data.frame
library(reshape2)
df_plot_fdr <- reshape2::melt(cor_matrix, na.rm = TRUE)
colnames(df_plot_fdr) <- c("Var1", "Var2", "cor")
df_plot_fdr$p <- reshape2::melt(p_adjusted, na.rm = TRUE)[, 3]
df_plot_fdr$asterisks <- reshape2::melt(asterisks_fdr, na.rm = TRUE)[, 3]

# Ordenar fatores para manter layout
df_plot_fdr$Var1 <- factor(df_plot_fdr$Var1, levels = colnames(cor_matrix))
df_plot_fdr$Var2 <- factor(df_plot_fdr$Var2, levels = colnames(cor_matrix))

# Triângulo inferior
df_plot_fdr <- df_plot_fdr[as.numeric(df_plot_fdr$Var1) > as.numeric(df_plot_fdr$Var2), ]

# Plot final
library(ggplot2)
p2 <- ggplot(df_plot_fdr, aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = asterisks), size = 3) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, limits = c(-1, 1), name = "Spearman"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Alpha Diversity and Inflamatory Markers (FDR adjusted) n= 76",
    x = "", y = ""
  )




# Salvar em alta resolução
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_cytokines_FDR_adjusted.png",
       plot = p2, width = 20, height = 15, units = "cm", dpi = 300)



#=====================================#
#===== correlacao para Dietas =====#
#=====================================#


#===== Correlação com p-valores brutos (raw) ===#

library(reshape2)
library(ggplot2)

# Selecionar apenas colunas numéricas relevantes (alpha + dieta)
metadados.dieta.alpha <- metadados.dieta.alpha[, sapply(metadados.dieta.alpha, is.numeric)]

# Calcular matrizes
n <- ncol(metadados.dieta.alpha)
cor_matrix <- matrix(NA, nrow = n, ncol = n)
p_matrix <- matrix(NA, nrow = n, ncol = n)
colnames(cor_matrix) <- rownames(cor_matrix) <- colnames(metadados.dieta.alpha)
colnames(p_matrix) <- rownames(p_matrix) <- colnames(metadados.dieta.alpha)

# Loop de correlação
for (i in 2:n) {
  for (j in 1:(i - 1)) {
    x <- as.numeric(metadados.dieta.alpha[[i]])
    y <- as.numeric(metadados.dieta.alpha[[j]])
    idx <- complete.cases(x, y)
    x_clean <- x[idx]
    y_clean <- y[idx]
    if (length(x_clean) > 2) {
      cor_val <- cor(x_clean, y_clean, method = "spearman")
      test <- cor.test(x_clean, y_clean, method = "spearman")
      cor_matrix[i, j] <- cor_val
      p_matrix[i, j] <- test$p.value
    }
  }
}

# Asteriscos com p bruto
asterisks <- ifelse(p_matrix < 0.001, "***",
                    ifelse(p_matrix < 0.01, "**",
                           ifelse(p_matrix < 0.05, "*", "")))

# Derreter e filtrar triângulo inferior
df_plot <- melt(cor_matrix, na.rm = TRUE)
colnames(df_plot) <- c("Var1", "Var2", "cor")
df_plot$p <- melt(p_matrix, na.rm = TRUE)[, 3]
df_plot$asterisks <- melt(asterisks, na.rm = TRUE)[, 3]
df_plot <- df_plot[as.numeric(df_plot$Var1) > as.numeric(df_plot$Var2), ]

# Heatmap
d <- ggplot(df_plot, aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = asterisks), size = 2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6)) +
  labs(title = "Alpha Diversity and Diet ( Raw p-value) n=105", x = "", y = "")

# Salvar em alta resolução
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_dietXalpha_p_raw.png",
       plot = d, width = 20, height = 15, units = "cm", dpi = 300)

#===== Correlação com p-valores brutos FDR ajusted ===#


# Ajuste FDR
p_adjusted <- matrix(NA, nrow = n, ncol = n)
p_adjusted[lower.tri(p_matrix)] <- p.adjust(p_matrix[lower.tri(p_matrix)], method = "fdr")

# Asteriscos FDR
asterisks_fdr <- ifelse(p_adjusted < 0.001, "***",
                        ifelse(p_adjusted < 0.01, "**",
                               ifelse(p_adjusted < 0.05, "*", "")))

# Derreter e triângulo inferior
df_plot_fdr <- melt(cor_matrix, na.rm = TRUE)
colnames(df_plot_fdr) <- c("Var1", "Var2", "cor")
df_plot_fdr$p <- melt(p_adjusted, na.rm = TRUE)[, 3]
df_plot_fdr$asterisks <- melt(asterisks_fdr, na.rm = TRUE)[, 3]
df_plot_fdr <- df_plot_fdr[as.numeric(df_plot_fdr$Var1) > as.numeric(df_plot_fdr$Var2), ]

# Heatmap FDR
d2 <- ggplot(df_plot_fdr, aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = asterisks), size = 2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6)) +
  labs(title = "Alpha Diversity and Diet (FDR adjusted) n=105", x = "", y = "")


# Salvar em alta resolução
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_dietXalpha_FDRadjusted.png",
       plot = d2, width = 20, height = 15, units = "cm", dpi = 300)


#=====================================#
#===== correlacao com NOVA.Index =====#
#=====================================#

metadados.NOVA.index.alpha <- metadados.NOVA.index.alpha %>%
  filter(!is.na(Percentual_NOVA_group_3))


#===== Correlação com p-valores brutos (raw) ===#

library(reshape2)
library(ggplot2)

# Selecionar apenas colunas numéricas relevantes (alpha + dieta)
metadados.NOVA.index.alpha <- metadados.NOVA.index.alpha[, sapply(metadados.NOVA.index.alpha, is.numeric)]

# Calcular matrizes
n <- ncol(metadados.NOVA.index.alpha)
cor_matrix <- matrix(NA, nrow = n, ncol = n)
p_matrix <- matrix(NA, nrow = n, ncol = n)
colnames(cor_matrix) <- rownames(cor_matrix) <- colnames(metadados.NOVA.index.alpha)
colnames(p_matrix) <- rownames(p_matrix) <- colnames(metadados.NOVA.index.alpha)

# Loop de correlação
for (i in 2:n) {
  for (j in 1:(i - 1)) {
    x <- as.numeric(metadados.NOVA.index.alpha[[i]])
    y <- as.numeric(metadados.NOVA.index.alpha[[j]])
    idx <- complete.cases(x, y)
    x_clean <- x[idx]
    y_clean <- y[idx]
    if (length(x_clean) > 2) {
      cor_val <- cor(x_clean, y_clean, method = "spearman")
      test <- cor.test(x_clean, y_clean, method = "spearman")
      cor_matrix[i, j] <- cor_val
      p_matrix[i, j] <- test$p.value
    }
  }
}

# Asteriscos com p bruto
asterisks <- ifelse(p_matrix < 0.001, "***",
                    ifelse(p_matrix < 0.01, "**",
                           ifelse(p_matrix < 0.05, "*", "")))

# Derreter e filtrar triângulo inferior
df_plot <- melt(cor_matrix, na.rm = TRUE)
colnames(df_plot) <- c("Var1", "Var2", "cor")
df_plot$p <- melt(p_matrix, na.rm = TRUE)[, 3]
df_plot$asterisks <- melt(asterisks, na.rm = TRUE)[, 3]
df_plot <- df_plot[as.numeric(df_plot$Var1) > as.numeric(df_plot$Var2), ]

# Heatmap
e <- ggplot(df_plot, aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = asterisks), size = 2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6)) +
  labs(title = "Alpha Diversity and NOVA INDEX (Raw p-value) n=95", x = "", y = "")

# Salvar em alta resolução
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_NOVAXalpha_p_raw.png",
       plot = e, width = 20, height = 15, units = "cm", dpi = 300)

#===== Correlação com p-valores brutos FDR ajusted ===#


# Ajuste FDR
p_adjusted <- matrix(NA, nrow = n, ncol = n)
p_adjusted[lower.tri(p_matrix)] <- p.adjust(p_matrix[lower.tri(p_matrix)], method = "fdr")

# Asteriscos FDR
asterisks_fdr <- ifelse(p_adjusted < 0.001, "***",
                        ifelse(p_adjusted < 0.01, "**",
                               ifelse(p_adjusted < 0.05, "*", "")))

# Derreter e triângulo inferior
df_plot_fdr <- melt(cor_matrix, na.rm = TRUE)
colnames(df_plot_fdr) <- c("Var1", "Var2", "cor")
df_plot_fdr$p <- melt(p_adjusted, na.rm = TRUE)[, 3]
df_plot_fdr$asterisks <- melt(asterisks_fdr, na.rm = TRUE)[, 3]
df_plot_fdr <- df_plot_fdr[as.numeric(df_plot_fdr$Var1) > as.numeric(df_plot_fdr$Var2), ]

# Heatmap FDR
e2 <- ggplot(df_plot_fdr, aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = asterisks), size = 2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6)) +
  labs(title = "Alpha Diversity and NOVA INDEX (FDR adjusted) n=95", x = "", y = "")


# Salvar em alta resolução
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_NOVAXalpha_FDRadjusted.png",
       plot = e2, width = 20, height = 15, units = "cm", dpi = 300)


#=============================================#
#===== correlacao com Metabolic Markers =====#
#=============================================#

#===== Correlação com p-valores brutos (raw) ===#

library(reshape2)
library(ggplot2)

# Selecionar apenas colunas numéricas relevantes (alpha + dieta)
metadados.HU.alpha <- metadados.HU.alpha[, sapply(metadados.HU.alpha, is.numeric)]

# Calcular matrizes
n <- ncol(metadados.HU.alpha)
cor_matrix <- matrix(NA, nrow = n, ncol = n)
p_matrix <- matrix(NA, nrow = n, ncol = n)
colnames(cor_matrix) <- rownames(cor_matrix) <- colnames(metadados.HU.alpha)
colnames(p_matrix) <- rownames(p_matrix) <- colnames(metadados.HU.alpha)

# Loop de correlação
for (i in 2:n) {
  for (j in 1:(i - 1)) {
    x <- as.numeric(metadados.HU.alpha[[i]])
    y <- as.numeric(metadados.HU.alpha[[j]])
    idx <- complete.cases(x, y)
    x_clean <- x[idx]
    y_clean <- y[idx]
    if (length(x_clean) > 2) {
      cor_val <- cor(x_clean, y_clean, method = "spearman")
      test <- cor.test(x_clean, y_clean, method = "spearman")
      cor_matrix[i, j] <- cor_val
      p_matrix[i, j] <- test$p.value
    }
  }
}

# Asteriscos com p bruto
asterisks <- ifelse(p_matrix < 0.001, "***",
                    ifelse(p_matrix < 0.01, "**",
                           ifelse(p_matrix < 0.05, "*", "")))

# Derreter e filtrar triângulo inferior
df_plot <- melt(cor_matrix, na.rm = TRUE)
colnames(df_plot) <- c("Var1", "Var2", "cor")
df_plot$p <- melt(p_matrix, na.rm = TRUE)[, 3]
df_plot$asterisks <- melt(asterisks, na.rm = TRUE)[, 3]
df_plot <- df_plot[as.numeric(df_plot$Var1) > as.numeric(df_plot$Var2), ]

# Heatmap
f <- ggplot(df_plot, aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = asterisks), size = 2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6)) +
  labs(title = "Alpha Diversity and Metabolic Markers 1 (p-valor bruto) - n=73", x = "", y = "")

# Salvar em alta resolução
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_MetabolicMarkers1Xalpha_p_raw.png",
       plot = f, width = 20, height = 15, units = "cm", dpi = 300)

#===== Correlação com p-valores brutos FDR ajusted ===#


# Ajuste FDR
p_adjusted <- matrix(NA, nrow = n, ncol = n)
p_adjusted[lower.tri(p_matrix)] <- p.adjust(p_matrix[lower.tri(p_matrix)], method = "fdr")

# Asteriscos FDR
asterisks_fdr <- ifelse(p_adjusted < 0.001, "***",
                        ifelse(p_adjusted < 0.01, "**",
                               ifelse(p_adjusted < 0.05, "*", "")))

# Derreter e triângulo inferior
df_plot_fdr <- melt(cor_matrix, na.rm = TRUE)
colnames(df_plot_fdr) <- c("Var1", "Var2", "cor")
df_plot_fdr$p <- melt(p_adjusted, na.rm = TRUE)[, 3]
df_plot_fdr$asterisks <- melt(asterisks_fdr, na.rm = TRUE)[, 3]
df_plot_fdr <- df_plot_fdr[as.numeric(df_plot_fdr$Var1) > as.numeric(df_plot_fdr$Var2), ]

# Heatmap FDR
f2 <- ggplot(df_plot_fdr, aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = asterisks), size = 2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6)) +
  labs(title = "Alpha Diversity and Metabolic Markers 1 (FDR adjusted) - n=73" , x = "", y = "")


# Salvar em alta resolução
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_MetabolicMarkers1_FDRadjusted.png",
       plot = f2, width = 20, height = 15, units = "cm", dpi = 300)


#========================================================================#
#===== correlacao com Alpha Diversity and Clinical Biomarkers n=97 =====#
#=======================================================================#

#===== Correlação com p-valores brutos (raw) ===#

library(reshape2)
library(ggplot2)

# Selecionar apenas colunas numéricas relevantes (alpha + dieta)
metadados.bioquimicos.alpha <- metadados.bioquimicos.alpha[, sapply(metadados.bioquimicos.alpha, is.numeric)]

# Calcular matrizes
n <- ncol(metadados.bioquimicos.alpha)
cor_matrix <- matrix(NA, nrow = n, ncol = n)
p_matrix <- matrix(NA, nrow = n, ncol = n)
colnames(cor_matrix) <- rownames(cor_matrix) <- colnames(metadados.bioquimicos.alpha)
colnames(p_matrix) <- rownames(p_matrix) <- colnames(metadados.bioquimicos.alpha)

# Loop de correlação
for (i in 2:n) {
  for (j in 1:(i - 1)) {
    x <- as.numeric(metadados.bioquimicos.alpha[[i]])
    y <- as.numeric(metadados.bioquimicos.alpha[[j]])
    idx <- complete.cases(x, y)
    x_clean <- x[idx]
    y_clean <- y[idx]
    if (length(x_clean) > 2) {
      cor_val <- cor(x_clean, y_clean, method = "spearman")
      test <- cor.test(x_clean, y_clean, method = "spearman")
      cor_matrix[i, j] <- cor_val
      p_matrix[i, j] <- test$p.value
    }
  }
}

# Asteriscos com p bruto
asterisks <- ifelse(p_matrix < 0.001, "***",
                    ifelse(p_matrix < 0.01, "**",
                           ifelse(p_matrix < 0.05, "*", "")))

# Derreter e filtrar triângulo inferior
df_plot <- melt(cor_matrix, na.rm = TRUE)
colnames(df_plot) <- c("Var1", "Var2", "cor")
df_plot$p <- melt(p_matrix, na.rm = TRUE)[, 3]
df_plot$asterisks <- melt(asterisks, na.rm = TRUE)[, 3]
df_plot <- df_plot[as.numeric(df_plot$Var1) > as.numeric(df_plot$Var2), ]

# Heatmap
g <- ggplot(df_plot, aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = asterisks), size = 2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6)) +
  labs(title = "Alpha Diversity and Clinical Biomarkers (raw p-value) - n=97", x = "", y = "")

# Salvar em alta resolução
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_Clinical BiomarkersXalpha_p_raw.png",
       plot = g, width = 20, height = 15, units = "cm", dpi = 300)

#===== Correlação com p-valores brutos FDR ajusted ===#


# Ajuste FDR
p_adjusted <- matrix(NA, nrow = n, ncol = n)
p_adjusted[lower.tri(p_matrix)] <- p.adjust(p_matrix[lower.tri(p_matrix)], method = "fdr")

# Asteriscos FDR
asterisks_fdr <- ifelse(p_adjusted < 0.001, "***",
                        ifelse(p_adjusted < 0.01, "**",
                               ifelse(p_adjusted < 0.05, "*", "")))

# Derreter e triângulo inferior
df_plot_fdr <- melt(cor_matrix, na.rm = TRUE)
colnames(df_plot_fdr) <- c("Var1", "Var2", "cor")
df_plot_fdr$p <- melt(p_adjusted, na.rm = TRUE)[, 3]
df_plot_fdr$asterisks <- melt(asterisks_fdr, na.rm = TRUE)[, 3]
df_plot_fdr <- df_plot_fdr[as.numeric(df_plot_fdr$Var1) > as.numeric(df_plot_fdr$Var2), ]

# Heatmap FDR
g2 <- ggplot(df_plot_fdr, aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = asterisks), size = 2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6)) +
  labs(title = "Alpha Diversity and Metabolic Markers 2 (FDR adjusted) - n=97" , x = "", y = "")


# Salvar em alta resolução
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_Clinical Biomarkers_FDRadjusted.png",
       plot = g2, width = 20, height = 15, units = "cm", dpi = 300)




#==============================================================================#

#                              Scatterplots 

#=============================================================================#

library(ggplot2)
library(ggpubr)
library(ggplot2)
library(gridExtra)

# Crie os scatterplots individuais

#==== TyG ====#

# Gráfico com rho e p-valor
p1 <- ggplot(metadados.bioquimicos.alpha, aes(x = TyG, y = chao1)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~")), label.x = 4.5, label.y = 400) +
  labs(title = "chao1 vs TyG", x = "TyG", y = "chao1") +
  theme_minimal()

  # Gráfico com rho e p-valor
  p2 <- ggplot(metadados.bioquimicos.alpha, aes(x = TyG, y = observed_features)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~")), label.x = 4.5, label.y = 400) +
  labs(title = "Observed Features vs TyG", x = "TyG", y = "Observed Features") +
  theme_minimal()

# Salvar em PNG com alta resolução
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/TyG_vs_alpha.png",
       arrangeGrob(p1, p2, ncol = 2),
       width = 45, height = 15, units = "cm", dpi = 300)


#==== IMC =====#

p1 <- ggplot(metadados.bioquimicos.alpha, aes(x = IMC, y = chao1)) +
  geom_point() + geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~")),
           label.x = 20, label.y = 400) +
  labs(title = "chao1 vs IMC", x = "IMC", y = "chao1") + theme_minimal()

p2 <- ggplot(metadados.bioquimicos.alpha, aes(x = IMC, y = observed_features)) +
  geom_point() + geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~")),
           label.x = 20, label.y = 400) +
  labs(title = "Observed Features vs IMC", x = "IMC", y = "Observed Features") + theme_minimal()

ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/IMC_vs_alpha.png",
       arrangeGrob(p1, p2, ncol = 2), width = 45, height = 15, units = "cm", dpi = 300)



 #=== TGP ===#

library(ggplot2)
library(ggpubr)
library(gridExtra)

# Gráfico 1
p1 <- ggplot(metadados.bioquimicos.alpha, aes(x = TGP, y = faith_pd)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "spearman", 
           aes(label = paste(..r.label.., ..p.label.., sep = "~")), 
           label.x = 10, label.y = 25) +
  labs(title = "Faith PD vs TGP", x = "TGP", y = "Faith PD") +
  theme_minimal()

# Gráfico 2
p2 <- ggplot(metadados.bioquimicos.alpha, aes(x = TGP, y = chao1)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "spearman", 
           aes(label = paste(..r.label.., ..p.label.., sep = "~")), 
           label.x = 10, label.y = 400) +
  labs(title = "Chao1 vs TGP", x = "TGP", y = "Chao1") +
  theme_minimal()

# Gráfico 3
p3 <- ggplot(metadados.bioquimicos.alpha, aes(x = TGP, y = observed_features)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "spearman", 
           aes(label = paste(..r.label.., ..p.label.., sep = "~")), 
           label.x = 10, label.y = 400) +
  labs(title = "Observed Features vs TGP", x = "TGP", y = "Observed Features") +
  theme_minimal()

# Gráfico 4
p4 <- ggplot(metadados.bioquimicos.alpha, aes(x = TGP, y = shannon_entropy)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "spearman", 
           aes(label = paste(..r.label.., ..p.label.., sep = "~")), 
           label.x = 10, label.y = 7) +
  labs(title = "Shannon Entropy vs TGP", x = "TGP", y = "Shannon Entropy") +
  theme_minimal()

# Salvar em PNG
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/TGP_vs_alpha_diversity.png",
       arrangeGrob(p1, p2, p3, p4, ncol = 2),
       width = 90, height = 50, units = "cm", dpi = 300)


#====  GGT =====#

library(ggplot2)
library(ggpubr)

# Gráfico único: GGT vs shannon_entropy
p <- ggplot(metadados.bioquimicos.alpha, aes(x = GGT, y = shannon_entropy)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "spearman",
           aes(label = paste(..r.label.., ..p.label.., sep = "~")),
           label.x = max(metadados.bioquimicos.alpha$GGT, na.rm = TRUE) * 0.6,
           label.y = max(metadados.bioquimicos.alpha$shannon_entropy, na.rm = TRUE) * 0.95) +
  labs(title = "Shannon entropy vs GGT", x = "GGT", y = "Shannon entropy") +
  theme_minimal()

# Salvar em PNG
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/GGT_vs_shannon.png",
       plot = p, width = 22.5, height = 15, units = "cm", dpi = 300)


#=== TRIGLICERIDES===#

library(ggplot2)
library(ggpubr)

# Gráfico 1: TRIGLICERIDES vs faith_pd
p1 <- ggplot(metadados.bioquimicos.alpha, aes(x = TRIGLICERIDES, y = faith_pd)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "spearman",
           aes(label = paste(..r.label.., ..p.label.., sep = "~")),
           label.x = max(metadados.bioquimicos.alpha$TRIGLICERIDES, na.rm = TRUE) * 0.6,
           label.y = max(metadados.bioquimicos.alpha$faith_pd, na.rm = TRUE) * 0.95) +
  labs(title = "Faith PD vs Triglycerides", x = "Triglycerides", y = "Faith PD") +
  theme_minimal()

# Gráfico 2: TRIGLICERIDES vs chao1
p2 <- ggplot(metadados.bioquimicos.alpha, aes(x = TRIGLICERIDES, y = chao1)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "spearman",
           aes(label = paste(..r.label.., ..p.label.., sep = "~")),
           label.x = max(metadados.bioquimicos.alpha$TRIGLICERIDES, na.rm = TRUE) * 0.6,
           label.y = max(metadados.bioquimicos.alpha$chao1, na.rm = TRUE) * 0.95) +
  labs(title = "Chao1 vs Triglycerides", x = "Triglycerides", y = "Chao1") +
  theme_minimal()

# Gráfico 3: TRIGLICERIDES vs observed_features
p3 <- ggplot(metadados.bioquimicos.alpha, aes(x = TRIGLICERIDES, y = observed_features)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "spearman",
           aes(label = paste(..r.label.., ..p.label.., sep = "~")),
           label.x = max(metadados.bioquimicos.alpha$TRIGLICERIDES, na.rm = TRUE) * 0.6,
           label.y = max(metadados.bioquimicos.alpha$observed_features, na.rm = TRUE) * 0.95) +
  labs(title = "Observed Features vs Triglycerides", x = "Triglycerides", y = "Observed Features") +
  theme_minimal()

# Salvar todos em um único PNG
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Triglycerides_vs_alpha.png",
       plot = ggpubr::ggarrange(p1, p2, p3, ncol = 3),
       width = 60, height = 15, units = "cm", dpi = 300)


#===== VLDL ======#

library(ggplot2)
library(ggpubr)
library(gridExtra)

# chao1 vs VLDL
p1 <- ggplot(metadados.bioquimicos.alpha, aes(x = VLDL, y = chao1)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "spearman",
           aes(label = paste(..r.label.., ..p.label.., sep = "~")),
           label.x = max(metadados.bioquimicos.alpha$VLDL, na.rm = TRUE) * 0.6,
           label.y = max(metadados.bioquimicos.alpha$chao1, na.rm = TRUE) * 0.95) +
  labs(title = "Chao1 vs VLDL", x = "VLDL", y = "Chao1") +
  theme_minimal()

# observed_features vs VLDL
p2 <- ggplot(metadados.bioquimicos.alpha, aes(x = VLDL, y = observed_features)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "spearman",
           aes(label = paste(..r.label.., ..p.label.., sep = "~")),
           label.x = max(metadados.bioquimicos.alpha$VLDL, na.rm = TRUE) * 0.6,
           label.y = max(metadados.bioquimicos.alpha$observed_features, na.rm = TRUE) * 0.95) +
  labs(title = "Observed Features vs VLDL", x = "VLDL", y = "Observed Features") +
  theme_minimal()

# Salvar os dois gráficos em um único PNG
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/VLDL_vs_alpha.png",
       plot = grid.arrange(p1, p2, ncol = 2),
       width = 45, height = 15, units = "cm", dpi = 300)




#=================================#
#Alpha com TyG - Tercis
#=================================#

metadados.bioquimicos.alpha$TyG_tertile <- cut(
  metadados.bioquimicos.alpha$TyG,
  breaks = quantile(metadados.bioquimicos.alpha$TyG, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
  labels = c("Low", "Medium", "High"),
  include.lowest = TRUE
)

table(metadados.bioquimicos.alpha$TyG_tertile)

quantile(metadados.bioquimicos.alpha$TyG, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)


# Gráfico Chao1 x TyG tercis
p1 <- ggplot(metadados.bioquimicos.alpha, aes(x = TyG_tertile, y = chao1, fill = TyG_tertile)) +
  geom_boxplot() +
  stat_compare_means(method = "kruskal.test") +
  labs(
    title = "Chao1 by TyG Tertiles",
    subtitle = "Low: 3.98–4.49 | Medium: 4.50–4.76 | High: 4.77–9.59",
    x = "TyG Tertiles", y = "Chao1"
  ) +
  theme_minimal()

# Observed Features
p2 <- ggplot(metadados.bioquimicos.alpha, aes(x = TyG_tertile, y = observed_features, fill = TyG_tertile)) +
  geom_boxplot() +
  stat_compare_means(method = "kruskal.test") +
  labs(
    title = "Observed Features by TyG Tertiles",
    subtitle = "Low: 3.98–4.49 | Medium: 4.50–4.76 | High: 4.77–9.59",
    x = "TyG Tertiles", y = "Observed Features"
  ) +
  theme_minimal()

# Salvar
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Box_plot_TyG_tertiles_alpha.png", arrangeGrob(p1, p2, ncol = 2), width = 12, height = 6, dpi = 300)



#====================================#
#Alpha com TyG - 2 grupos pela mediana
#====================================#

library(ggplot2)
library(dplyr)
library(ggpubr)

# Criar grupo binário com base na mediana de TyG
metadados.bioquimicos.alpha <- metadados.bioquimicos.alpha %>%
  mutate(TyG_group = ifelse(TyG <= median(TyG, na.rm = TRUE), "Low", "High"))

# Chao1 vs TyG_group
p1 <- ggplot(metadados.bioquimicos.alpha, aes(x = TyG_group, y = chao1, fill = TyG_group)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label.y = max(metadados.bioquimicos.alpha$chao1, na.rm = TRUE) + 30) +
  labs(title = "Chao1 by TyG Group (Median Split)", x = "TyG Group", y = "Chao1") +
  scale_fill_manual(values = c("Low" = "#E57373", "High" = "#64B5F6")) +
  theme_minimal()

# Observed Features vs TyG_group
p2 <- ggplot(metadados.bioquimicos.alpha, aes(x = TyG_group, y = observed_features, fill = TyG_group)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label.y = max(metadados.bioquimicos.alpha$observed_features, na.rm = TRUE) + 30) +
  labs(title = "Observed Features by TyG Group (Median Split)", x = "TyG Group", y = "Observed Features") +
  scale_fill_manual(values = c("Low" = "#E57373", "High" = "#64B5F6")) +
  theme_minimal()

# Salvar os dois gráficos lado a lado
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Box_plot_TyG_MedianSplit_vs_alpha.png",
       plot = ggpubr::ggarrange(p1, p2, ncol = 2),
       width = 40, height = 15, units = "cm", dpi = 300)



#==================================#

#     Alpha TyG com Extremos

#===================================#




library(ggplot2)
library(ggpubr)
library(dplyr)

# Filtra apenas os grupos Low e High
dados_filtrados <- metadados.bioquimicos.alpha %>%
  filter(TyG_tertile %in% c("Low", "High"))

# Boxplot Chao1
p1 <- ggplot(dados_filtrados, aes(x = TyG_tertile, y = chao1, fill = TyG_tertile)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label.y = max(dados_filtrados$chao1, na.rm = TRUE) + 30) +
  labs(title = "Chao1 by TyG (Tertile 1 vs 3) n=66", x = "TyG Group", y = "Chao1") +
  scale_fill_manual(values = c("Low" = "#E57373", "High" = "#64B5F6")) +
  theme_minimal()

# Boxplot Chao1
p2 <- ggplot(dados_filtrados, aes(x = TyG_tertile, y = observed_features, fill = TyG_tertile)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label.y = max(dados_filtrados$chao1, na.rm = TRUE) + 30) +
  labs(title = "Observed Features by TyG (Tertile 1 vs 3) n=66", x = "TyG Group", y = "Observed Features") +
  scale_fill_manual(values = c("Low" = "#E57373", "High" = "#64B5F6")) +
  theme_minimal()

# Salvar os dois gráficos lado a lado
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Boxplot_TyG_T1vsT3.png",
       plot = ggpubr::ggarrange(p1, p2, ncol = 2),
       width = 40, height = 15, units = "cm", dpi = 300)




#=================================#
#Alpha com BHEI-R - Tercis
#=================================#

metadados.dieta.alpha$BHEI_tertile <- cut(
  metadados.dieta.alpha$BHEI_R_Score_Total,
  breaks = quantile(metadados.dieta.alpha$BHEI_R_Score_Total, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
  labels = c("Low", "Medium", "High"),
  include.lowest = TRUE
)

quantile(metadados.dieta.alpha$BHEI_R_Score_Total, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

table(metadados.dieta.alpha$BHEI_tertile)

# Gráfico Chao1 x BHEI tercis
p1 <- ggplot(metadados.dieta.alpha, aes(x = BHEI_tertile, y = chao1, fill = BHEI_tertile)) +
  geom_boxplot() +
  stat_compare_means(method = "kruskal.test") +
  labs(
    title = "Chao1 by BHEI Tertiles",
    subtitle = "Low: 21.49–49.21 | Medium: 49.22–58.40 | High: 58.41–77.80",
    x = "BHEI Tertiles", y = "Chao1"
  ) +
  theme_minimal()

# Observed Features
p2 <- ggplot(metadados.dieta.alpha, aes(x = BHEI_tertile, y = observed_features, fill = BHEI_tertile)) +
  geom_boxplot() +
  stat_compare_means(method = "kruskal.test") +
  labs(
    title = "Observed Features by BHEI Tertiles",
    subtitle = "Low: 21.49–49.21 | Medium: 49.22–58.40 | High: 58.41–77.80",
    x = "BHEI Tertiles", y = "Observed Features"
  ) +
  theme_minimal()
# Salvar
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Boxplot_BHEI_tertiles_alpha.png", arrangeGrob(p1, p2, ncol = 2), width = 12, height = 6, dpi = 300)






#=================================#
#Alpha com BHEI-R Extremos
#=================================#


library(ggplot2)
library(ggpubr)
library(dplyr)

# Filtra apenas os grupos Low e High
dados_filtrados <- metadados.dieta.alpha %>%
  filter(BHEI_tertile %in% c("Low", "High"))

# Boxplot Chao1
p1 <- ggplot(dados_filtrados, aes(x = BHEI_tertile, y = chao1, fill = BHEI_tertile)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label.y = max(dados_filtrados$chao1, na.rm = TRUE) + 30) +
  labs(title = "Chao1 by BHEI (Tertile 1 vs 3) n=70", subtitle = "Low: 21.49–49.20 | High: 58.40–77.80", x = "BHEI Group", y = "Chao1") +
  scale_fill_manual(values = c("Low" = "#E57373", "High" = "#64B5F6")) +
  theme_minimal()

# Boxplot Chao1
p2 <- ggplot(dados_filtrados, aes(x = BHEI_tertile, y = observed_features, fill = BHEI_tertile)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label.y = max(dados_filtrados$chao1, na.rm = TRUE) + 30) +
  labs(title = "Observed Features by BHEI (Tertile 1 vs 3) n=70", subtitle = "Low: 21.49–49.20 | High: 58.40–77.80", x = "BHEI Group", y = "Observed Features") +
  scale_fill_manual(values = c("Low" = "#E57373", "High" = "#64B5F6")) +
  theme_minimal()

# Salvar os dois gráficos lado a lado
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Boxplot_BHEI_T1vsT3.png",
       plot = ggpubr::ggarrange(p1, p2, ncol = 2),
       width = 40, height = 15, units = "cm", dpi = 300)


#===================================#
#     Dieta Residual Tercis
#===================================#

metadados.dieta.alpha <- metadados.dieta.alpha %>%
  mutate(
    tercil_saturados = ntile(acidos_graxos_saturados_g, 3),
    tercil_trans = ntile(acidos_graxos_trans_g, 3),
    tercil_colesterol = ntile(colesterol_mg, 3)
  ) %>%
  mutate(
    tercil_saturados = factor(tercil_saturados, labels = c("Low", "Medium", "High")),
    tercil_trans = factor(tercil_trans, labels = c("Low", "Medium", "High")),
    tercil_colesterol = factor(tercil_colesterol, labels = c("Low", "Medium", "High"))
  )



#===================================================#
#             Alpha e Saturated Fat
#===================================================#



library(ggpubr)
library(dplyr)

# Função para gerar gráfico para uma variável e índice
make_plot <- function(df, xvar, yvar, ylab, title_prefix, palette = "viridis", ypos_start = 0.9) {
  # Kruskal-Wallis global
  p_kw <- kruskal.test(reformulate(xvar, yvar), data = df)$p.value
  
  # Comparações entre pares com FDR
  pares <- compare_means(as.formula(paste(yvar, "~", xvar)),
                         data = df, method = "wilcox.test", p.adjust.method = "fdr") %>%
    filter(p.adj <= 0.05)
  
  # Adicionar posições Y para os p-valor dos pares
  if (nrow(pares) > 0) {
    max_y <- max(df[[yvar]], na.rm = TRUE)
    pares$y.position <- seq(from = max_y + 0.2, by = 0.2, length.out = nrow(pares))
  }
  # Plot
  p <- ggboxplot(df, x = xvar, y = yvar, fill = xvar, palette = palette) +
    labs(
      title = paste0(title_prefix, " (Kruskal-Wallis p = ", signif(p_kw, 3), ")"),
      x = "Saturated Fat Intake (terciles)", y = ylab
    ) +
    theme_minimal()
  
  if (nrow(pares) > 0) {
    p <- p + stat_pvalue_manual(pares, label = "p.signif", tip.length = 0.01)
  }
  
  return(p)
}

# Criar os 4 gráficos
p1_saturado <- make_plot(metadados.dieta.alpha, "tercil_saturados", "shannon_entropy", "Shannon Entropy", "A. Shannon", ypos_start = 7.0)
p2_saturado <- make_plot(metadados.dieta.alpha, "tercil_saturados", "pielou_evenness", "Pielou Evenness", "B. Pielou", ypos_start = 0.9)
p3_saturado <- make_plot(metadados.dieta.alpha, "tercil_saturados", "chao1", "Chao1 Richness", "C. Chao1", ypos_start = 420)
p4_saturado <- make_plot(metadados.dieta.alpha, "tercil_saturados", "faith_pd", "Faith's PD", "D. Faith", ypos_start = 25)

# Juntar os plots
painel_saturado <- ggarrange(p1_saturado, p2_saturado, p3_saturado, p4_saturado,
                             ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

# Salvar
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/painel_saturated_fat.png", painel_saturado, width = 12, height = 8, dpi = 300)


#============#
Gordura Trans
#============#

# Kruskal-Wallis geral
kw_shannon <- kruskal.test(shannon_entropy ~ tercil_trans, data = metadados.dieta.alpha)$p.value
kw_pielou  <- kruskal.test(pielou_evenness ~ tercil_trans, data = metadados.dieta.alpha)$p.value
kw_chao1   <- kruskal.test(chao1 ~ tercil_trans, data = metadados.dieta.alpha)$p.value
kw_faith   <- kruskal.test(faith_pd ~ tercil_trans, data = metadados.dieta.alpha)$p.value

comparacoes_trans_shannon <- compare_means(shannon_entropy ~ tercil_trans, data = metadados.dieta.alpha, method = "wilcox.test", p.adjust.method = "fdr") %>%
  filter(p.adj <= 0.05)

# Se ainda houver pelo menos uma comparação significativa:
if (nrow(comparacoes_trans_shannon) > 0) {
  comparacoes_trans_shannon$y.position <- seq(7.1, 7.1 + 0.2 * (nrow(comparacoes_trans_shannon) - 1), by = 0.2)
}



# Comparações par a par com FDR


comparacoes_trans_pielou <- compare_means(pielou_evenness ~ tercil_trans, 
                                          data = metadados.dieta.alpha, 
                                          method = "wilcox.test", 
                                          p.adjust.method = "fdr") %>%
  filter(p.adj <= 0.05)

if (nrow(comparacoes_trans_pielou) > 0) {
  comparacoes_trans_pielou$y.position <- seq(0.88, 0.88 + 0.03 * (nrow(comparacoes_trans_pielou) - 1), by = 0.03)
}


comparacoes_trans_chao1 <- compare_means(chao1 ~ tercil_trans, 
                                         data = metadados.dieta.alpha, 
                                         method = "wilcox.test", 
                                         p.adjust.method = "fdr") %>%
  filter(p.adj <= 0.05)

if (nrow(comparacoes_trans_chao1) > 0) {
  comparacoes_trans_chao1$y.position <- seq(420, 420 + 20 * (nrow(comparacoes_trans_chao1) - 1), by = 20)
}

comparacoes_trans_faith <- compare_means(faith_pd ~ tercil_trans, 
                                         data = metadados.dieta.alpha, 
                                         method = "wilcox.test", 
                                         p.adjust.method = "fdr") %>%
  filter(p.adj <= 0.05)

if (nrow(comparacoes_trans_faith) > 0) {
  comparacoes_trans_faith$y.position <- seq(24, 24 + 2 * (nrow(comparacoes_trans_faith) - 1), by = 2)
}


# Plots
p1 <- ggboxplot(metadados.dieta.alpha, x = "tercil_trans", y = "shannon_entropy",
                fill = "tercil_trans", palette = "viridis") +
  labs(title = paste0("A. Shannon Entropy (Kruskal-Wallis p = ", signif(kw_shannon, 3), ")"),
       x = "Trans Fat Intake", y = "Shannon Entropy") +
  stat_pvalue_manual(comparacoes_trans_shannon, label = "p.signif", tip.length = 0.01) +
  theme_minimal()

p2 <- ggboxplot(metadados.dieta.alpha, x = "tercil_trans", y = "pielou_evenness",
                fill = "tercil_trans", palette = "viridis") +
  labs(title = paste0("B. Pielou Index (Kruskal-Wallis p = ", signif(kw_pielou, 3), ")"),
       x = "Trans Fat Intake", y = "Pielou Evenness") +
  stat_pvalue_manual(comparacoes_trans_pielou, label = "p.signif", tip.length = 0.01) +
  theme_minimal()

p3 <- ggboxplot(metadados.dieta.alpha, x = "tercil_trans", y = "chao1",
                fill = "tercil_trans", palette = "viridis") +
  labs(title = paste0("C. Chao1 Richness (Kruskal-Wallis p = ", signif(kw_chao1, 3), ")"),
       x = "Trans Fat Intake", y = "Chao1 Richness") +
  stat_pvalue_manual(comparacoes_trans_chao1, label = "p.signif", tip.length = 0.01) +
  theme_minimal()

p4 <- ggboxplot(metadados.dieta.alpha, x = "tercil_trans", y = "faith_pd",
                fill = "tercil_trans", palette = "viridis") +
  labs(title = paste0("D. Faith's PD (Kruskal-Wallis p = ", signif(kw_faith, 3), ")"),
       x = "Trans Fat Intake", y = "Faith’s Phylogenetic Diversity") +
  stat_pvalue_manual(comparacoes_trans_faith, label = "p.signif", tip.length = 0.01) +
  theme_minimal()

# Juntar
painel_trans <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/painel_trans_significativo.png", painel_trans, width = 12, height = 8, dpi = 300)



#=====================================#
#Colesterol
#====================================#

# Função para comparar com FDR e retornar apenas significativos com posição y
get_comparacoes <- function(var, y_pos) {
  comp <- compare_means(as.formula(paste0(var, " ~ tercil_colesterol")),
                        data = metadados.dieta.alpha,
                        method = "wilcox.test", p.adjust.method = "fdr")
  comp_sig <- comp %>% filter(p.adj <= 0.05)
  if (nrow(comp_sig) > 0) {
    comp_sig$y.position <- y_pos[1:nrow(comp_sig)]
    return(comp_sig)
  } else {
    return(NULL)
  }
}

# Kruskal-Wallis geral
kw_colesterol_shannon <- kruskal.test(shannon_entropy ~ tercil_colesterol, data = metadados.dieta.alpha)$p.value
kw_colesterol_pielou  <- kruskal.test(pielou_evenness ~ tercil_colesterol, data = metadados.dieta.alpha)$p.value
kw_colesterol_chao1   <- kruskal.test(chao1 ~ tercil_colesterol, data = metadados.dieta.alpha)$p.value
kw_colesterol_faith   <- kruskal.test(faith_pd ~ tercil_colesterol, data = metadados.dieta.alpha)$p.value

# Comparações com Wilcoxon
comparacoes_col_shannon <- get_comparacoes("shannon_entropy", c(7.2, 7.4, 7.6))
comparacoes_col_pielou  <- get_comparacoes("pielou_evenness", c(0.88, 0.91, 0.94))
comparacoes_col_chao1   <- get_comparacoes("chao1", c(420, 440, 460))
comparacoes_col_faith   <- get_comparacoes("faith_pd", c(24, 26, 28))

# Gráfico Shannon
p1_col <- ggboxplot(metadados.dieta.alpha, x = "tercil_colesterol", y = "shannon_entropy",
                    fill = "tercil_colesterol", palette = "viridis") +
  labs(title = paste0("A. Shannon Entropy (Kruskal-Wallis p = ", signif(kw_colesterol_shannon, 3), ")"),
       x = "Cholesterol Intake (tercile)", y = "Shannon Entropy") +
  theme_minimal()
if (!is.null(comparacoes_col_shannon)) {
  p1_col <- p1_col + stat_pvalue_manual(comparacoes_col_shannon, label = "p.signif", tip.length = 0.01)
}

# Gráfico Pielou
p2_col <- ggboxplot(metadados.dieta.alpha, x = "tercil_colesterol", y = "pielou_evenness",
                    fill = "tercil_colesterol", palette = "viridis") +
  labs(title = paste0("B. Pielou Evenness (Kruskal-Wallis p = ", signif(kw_colesterol_pielou, 3), ")"),
       x = "Cholesterol Intake (tercile)", y = "Pielou Index") +
  theme_minimal()
if (!is.null(comparacoes_col_pielou)) {
  p2_col <- p2_col + stat_pvalue_manual(comparacoes_col_pielou, label = "p.signif", tip.length = 0.01)
}

# Gráfico Chao1
p3_col <- ggboxplot(metadados.dieta.alpha, x = "tercil_colesterol", y = "chao1",
                    fill = "tercil_colesterol", palette = "viridis") +
  labs(title = paste0("C. Chao1 Richness (Kruskal-Wallis p = ", signif(kw_colesterol_chao1, 3), ")"),
       x = "Cholesterol Intake (tercile)", y = "Chao1 Richness") +
  theme_minimal()
if (!is.null(comparacoes_col_chao1)) {
  p3_col <- p3_col + stat_pvalue_manual(comparacoes_col_chao1, label = "p.signif", tip.length = 0.01)
}

# Gráfico Faith's PD
p4_col <- ggboxplot(metadados.dieta.alpha, x = "tercil_colesterol", y = "faith_pd",
                    fill = "tercil_colesterol", palette = "viridis") +
  labs(title = paste0("D. Faith's PD (Kruskal-Wallis p = ", signif(kw_colesterol_faith, 3), ")"),
       x = "Cholesterol Intake (tercile)", y = "Faith's Phylogenetic Diversity") +
  theme_minimal()
if (!is.null(comparacoes_col_faith)) {
  p4_col <- p4_col + stat_pvalue_manual(comparacoes_col_faith, label = "p.signif", tip.length = 0.01)
}

# Juntar os gráficos
painel_colesterol_final <- ggarrange(p1_col, p2_col, p3_col, p4_col, 
                                     ncol = 2, nrow = 2, 
                                     common.legend = TRUE, legend = "bottom")

# Salvar
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/painel_colesterol_significativo.png", painel_colesterol_final, width = 12, height = 8, dpi = 300)



# Shannon vs BHEI
p1_BHEI <- ggplot(metadados.dieta.alpha, aes(x = BHEI_R_Score_Total, y = shannon_entropy)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "spearman") +
  labs(title = "Shannon Entropy vs BHEI Total", x = "BHEI Total", y = "Shannon Entropy") +
  theme_minimal()

# Chao1 vs BHEI
p2_BHEI <- ggplot(metadados.dieta.alpha, aes(x = BHEI_R_Score_Total, y = chao1)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "spearman") +
  labs(title = "Chao1 Richness vs BHEI Total", x = "BHEI Total", y = "Chao1 Richness") +
  theme_minimal()


# Juntar os gráficos
scatterplot_BHEI_total <- ggarrange(p1_BHEI, p2_BHEI, 
                                     ncol = 2, nrow = 1, 
                                     common.legend = TRUE, legend = "bottom")

# Salvar
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/scatterplot_BHEI_total.png", scatterplot_BHEI_total, width = 24, height = 8, dpi = 300)

#==========================================================#
#                   Modelos Lineares Multiplos
#=========================================================#


#====== CHAO1 ===========#

modelo_chao1 <- lm(chao1 ~ BHEI_R_Score_Total + TyG + IMC + Age + Sex + Region_type, data = metadados.all.regressão)
summary(modelo_chao1)

par(mfrow=c(2,2))
plot(modelo_chao1)


# Instale o pacote se ainda não tiver
install.packages("lmtest")  # apenas uma vez

# Carregue o pacote
library(lmtest)

# Teste de Breusch-Pagan
bptest(modelo_chao1)  #Avalias se há homocedastidade

library(MASS)

modelo_chao1_robusto <- rlm(
  chao1 ~ BHEI_R_Score_Total + TyG + IMC + Age + Sex + Region_type,
  data = metadados.all.regressão
)

summary(modelo_chao1_robusto)

library(lmtest)
coeftest(modelo_chao1_robusto)

#Modelos com seleção automática de variáveis
##Você pode aplicar o método step() no modelo completo para identificar um subconjunto ótimo de variáveis com base no critério de informação de Akaike (AIC):

modelo_completo <- lm(chao1 ~ BHEI_R_Score_Total + TyG + IMC + Age + Sex + Region_type, data = metadados.all.regressão)

modelo_stepwise <- step(modelo_completo, direction = "both")

summary(modelo_stepwise)


# Carregar pacotes necessários
library(MASS)
library(lmtest)
library(ggplot2)
library(ggfortify)

# Ajustar modelo robusto
modelo_chao1_stepwise_robusto <- rlm(chao1 ~ TyG + IMC + Age, data = metadados.all.regressão)

# Ver resumo dos coeficientes com teste z
coeftest(modelo_chao1_stepwise_robusto)

# Plotar diagnóstico dos resíduos
autoplot(modelo_chao1_stepwise_robusto, which = 1:4, ncol = 2)


#====== observed_features ===========#

modelo_observed_features <- lm(observed_features ~ BHEI_R_Score_Total + TyG + IMC + Age + Sex + Region_type, data = metadados.all.regressão)
summary(modelo_observed_features)

par(mfrow=c(2,2))
plot(modelo_observed_features)


# Carregue o pacote
library(lmtest)

# Teste de Breusch-Pagan
bptest(modelo_observed_features)  #Avalias se há homocedastidade

library(MASS)

modelo_observed_features_robusto <- rlm(
  observed_features ~ BHEI_R_Score_Total + TyG + IMC + Age + Sex + Region_type,
  data = metadados.all.regressão
)

summary(modelo_observed_features_robusto)

library(lmtest)
coeftest(modelo_observed_features_robusto)



#
###================================= Statistical Power Analisys ============================#

#----------------------#
#   Statistical Power Analisys


# Criar grupos de TyG por tercis (baixo = tercil 1, alto = tercil 3)
tercis <- quantile(metadados.all.filtrado$TyG, probs = c(1/3, 2/3), na.rm = TRUE)

metadados.all.filtrado$TyG_group <- with(metadados.all.filtrado,
                                         ifelse(TyG <= tercis[1], "Baixo",
                                                ifelse(TyG > tercis[2], "Alto", NA)))  # NA para o grupo intermediário

# Teste de Mann-Whitney
wilcox.test(shannon_entropy ~ TyG_group, data = metadados.all.filtrado)

#install.packages("pwr")

library(pwr)

# Filtrar apenas os grupos baixo e alto
dados_tyG <- subset(metadados.all.filtrado, TyG_group %in% c("Baixo", "Alto"))

# Calcular média e desvio padrão de cada grupo
m_baixo <- mean(dados_tyG$shannon_entropy[dados_tyG$TyG_group == "Baixo"], na.rm = TRUE)
m_alto <- mean(dados_tyG$shannon_entropy[dados_tyG$TyG_group == "Alto"], na.rm = TRUE)
sd_pooled <- sd(dados_tyG$shannon_entropy, na.rm = TRUE)

# Calcular tamanho do efeito
d <- abs(m_baixo - m_alto) / sd_pooled

# Verificar tamanho amostral por grupo
n <- table(dados_tyG$TyG_group)
n_baixo <- n["Baixo"]
n_alto  <- n["Alto"]

# Cálculo do poder
pwr_result <- pwr.t2n.test(n1 = n_baixo, n2 = n_alto, d = d, sig.level = 0.05)
pwr_result


#---- tyg mediana----#
#Código para criar dois grupos de TyG:

# Calcular a mediana de TyG
mediana_tyg <- median(metadados.all.filtrado$TyG, na.rm = TRUE)

# Criar variável de grupo binário
metadados.all.filtrado$TyG_group2 <- with(metadados.all.filtrado,
                                          ifelse(TyG <= mediana_tyg, "Baixo", "Alto"))


# Remove possíveis NAs e mantém só os dois grupos de interesse
dados_tyG_bin <- subset(metadados.all.filtrado, TyG_group2 %in% c("Baixo", "Alto"))

# Testa se há diferença significativa no índice de Shannon entre os grupos
wilcox.test(shannon_entropy ~ TyG_group2, data = dados_tyG_bin)

#====== Cálculo do tamanho do efeito (Cohen's d) e poder=====#

# Carrega o pacote necessário para análise de poder
#install.packages("pwr")  # Execute apenas se ainda não estiver instalado
library(pwr)

# Calcula a média de Shannon em cada grupo
media_baixo <- mean(dados_tyG_bin$shannon_entropy[dados_tyG_bin$TyG_group2 == "Baixo"], na.rm = TRUE)
media_alto  <- mean(dados_tyG_bin$shannon_entropy[dados_tyG_bin$TyG_group2 == "Alto"], na.rm = TRUE)

# Calcula o desvio padrão combinado (pooled)
sd_pooled <- sd(dados_tyG_bin$shannon_entropy, na.rm = TRUE)

# Calcula o tamanho do efeito de Cohen (d)
d <- abs(media_baixo - media_alto) / sd_pooled

# Conta o número de indivíduos em cada grupo
n_baixo <- sum(dados_tyG_bin$TyG_group2 == "Baixo")
n_alto  <- sum(dados_tyG_bin$TyG_group2 == "Alto")

# Estima o poder do teste t com esses valores
pwr.t2n.test(n1 = n_baixo, n2 = n_alto, d = d, sig.level = 0.05)


#============= Poder estatístico para Chao1 ==============#

# Teste não paramétrico para comparar chao1 entre grupos de TyG
wilcox.test(chao1 ~ TyG_group2, data = dados_tyG_bin)

# Média de Chao1 nos dois grupos
media_chao_baixo <- mean(dados_tyG_bin$chao1[dados_tyG_bin$TyG_group2 == "Baixo"], na.rm = TRUE)
media_chao_alto  <- mean(dados_tyG_bin$chao1[dados_tyG_bin$TyG_group2 == "Alto"], na.rm = TRUE)

# Desvio padrão combinado
sd_chao <- sd(dados_tyG_bin$chao1, na.rm = TRUE)

# Tamanho do efeito
d_chao <- abs(media_chao_baixo - media_chao_alto) / sd_chao

# Poder estatístico para Chao1
pwr.t2n.test(n1 = n_baixo, n2 = n_alto, d = d_chao, sig.level = 0.05)


#====== Poder estatístico para Observed =========#

#Teste de Mann–Whitney para Observed Features
# Teste não paramétrico para comparar número de ASVs observados
wilcox.test(observed_features ~ TyG_group2, data = dados_tyG_bin)

# Média de observed features nos dois grupos
media_obs_baixo <- mean(dados_tyG_bin$observed_features[dados_tyG_bin$TyG_group2 == "Baixo"], na.rm = TRUE)
media_obs_alto  <- mean(dados_tyG_bin$observed_features[dados_tyG_bin$TyG_group2 == "Alto"], na.rm = TRUE)

# Desvio padrão combinado
sd_obs <- sd(dados_tyG_bin$observed_features, na.rm = TRUE)

# Tamanho do efeito
d_obs <- abs(media_obs_baixo - media_obs_alto) / sd_obs

# Poder estatístico para Observed Features
pwr.t2n.test(n1 = n_baixo, n2 = n_alto, d = d_obs, sig.level = 0.05)



#======= AGORA COM TERCIS =============#

# Calcula os tercis de TyG
tercis <- quantile(metadados.all.filtrado$TyG, probs = c(1/3, 2/3), na.rm = TRUE)

# Cria uma coluna para identificar só os extremos
metadados.all.filtrado$TyG_tercil_extremos <- with(metadados.all.filtrado,
                                                   ifelse(TyG <= tercis[1], "Baixo",
                                                          ifelse(TyG > tercis[2], "Alto", NA)))  # tercil do meio vira NA

#Filtrar os dados com apenas os grupos extremos:
dados_extremos <- subset(metadados.all.filtrado, TyG_tercil_extremos %in% c("Baixo", "Alto"))


#Teste de Mann–Whitney para Chao1 nos grupos extremos:
wilcox.test(chao1 ~ TyG_tercil_extremos, data = dados_extremos)

#Calcular Cohen's d e poder para Chao1:
# Médias de Chao1
media_chao_b <- mean(dados_extremos$chao1[dados_extremos$TyG_tercil_extremos == "Baixo"], na.rm = TRUE)
media_chao_a <- mean(dados_extremos$chao1[dados_extremos$TyG_tercil_extremos == "Alto"], na.rm = TRUE)

# Desvio padrão combinado
sd_chao_ext <- sd(dados_extremos$chao1, na.rm = TRUE)

# Tamanho do efeito
d_chao_ext <- abs(media_chao_b - media_chao_a) / sd_chao_ext

# Tamanhos dos grupos
n_b <- sum(dados_extremos$TyG_tercil_extremos == "Baixo")
n_a <- sum(dados_extremos$TyG_tercil_extremos == "Alto")

# Cálculo do poder
library(pwr)
pwr.t2n.test(n1 = n_b, n2 = n_a, d = d_chao_ext, sig.level = 0.05)

#Repetir os mesmos passos para Observed Features
# Teste
wilcox.test(observed_features ~ TyG_tercil_extremos, data = dados_extremos)

# Médias
media_obs_b <- mean(dados_extremos$observed_features[dados_extremos$TyG_tercil_extremos == "Baixo"], na.rm = TRUE)
media_obs_a <- mean(dados_extremos$observed_features[dados_extremos$TyG_tercil_extremos == "Alto"], na.rm = TRUE)

# Desvio padrão
sd_obs_ext <- sd(dados_extremos$observed_features, na.rm = TRUE)

# Tamanho do efeito
d_obs_ext <- abs(media_obs_b - media_obs_a) / sd_obs_ext

# Poder
pwr.t2n.test(n1 = n_b, n2 = n_a, d = d_obs_ext, sig.level = 0.05)

