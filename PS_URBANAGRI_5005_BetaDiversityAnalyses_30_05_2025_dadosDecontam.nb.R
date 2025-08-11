#Intalando pacotes necessários



library(remotes)
library(devtools)
library(Rcpp)
library(Biostrings)
library(rhdf5)
library(phyloseq)
library(S4Vectors)
library(TreeSummarizedExperiment)

#devtools::install_github("jbisanz/qiime2R")

#install.packages(c("ggrepel", "gridExtra", "viridis", "ggpubr", "tidyverse", "vegan", "dplyr", "ggmisc", "ggplot2", "readr"))
#install.packages("ggpmisc")

library(qiime2R)
library(ggrepel)
library(gridExtra)
library(viridis)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(ggpmisc)
library(ggplot2)
library(readr)
# 1. Instalar (se necessário) e carregar o pacote
# install.packages("ape")  # Descomente esta linha se ainda não tiver o pacote
library(ape)
library(vegan)





#Limpeza de dados: Distancias


unweighted_unifrac.raw <- as.matrix(read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/qiime_decontam_jun2025/core-metrics-results/unweighted_unifrac_distance_matrix.qza")$data)

weighted_unifrac.raw <- as.matrix(read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/qiime_decontam_jun2025/core-metrics-results/weighted_unifrac_distance_matrix.qza")$data)



# Vetor com IDs a remover
ids_remover <- c("S30092.F00", "S40142.F00")

# Unweighted UniFrac
unweighted_unifrac <- unweighted_unifrac.raw[!(rownames(unweighted_unifrac.raw) %in% ids_remover),
                                                   !(colnames(unweighted_unifrac.raw) %in% ids_remover)]

# Weighted UniFrac
weighted_unifrac <- weighted_unifrac.raw[!(rownames(weighted_unifrac.raw) %in% ids_remover),
                                               !(colnames(weighted_unifrac.raw) %in% ids_remover)]


str(unweighted_unifrac)
class(unweighted_unifrac)


# Converte unifract objects 'dist' em um data.frame 
# remove controls
# filter out glycerol samples

#unweighted
unweighted_unifrac <- data.frame(as.matrix(unweighted_unifrac))

rownames(unweighted_unifrac) <- colnames(unweighted_unifrac)


#=== Limpar weighted.unifrac

# Converte unifract objects 'dist' em um data.frame 
# remove controls
# filter out glycerol samples

#weighted
weighted_unifrac <- data.frame(as.matrix(weighted_unifrac.raw))

rownames(weighted_unifrac) <- colnames(weighted_unifrac)


# -----------------------------
# PASSO A PASSO - PCoA Unweighted UniFrac
# -----------------------------



# 2. Verifique se a matriz está no formato adequado
# Ela deve ser simétrica e com nomes de linhas e colunas iguais (amostras)
# Aqui assumimos que unweighted_unifrac.raw.decontam está pronta e já foi filtrada

# 3. Converter para matriz de distâncias
# Isso é necessário para o pcoa() funcionar corretamente
dist_matrix.unw <- as.dist(unweighted_unifrac)

# 4. Calcular a PCoA
# O resultado terá os eixos principais (componentes) e os percentuais explicados
pcoa_result.unw <- pcoa(dist_matrix.unw)

# 5. Visualizar a porcentagem de variância explicada pelos dois primeiros eixos
# Isso ajuda a entender quanto da diversidade entre amostras está sendo representada
variance_explained.unw <- round(pcoa_result.unw$values$Relative_eig[1:2] * 100, 2)
cat("Variância explicada:\n")
cat("PCoA1:", variance_explained.unw[1], "%\n")
cat("PCoA2:", variance_explained.unw[2], "%\n")

# Preparar o dataframe com os dois primeiros eixos da PCoA
df_pcoa_unw <- as.data.frame(pcoa_result.unw$vectors[, 1:2])
colnames(df_pcoa_unw) <- c("PCoA1", "PCoA2")
df_pcoa_unw$SampleID <- rownames(df_pcoa_unw)

uw1 <- ggplot(df_pcoa_unw, aes(x = PCoA1, y = PCoA2, label = SampleID)) +
  geom_point(color = "blue", size = 3) +  # aumenta o tamanho do ponto
  geom_text(vjust = -0.7, size = 4) +     # aumenta o tamanho da fonte
  labs(
    title = "PCoA - Unweighted UniFrac - Decontam",
    x = paste0("PCoA1 (", round(variance_explained.unw[1], 2), "%)"),
    y = paste0("PCoA2 (", round(variance_explained.unw[2], 2), "%)")
  ) +
  theme_minimal(base_size = 14)  # aumenta o tamanho geral das fontes do tema

# Salvar com tamanho mais equilibrado
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Pcoa_unweighted_Decontam.png",
       plot = uw1, width = 12, height = 8, dpi = 300)


# -----------------------------
# PASSO A PASSO - PCoA weighted UniFrac
# -----------------------------



# 2. Verifique se a matriz está no formato adequado
# Ela deve ser simétrica e com nomes de linhas e colunas iguais (amostras)
# Aqui assumimos que unweighted_unifrac.raw.decontam está pronta e já foi filtrada

# 3. Converter para matriz de distâncias
# Isso é necessário para o pcoa() funcionar corretamente
dist_matrix.wei <- as.dist(weighted_unifrac)

# 4. Calcular a PCoA
# O resultado terá os eixos principais (componentes) e os percentuais explicados
pcoa_result.wei <- pcoa(dist_matrix.wei)

# 5. Visualizar a porcentagem de variância explicada pelos dois primeiros eixos
# Isso ajuda a entender quanto da diversidade entre amostras está sendo representada
variance_explained.wei <- round(pcoa_result.wei$values$Relative_eig[1:2] * 100, 2)
cat("Variância explicada:\n")
cat("PCoA1:", variance_explained.wei[1], "%\n")
cat("PCoA2:", variance_explained.wei[2], "%\n")

# Preparar o dataframe com os dois primeiros eixos da PCoA
df_pcoa_wei <- as.data.frame(pcoa_result.wei$vectors[, 1:2])
colnames(df_pcoa_wei) <- c("PCoA1", "PCoA2")
df_pcoa_wei$SampleID <- rownames(df_pcoa_wei)

w1 <- ggplot(df_pcoa_wei, aes(x = PCoA1, y = PCoA2, label = SampleID)) +
  geom_point(color = "blue", size = 3) +  # aumenta o tamanho do ponto
  geom_text(vjust = -0.7, size = 4) +     # aumenta o tamanho da fonte
  labs(
    title = "PCoA - weighted UniFrac - Decontam",
    x = paste0("PCoA1 (", round(variance_explained.wei[1], 2), "%)"),
    y = paste0("PCoA2 (", round(variance_explained.wei[2], 2), "%)")
  ) +
  theme_minimal(base_size = 14)  # aumenta o tamanho geral das fontes do tema

# Salvar com tamanho mais equilibrado
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Pcoa_weighted_Decontam.png",
       plot = w1, width = 12, height = 8, dpi = 300)


#=================================#
#     Juntar com metadados         #
#==================================#


#=======Unweighted



# 1. Ler metadados
#metadados <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados.raw.csv")

# 2. Juntar metadados com df_pcoa_unw
# Verificar se colunas de ID coincidem e ajustar se necessário
# Supondo que seja Sample.id no metadata e SampleID no df_pcoa_unw
df_pcoa_unw <- merge(df_pcoa_unw, metadados_grupos_saude_binario[, c("Sample.id", "Health_Status")],
                     by.x = "SampleID", by.y = "Sample.id")




# Isso ajuda a entender quanto da diversidade entre amostras está sendo representada
variance_explained.unw <- round(pcoa_result.unw$values$Relative_eig[1:2] * 100, 2)
cat("Variância explicada:\n")
cat("PCoA1:", variance_explained.unw[1], "%\n")
cat("PCoA2:", variance_explained.unw[2], "%\n")

# Rodar PERMANOVA
permanova_result_unwei <- adonis2(as.dist(unweighted_unifrac.limpo) ~ Health_Status, data = metadados_grupos_saude_binario)

io]
´~# Ver resultado
print(permanova_result_unwei)

#Extrair valor-p e R² para colocar no gráfico
p_value <- permanova_result_unwei$`Pr(>F)`[1]
r2_value <- permanova_result_unwei$R2[1]

# Criar label formatado
permanova_label_unwei <- paste0("PERMANOVA: R² = ", round(r2_value, 3), 
                          ", p = ", format.pval(p_value, digits = 3, eps = .001))




# 4. Criar gráfico com pontos coloridos por região e rótulos das regiões
uw1 <- ggplot(df_pcoa_unw, aes(x = PCoA1, y = PCoA2, color = Health_Status)) +
  geom_point(size = 3) +
  labs(
    title = "PCoA - Unweighted UniFrac",
    x = paste0("PCoA1 (", round(variance_explained.unw[1], 2), "%)"),
    y = paste0("PCoA2 (", round(variance_explained.unw[2], 2), "%)")
  ) +
  annotate("text", 
           x = 0.3, y = 0.25,   # ajuste conforme necessário
           label = permanova_label_unwei, 
           hjust = 1, size = 5) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

print(uw1)

# 5. Salvar o gráfico
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Pcoa_unweighted_Decontam_Health_Status_labels.png",
       plot = uw1, width = 12, height = 8, dpi = 300)

#==== Weighted


# 1. Ler metadados
#metadados <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados.raw.csv")

# 2. Juntar metadados com df_pcoa_wei
# Verificar se colunas de ID coincidem e ajustar se necessário
# Supondo que seja Sample.id no metadata e SampleID no df_pcoa_wei
df_pcoa_wei <- merge(df_pcoa_wei, metadados_grupos_saude_binario[, c("Sample.id", "Health_Status")],
                     by.x = "SampleID", by.y = "Sample.id")


# 5. Visualizar a porcentagem de variância explicada pelos dois primeiros eixos
# Isso ajuda a entender quanto da diversidade entre amostras está sendo representada
variance_explained.wei <- round(pcoa_result.wei$values$Relative_eig[1:2] * 100, 2)
cat("Variância explicada:\n")
cat("PCoA1:", variance_explained[1], "%\n")
cat("PCoA2:", variance_explained[2], "%\n")

# Rodar PERMANOVA
permanova_result_wei <- adonis2(as.dist(weighted_unifrac.limpo) ~ Health_Status, data = metadados_grupos_saude_binario)

# Ver resultado
print(permanova_result_wei)

#Extrair valor-p e R² para colocar no gráfico
p_value <- permanova_result_wei$`Pr(>F)`[1]
r2_value <- permanova_result_wei$R2[1]

# Criar label formatado
permanova_label <- paste0("PERMANOVA: R² = ", round(r2_value, 3), 
                          ", p = ", format.pval(p_value, digits = 3, eps = .001))




# 4. Criar gráfico com pontos coloridos por região e rótulos das regiões
w1 <- ggplot(df_pcoa_wei, aes(x = PCoA1, y = PCoA2, color = Health_Status)) +
  geom_point(size = 3) +
  labs(
    title = "PCoA - Weighted UniFrac",
    x = paste0("PCoA1 (", round(variance_explained.wei[1], 2), "%)"),
    y = paste0("PCoA2 (", round(variance_explained.wei[2], 2), "%)")
  ) +
  annotate("text", 
           x = 0.3, y = 0.25,   # ajuste conforme necessário
           label = permanova_label, 
           hjust = 1, size = 5) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

# 5. Salvar o gráfico
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Pcoa_Weighted_Decontam_Health_Statuslabels.png",
       plot = w1, width = 12, height = 8, dpi = 300)

#====== Beta Diversity com Bray curtis =========#

bray_curtis_raw <- as.matrix(read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/qiime_decontam_jun2025/core-metrics-results/bray_curtis_distance_matrix.qza")$data)

# Unweighted UniFrac
bray_curtis <- bray_curtis_raw[!(rownames(bray_curtis_raw) %in% ids_remover),
                                             !(colnames(bray_curtis_raw) %in% ids_remover)]

# Passo 1: Garantir que os metadados tenham os mesmos IDs da matriz
# (apenas os que sobraram após limpeza)
ids_validos <- rownames(bray.curtis.limpo)
metadados_filtrados <- metadados.metabolic.score[metadados.metabolic.score$Sample.id %in% ids_validos, ]

# Passo 2: Garantir que a ordem das linhas do metadata bata com a matriz
metadados_filtrados <- metadados_filtrados[match(ids_validos, metadados_filtrados$Sample.id), ]

# Verificar se está tudo certo
stopifnot(identical(rownames(bray.curtis.limpo), metadados_filtrados$Sample.id))

# Passo 3: Rodar PERMANOVA
set.seed(123)  # Para reprodutibilidade
resultado_permanova <- adonis2(bray_curtis ~ Health_Status, data = metadados_grupos_saude_binario)

# Mostrar resultado
print(resultado_permanova)

#=======> Jaccard com Health_Status

jaccard_raw <- as.matrix(read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/qiime_decontam_jun2025/core-metrics-results/jaccard_distance_matrix.qza")$data)

# Unweighted UniFrac
jaccard <- jaccard_raw[!(rownames(jaccard_raw) %in% ids_remover),
                               !(colnames(jaccard_raw) %in% ids_remover)]




set.seed(123)  # Para reprodutibilidade
resultado_permanova_jaccard <- adonis2(jaccard ~ Health_Status, data = metadados_grupos_saude_binario)

# Mostrar resultado
print(resultado_permanova_jaccard)

