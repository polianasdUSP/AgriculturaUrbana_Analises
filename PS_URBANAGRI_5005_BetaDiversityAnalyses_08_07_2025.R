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




#Limpeza de dados: Distancias


unweighted_unifrac.raw.decontam <- as.matrix(read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/qiime_novo_28_05_2025/core-metrics-results/unweighted_unifrac_distance_matrix.qza")$data)

weighted_unifrac.raw.decontam <- as.matrix(read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Qiime_Decontam/core-metrics-results_Decontam/weighted_unifrac_distance_matrix.qza")$data)



str(unweighted_unifrac.raw.decontam)
class(unweighted_unifrac.raw.decontam)


# Converte unifract objects 'dist' em um data.frame 
# remove controls
# filter out glycerol samples

#unweighted
unweighted_unifrac.raw.decontam <- data.frame(as.matrix(unweighted_unifrac.raw.decontam))

rownames(unweighted_unifrac.raw.decontam) <- colnames(unweighted_unifrac.raw.decontam)

# Substituir "X" inicial por "S" nos nomes das colunas
colnames(unweighted_unifrac.raw.decontam) <- gsub("^X", "S", colnames(unweighted_unifrac.raw.decontam))

# Substituir "X" inicial por "S" nos nomes das linhas
rownames(unweighted_unifrac.raw.decontam) <- gsub("^X", "S", rownames(unweighted_unifrac.raw.decontam))

# Remover tudo a partir de "_DNA" nos colnames
colnames(unweighted_unifrac.raw.decontam) <- gsub("_DNA.*", "", colnames(unweighted_unifrac.raw.decontam))

# Remover tudo a partir de "_DNA" nos rownames
rownames(unweighted_unifrac.raw.decontam) <- gsub("_DNA.*", "", rownames(unweighted_unifrac.raw.decontam))

# Substituir "_" por "." nos nomes das colunas
colnames(unweighted_unifrac.raw.decontam) <- gsub("_", ".", colnames(unweighted_unifrac.raw.decontam))

# Substituir "_" por "." nos nomes das linhas
rownames(unweighted_unifrac.raw.decontam) <- gsub("_", ".", rownames(unweighted_unifrac.raw.decontam))


# Filtrar colunas que terminam com ".F00"
keep_cols <- grepl("\\.F00$", colnames(unweighted_unifrac.raw.decontam))
unweighted_unifrac.raw.decontam <- unweighted_unifrac.raw.decontam[, keep_cols]

# Filtrar linhas que terminam com ".F00"
keep_rows <- grepl("\\.F00$", rownames(unweighted_unifrac.raw.decontam))
unweighted_unifrac.raw.decontam <- unweighted_unifrac.raw.decontam[keep_rows, ]


#=== Limpar weighted.unifrac

# Converte unifract objects 'dist' em um data.frame 
# remove controls
# filter out glycerol samples

#weighted
weighted_unifrac.raw.decontam <- data.frame(as.matrix(weighted_unifrac.raw.decontam))

rownames(weighted_unifrac.raw.decontam) <- colnames(weighted_unifrac.raw.decontam)

# Substituir "X" inicial por "S" nos nomes das colunas
colnames(weighted_unifrac.raw.decontam) <- gsub("^X", "S", colnames(weighted_unifrac.raw.decontam))

# Substituir "X" inicial por "S" nos nomes das linhas
rownames(weighted_unifrac.raw.decontam) <- gsub("^X", "S", rownames(weighted_unifrac.raw.decontam))

# Remover tudo a partir de "_DNA" nos colnames
colnames(weighted_unifrac.raw.decontam) <- gsub("_DNA.*", "", colnames(weighted_unifrac.raw.decontam))

# Remover tudo a partir de "_DNA" nos rownames
rownames(weighted_unifrac.raw.decontam) <- gsub("_DNA.*", "", rownames(weighted_unifrac.raw.decontam))

# Substituir "_" por "." nos nomes das colunas
colnames(weighted_unifrac.raw.decontam) <- gsub("_", ".", colnames(weighted_unifrac.raw.decontam))

# Substituir "_" por "." nos nomes das linhas
rownames(weighted_unifrac.raw.decontam) <- gsub("_", ".", rownames(weighted_unifrac.raw.decontam))


# Filtrar colunas que terminam com ".F00"
keep_cols <- grepl("\\.F00$", colnames(weighted_unifrac.raw.decontam))
weighted_unifrac.raw.decontam <- weighted_unifrac.raw.decontam[, keep_cols]

# Filtrar linhas que terminam com ".F00"
keep_rows <- grepl("\\.F00$", rownames(weighted_unifrac.raw.decontam))
weighted_unifrac.raw.decontam <- weighted_unifrac.raw.decontam[keep_rows, ]


# -----------------------------
# PASSO A PASSO - PCoA Unweighted UniFrac
# -----------------------------

# 1. Instalar (se necessário) e carregar o pacote
# install.packages("ape")  # Descomente esta linha se ainda não tiver o pacote
library(ape)

# 2. Verifique se a matriz está no formato adequado
# Ela deve ser simétrica e com nomes de linhas e colunas iguais (amostras)
# Aqui assumimos que unweighted_unifrac.raw.decontam está pronta e já foi filtrada

# 3. Converter para matriz de distâncias
# Isso é necessário para o pcoa() funcionar corretamente
dist_matrix.unw <- as.dist(unweighted_unifrac.raw.decontam)

# 4. Calcular a PCoA
# O resultado terá os eixos principais (componentes) e os percentuais explicados
pcoa_result.unw <- pcoa(dist_matrix.unw)

# 5. Visualizar a porcentagem de variância explicada pelos dois primeiros eixos
# Isso ajuda a entender quanto da diversidade entre amostras está sendo representada
variance_explained.unw <- round(pcoa_result.unw$values$Relative_eig[1:2] * 100, 2)
cat("Variância explicada:\n")
cat("PCoA1:", variance_explained[1], "%\n")
cat("PCoA2:", variance_explained[2], "%\n")

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

# 1. Instalar (se necessário) e carregar o pacote
# install.packages("ape")  # Descomente esta linha se ainda não tiver o pacote
library(ape)

# 2. Verifique se a matriz está no formato adequado
# Ela deve ser simétrica e com nomes de linhas e colunas iguais (amostras)
# Aqui assumimos que unweighted_unifrac.raw.decontam está pronta e já foi filtrada

# 3. Converter para matriz de distâncias
# Isso é necessário para o pcoa() funcionar corretamente
dist_matrix.wei <- as.dist(weighted_unifrac.raw.decontam)

# 4. Calcular a PCoA
# O resultado terá os eixos principais (componentes) e os percentuais explicados
pcoa_result.wei <- pcoa(dist_matrix.wei)

# 5. Visualizar a porcentagem de variância explicada pelos dois primeiros eixos
# Isso ajuda a entender quanto da diversidade entre amostras está sendo representada
variance_explained.wei <- round(pcoa_result.wei$values$Relative_eig[1:2] * 100, 2)
cat("Variância explicada:\n")
cat("PCoA1:", variance_explained[1], "%\n")
cat("PCoA2:", variance_explained[2], "%\n")

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



#==== Juntar com metadados


#=======Unweighted

# Carregar pacotes necessários
library(ggplot2)
library(ggrepel)

# 1. Ler metadados
metadados <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados.raw.csv")

# 2. Juntar metadados com df_pcoa_unw
# Verificar se colunas de ID coincidem e ajustar se necessário
# Supondo que seja Sample.id no metadata e SampleID no df_pcoa_unw
df_pcoa_unw <- merge(df_pcoa_unw, metadados[, c("Sample.id", "Region")],
                     by.x = "SampleID", by.y = "Sample.id")

# 3. Calcular centro de cada região para posicionar o rótulo
centroides <- aggregate(cbind(PCoA1, PCoA2) ~ Region, data = df_pcoa_unw, mean)

# 4. Criar gráfico com pontos coloridos por região e rótulos das regiões
uw1 <- ggplot(df_pcoa_unw, aes(x = PCoA1, y = PCoA2, color = Region)) +
  geom_point(size = 3) +
  geom_text_repel(data = centroides,
                  aes(x = PCoA1, y = PCoA2, label = Region),
                  color = "black",
                  fontface = "bold",
                  size = 5,
                  box.padding = 0.5,
                  segment.color = NA) +
  labs(
    title = "PCoA - Unweighted UniFrac",
    x = paste0("PCoA1 (", round(variance_explained.unw[1], 2), "%)"),
    y = paste0("PCoA2 (", round(variance_explained.unw[2], 2), "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

# 5. Salvar o gráfico
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Pcoa_unweighted_Decontam_regionlabels.png",
       plot = uw1, width = 12, height = 8, dpi = 300)

#==== Weighted

# Carregar pacotes necessários
library(ggplot2)
library(ggrepel)

# 1. Ler metadados
metadados <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados.raw.csv")

# 2. Juntar metadados com df_pcoa_wei
# Verificar se colunas de ID coincidem e ajustar se necessário
# Supondo que seja Sample.id no metadata e SampleID no df_pcoa_wei
df_pcoa_wei <- merge(df_pcoa_wei, metadados[, c("Sample.id", "Region")],
                     by.x = "SampleID", by.y = "Sample.id")

# 3. Calcular centro de cada região para posicionar o rótulo
centroides <- aggregate(cbind(PCoA1, PCoA2) ~ Region, data = df_pcoa_wei, mean)

# 4. Criar gráfico com pontos coloridos por região e rótulos das regiões
uw1 <- ggplot(df_pcoa_wei, aes(x = PCoA1, y = PCoA2, color = Region)) +
  geom_point(size = 3) +
  geom_text_repel(data = centroides,
                  aes(x = PCoA1, y = PCoA2, label = Region),
                  color = "black",
                  fontface = "bold",
                  size = 5,
                  box.padding = 0.5,
                  segment.color = NA) +
  labs(
    title = "PCoA - Unweighted UniFrac",
    x = paste0("PCoA1 (", round(variance_explained.wei[1], 2), "%)"),
    y = paste0("PCoA2 (", round(variance_explained.wei[2], 2), "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

# 5. Salvar o gráfico
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Pcoa_Weighted_Decontam_regionlabels.png",
       plot = uw1, width = 12, height = 8, dpi = 300)

#====== Beta Diversity com Bray curtis =========#

library(vegan)

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
resultado_permanova <- adonis2(bray.curtis.limpo ~ Health_Status, data = metadados.metabolic.score)

# Mostrar resultado
print(resultado_permanova)

#=======> Jaccard com Health_Status

# Garantir os metadados com os mesmos IDs
ids_validos <- rownames(jaccard.limpo)
metadados_filtrados_jaccard <- metadados.metabolic.score[metadados.metabolic.score$Sample.id %in% ids_validos, ]

# Reordenar para garantir correspondência perfeita
metadados_filtrados_jaccard <- metadados_filtrados_jaccard[match(ids_validos, metadados_filtrados_jaccard$Sample.id), ]

# Verificar correspondência
stopifnot(identical(rownames(jaccard.limpo), metadados_filtrados_jaccard$Sample.id))



set.seed(123)  # Para reprodutibilidade
resultado_permanova_jaccard <- adonis2(jaccard.limpo ~ Health_Status, data = metadados_filtrados_jaccard)

# Mostrar resultado
print(resultado_permanova_jaccard)

