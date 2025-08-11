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


unweighted_unifrac.raw.BrasilxSuica <- as.matrix(read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/core-metrics-results_BrasilxSuica/unweighted_unifrac_distance_matrix.qza")$data)

weighted_unifrac.raw.BrasilxSuica <- as.matrix(read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/core-metrics-results_BrasilxSuica/weighted_unifrac_distance_matrix.qza")$data)



str(unweighted_unifrac.raw.BrasilxSuica)
class(unweighted_unifrac.raw.BrasilxSuica)


# Converte unifract objects 'dist' em um data.frame 
# remove controls
# filter out glycerol samples

#unweighted
unweighted_unifrac.raw.BrasilxSuica <- data.frame(as.matrix(unweighted_unifrac.raw.BrasilxSuica))

rownames(unweighted_unifrac.raw.BrasilxSuica) <- colnames(unweighted_unifrac.raw.BrasilxSuica)

# Substituir "-" inicial por "." nos nomes das colunas
colnames(unweighted_unifrac.raw.BrasilxSuica) <- gsub("-", ".", colnames(unweighted_unifrac.raw.BrasilxSuica))

# Substituir "-" inicial por "." nos nomes das linhas
rownames(unweighted_unifrac.raw.BrasilxSuica) <- gsub("-", ".", rownames(unweighted_unifrac.raw.BrasilxSuica))

# Corrigir nomes de colunas
colnames(unweighted_unifrac.raw.BrasilxSuica) <- colnames(unweighted_unifrac.raw.BrasilxSuica) %>%
  # Remover tudo depois de ".F00", ".T1" ou "_F00", "_T1" etc
  gsub("([SB]\\d+)[._]?(F00|T1)[^\\d]*.*", "\\1.\\2", .) %>%
  # Também remover possíveis terminações adicionais como _L001, _RunX etc
  gsub("([SB]\\d+)[._]?(F00|T1).*", "\\1.\\2", .)

# Corrigir nomes de colunas
rownames(unweighted_unifrac.raw.BrasilxSuica) <- rownames(unweighted_unifrac.raw.BrasilxSuica) %>%
  # Remover tudo depois de ".F00", ".T1" ou "_F00", "_T1" etc
  gsub("([SB]\\d+)[._]?(F00|T1)[^\\d]*.*", "\\1.\\2", .) %>%
  # Também remover possíveis terminações adicionais como _L001, _RunX etc
  gsub("([SB]\\d+)[._]?(F00|T1).*", "\\1.\\2", .)


#=== Limpar weighted.unifrac

# Converte unifract objects 'dist' em um data.frame 
# remove controls
# filter out glycerol samples

#weighted
weighted_unifrac.raw.BrasilxSuica <- data.frame(as.matrix(weighted_unifrac.raw.BrasilxSuica))

rownames(weighted_unifrac.raw.BrasilxSuica) <- colnames(weighted_unifrac.raw.BrasilxSuica)

# Substituir "-" inicial por "." nos nomes das colunas
colnames(weighted_unifrac.raw.BrasilxSuica) <- gsub("-", ".", colnames(weighted_unifrac.raw.BrasilxSuica))

# Substituir "-" inicial por "." nos nomes das linhas
rownames(weighted_unifrac.raw.BrasilxSuica) <- gsub("-", ".", rownames(weighted_unifrac.raw.BrasilxSuica))



# Corrigir nomes de colunas
colnames(weighted_unifrac.raw.BrasilxSuica) <- colnames(weighted_unifrac.raw.BrasilxSuica) %>%
  # Remover tudo depois de ".F00", ".T1" ou "_F00", "_T1" etc
  gsub("([SB]\\d+)[._]?(F00|T1)[^\\d]*.*", "\\1.\\2", .) %>%
  # Também remover possíveis terminações adicionais como _L001, _RunX etc
  gsub("([SB]\\d+)[._]?(F00|T1).*", "\\1.\\2", .)

# Corrigir nomes de colunas
rownames(weighted_unifrac.raw.BrasilxSuica) <- rownames(weighted_unifrac.raw.BrasilxSuica) %>%
  # Remover tudo depois de ".F00", ".T1" ou "_F00", "_T1" etc
  gsub("([SB]\\d+)[._]?(F00|T1)[^\\d]*.*", "\\1.\\2", .) %>%
  # Também remover possíveis terminações adicionais como _L001, _RunX etc
  gsub("([SB]\\d+)[._]?(F00|T1).*", "\\1.\\2", .)

# Filtrar colunas que terminam com ".F00"
keep_cols <- grepl("\\.F00$", colnames(weighted_unifrac.raw.BrasilxSuica))
weighted_unifrac.raw.BrasilxSuica <- weighted_unifrac.raw.BrasilxSuica[, keep_cols]

# Filtrar linhas que terminam com ".F00"
keep_rows <- grepl("\\.F00$", rownames(weighted_unifrac.raw.BrasilxSuica))
weighted_unifrac.raw.BrasilxSuica <- weighted_unifrac.raw.BrasilxSuica[keep_rows, ]


# -----------------------------
# PASSO A PASSO - PCoA Unweighted UniFrac
# -----------------------------

# 1. Instalar (se necessário) e carregar o pacote
# install.packages("ape")  # Descomente esta linha se ainda não tiver o pacote
library(ape)

# 2. Verifique se a matriz está no formato adequado
# Ela deve ser simétrica e com nomes de linhas e colunas iguais (amostras)
# Aqui assumimos que unweighted_unifrac.raw.BrasilxSuica está pronta e já foi filtrada

# 3. Converter para matriz de distâncias
# Isso é necessário para o pcoa() funcionar corretamente
dist_matrix.unw <- as.dist(unweighted_unifrac.raw.BrasilxSuica)

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
    title = "PCoA - Unweighted UniFrac - BrasilxSuica",
    x = paste0("PCoA1 (", round(variance_explained.unw[1], 2), "%)"),
    y = paste0("PCoA2 (", round(variance_explained.unw[2], 2), "%)")
  ) +
  theme_minimal(base_size = 14)  # aumenta o tamanho geral das fontes do tema

# Salvar com tamanho mais equilibrado
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Pcoa_unweighted_BrasilxSuica.png",
       plot = uw1, width = 12, height = 8, dpi = 300)


#Destacar amostras contaminadas

library(ggplot2)
library(dplyr)

# Lista de amostras para destacar
amostras_destacadas <- c(
  "S10011.F00", "B10011.F00", "B30021.F00", "S30021.F00", "S10111.F00", "B10111.F00",  "S40601.F00", "B40601.F00", "B40401.F00", "S40401.F00", "S10061.F00", "B10061.F00", "S10201.F00", "B10201.F00" , "S30091.F00", "B30091.F00",  "B40241.F00", "S40241.F00", "B40651.F00", "S40651.F00", "B40481.F00", "S40481.F00",
  "B40371.F00", "S40371.F00"
)

# Criar colunas auxiliares para destaque
df_pcoa_unw <- df_pcoa_unw %>%
  mutate(
    Destacar = ifelse(SampleID %in% amostras_destacadas, "Sim", "Nao"),
    Cor = case_when(
      SampleID %in% amostras_destacadas & grepl("^B", SampleID) ~ "verde",
      SampleID %in% amostras_destacadas & grepl("^S", SampleID) ~ "vermelho",
      TRUE ~ "cinza"
    )
  )

# Gráfico
uw1_destaque <- ggplot(df_pcoa_unw, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = Cor), size = 3, alpha = 0.8) +
  scale_color_manual(
    values = c("verde" = "forestgreen", "vermelho" = "firebrick", "cinza" = "gray70"),
    guide = "none"
  ) +
  geom_text(
    data = filter(df_pcoa_unw, Destacar == "Sim"),
    aes(label = SampleID),
    vjust = -0.7, size = 4
  ) +
  labs(
    title = "PCoA - Unweighted UniFrac - BrasilxSuica",
    x = paste0("PCoA1 (", round(variance_explained.unw[1], 2), "%)"),
    y = paste0("PCoA2 (", round(variance_explained.unw[2], 2), "%)")
  ) +
  theme_minimal(base_size = 14)

#====> Amostras contaminadas aos pares

library(ggplot2)
library(dplyr)

# Definir pares como grupos para colorir e conectar
pares <- data.frame(
  SampleID = c(
    "S10011.F00", "B10011.F00", "B30021.F00", "S30021.F00", "S10111.F00", "B10111.F00",
    "S40601.F00", "B40601.F00", "B40401.F00", "S40401.F00", "S10061.F00", "B10061.F00",
    "S10201.F00", "B10201.F00", "S30091.F00", "B30091.F00", "B40241.F00", "S40241.F00",
    "B40651.F00", "S40651.F00", "B40481.F00", "S40481.F00", "B40371.F00", "S40371.F00"
  ),
  PairGroup = rep(paste0("Par", 1:12), each = 2)
)

# Juntar ao dataframe de PCoA
df_pcoa_unw_ligado <- df_pcoa_unw %>%
  left_join(pares, by = "SampleID")

# Gráfico
ggplot(df_pcoa_unw_ligado, aes(x = PCoA1, y = PCoA2)) +
  geom_point(data = df_pcoa_unw_ligado %>% filter(is.na(PairGroup)), color = "gray", size = 2, alpha = 0.6) +
  geom_point(data = df_pcoa_unw_ligado %>% filter(!is.na(PairGroup)), aes(color = PairGroup), size = 3) +
  geom_text(data = df_pcoa_unw_ligado %>% filter(!is.na(PairGroup)),
            aes(label = SampleID, color = PairGroup), vjust = -0.8, size = 3.5) +
  geom_line(data = df_pcoa_unw_ligado %>% filter(!is.na(PairGroup)) %>% group_by(PairGroup) %>% filter(n() == 2),
            aes(group = PairGroup, color = PairGroup), linewidth = 1) +
  labs(
    title = "PCoA - Unweighted UniFrac - BrasilxSuica",
    x = paste0("PCoA1 (", round(variance_explained.unw[1], 2), "%)"),
    y = paste0("PCoA2 (", round(variance_explained.unw[2], 2), "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


# -----------------------------
# PASSO A PASSO - PCoA weighted UniFrac
# -----------------------------

# 1. Instalar (se necessário) e carregar o pacote
# install.packages("ape")  # Descomente esta linha se ainda não tiver o pacote
library(ape)

# 2. Verifique se a matriz está no formato adequado
# Ela deve ser simétrica e com nomes de linhas e colunas iguais (amostras)
# Aqui assumimos que unweighted_unifrac.raw.BrasilxSuica está pronta e já foi filtrada

# 3. Converter para matriz de distâncias
# Isso é necessário para o pcoa() funcionar corretamente
dist_matrix.wei <- as.dist(weighted_unifrac.raw.BrasilxSuica)

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
    title = "PCoA - weighted UniFrac - BrasilxSuica",
    x = paste0("PCoA1 (", round(variance_explained.wei[1], 2), "%)"),
    y = paste0("PCoA2 (", round(variance_explained.wei[2], 2), "%)")
  ) +
  theme_minimal(base_size = 14)  # aumenta o tamanho geral das fontes do tema

# Salvar com tamanho mais equilibrado
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Pcoa_weighted_BrasilxSuica.png",
       plot = w1, width = 12, height = 8, dpi = 300)


#Amostra destacada

library(ggplot2)
library(dplyr)

# Lista de amostras a destacar
amostras_destacadas <- c(
  "S40601.F00", "B40601.F00", "B40401.F00", "S40401.F00",
  "B40651.F00", "S40651.F00", "B40481.F00", "S40481.F00",
  "B40371.F00", "S40371.F00"
)

amostras_destacadas <- c(
  "S10011.F00", "B10011.F00", "B30021.F00", "S30021.F00", "S10111.F00", "B10111.F00",  "S40601.F00", "B40601.F00", "B40401.F00", "S40401.F00", "S10061.F00", "B10061.F00", "S10201.F00", "B10201.F00" , "S30091.F00", "B30091.F00",  "B40241.F00", "S40241.F00", "B40651.F00", "S40651.F00", "B40481.F00", "S40481.F00",
  "B40371.F00", "S40371.F00"
)

# Criar colunas auxiliares
df_pcoa_wei <- df_pcoa_wei %>%
  mutate(
    Destacar = ifelse(SampleID %in% amostras_destacadas, "Sim", "Nao"),
    Cor = case_when(
      SampleID %in% amostras_destacadas & grepl("^B", SampleID) ~ "verde",
      SampleID %in% amostras_destacadas & grepl("^S", SampleID) ~ "vermelho",
      TRUE ~ "cinza"
    )
  )

# Gráfico
w1_destacada2 <- ggplot(df_pcoa_wei, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = Cor), size = 3, alpha = 0.8) +
  scale_color_manual(
    values = c("verde" = "forestgreen", "vermelho" = "firebrick", "cinza" = "gray70"),
    guide = "none"
  ) +
  geom_text(
    data = filter(df_pcoa_wei, Destacar == "Sim"),
    aes(label = SampleID),
    vjust = -0.7, size = 4
  ) +
  labs(
    title = "PCoA - Weighted UniFrac - BrasilxSuica",
    x = paste0("PCoA1 (", round(variance_explained.wei[1], 2), "%)"),
    y = paste0("PCoA2 (", round(variance_explained.wei[2], 2), "%)")
  ) +
  theme_minimal(base_size = 14)



#====> Amostras contaminadas aos pares Wei

library(ggplot2)
library(dplyr)

# Definir pares como grupos para colorir e conectar
pares <- data.frame(
  SampleID = c(
    "S10011.F00", "B10011.F00", "B30021.F00", "S30021.F00", "S10111.F00", "B10111.F00",
    "S40601.F00", "B40601.F00", "B40401.F00", "S40401.F00", "S10061.F00", "B10061.F00",
    "S10201.F00", "B10201.F00", "S30091.F00", "B30091.F00", "B40241.F00", "S40241.F00",
    "B40651.F00", "S40651.F00", "B40481.F00", "S40481.F00", "B40371.F00", "S40371.F00"
  ),
  PairGroup = rep(paste0("Par", 1:12), each = 2)
)

# Juntar ao dataframe de PCoA
df_pcoa_wei_ligado <- df_pcoa_wei %>%
  left_join(pares, by = "SampleID")

# Gráfico
ggplot(df_pcoa_wei_ligado, aes(x = PCoA1, y = PCoA2)) +
  geom_point(data = df_pcoa_wei_ligado %>% filter(is.na(PairGroup)), color = "gray", size = 2, alpha = 0.6) +
  geom_point(data = df_pcoa_wei_ligado %>% filter(!is.na(PairGroup)), aes(color = PairGroup), size = 3) +
  geom_text(data = df_pcoa_wei_ligado %>% filter(!is.na(PairGroup)),
            aes(label = SampleID, color = PairGroup), vjust = -0.8, size = 3.5) +
  geom_line(data = df_pcoa_wei_ligado %>% filter(!is.na(PairGroup)) %>% group_by(PairGroup) %>% filter(n() == 2),
            aes(group = PairGroup, color = PairGroup), linewidth = 1) +
  labs(
    title = "PCoA - weighted UniFrac - BrasilxSuica",
    x = paste0("PCoA1 (", round(variance_explained.wei[1], 2), "%)"),
    y = paste0("PCoA2 (", round(variance_explained.wei[2], 2), "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


#====== Para Heatmap ===========#

# Leitura correta do arquivo com ASVs na primeira coluna
SVs_rel <- read.table(
  "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/core-metrics-results_BrasilxSuica/rel_freq_exported/feature_table_rel.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  skip = 1,               # <--- pula a linha com '# Constructed from biom file'
  comment.char = "",
  check.names = FALSE
)



# Padronizar nomes
colnames(SVs_rel) <- gsub("-", ".", colnames(SVs_rel))
colnames(SVs_rel) <- gsub("([SB]\\d+)[._]?(F00|T1)[^\\d]*.*", "\\1.\\2", colnames(SVs_rel))
colnames(SVs_rel) <- gsub("([SB]\\d+)[._]?(F00|T1).*", "\\1.\\2", colnames(SVs_rel))



# Transformar em matriz numérica
SVs_rel <- as.matrix(apply(SVs_rel, 2, as.numeric))

#transpor

SVs_rel_t <- t(SVs_rel)


#=========> PCOA com simbolos




# Para UNWEIGHTED
df_pcoa_wei <- df_pcoa_wei %>%
  mutate(SamplePlate = case_when(
    SampleID %in% c(  # <- AMOSTRAS ADICIONADAS MANUALMENTE
      "S20011.F00", "S20041.F00", "S20042.F00", "S20061.F00", "S20081.F00", "S20091.F00", "S20092.F00",
      "S20101.F00", "S30011.F00", "S30021.F00", "S30031.F00", "S30032.F00", "S30041.F00", "S30042.F00",
      "S30051.F00", "S30052.F00", "S30061.F00", "S30081.F00", "S30091.F00", "S30092.F00", "S30101.F00",
      "S30112.F00", "S30121.F00", "S30122.F00", "S30131.F00", "S30171.F00", "S30181.F00", "S30211.F00",
      "S30231.F00", "S30241.F00", "S30261.F00"
    ) ~ "Sampleplate2",
    
    grepl("^S1\\d{4}", SampleID) & SampleID <= "S10171.F00" ~ "Sampleplate1",
    grepl("^S1\\d{4}", SampleID) & SampleID > "S10171.F00" & SampleID <= "S40371.F00" ~ "Sampleplate2",
    grepl("^S4\\d{4}", SampleID) ~ "Sampleplate3",
    grepl("^B", SampleID) ~ "Sampleplate4",
    
    SampleID %in% c("SNEG55_Run5_L001", "SNEG5_Run5_DL001") ~ "Sampleplate3",
    SampleID %in% c("SNEG1_Run1_DL001", "SNEG2_Run2_DL001", "SNEG3_Run3_DL001", "SNEG4_Run4_DL001", "SPOSZymo_RuL001", "S20240319_PCRNEG") ~ "Sampleplate2",
    SampleID %in% c("SPooled.ExtL001", "SPooled_Extractioncontrol_POS_MOCK", "S20240125.PL001", "S20240131_PCRNEG") ~ "Sampleplate1",
    
    TRUE ~ "Desconhecido"
  ))


library(ggplot2)
library(dplyr)

# Adiciona os pares
pares <- data.frame(
  SampleID = c(
    "S10011.F00", "B10011.F00", "B30021.F00", "S30021.F00", "S10111.F00", "B10111.F00",
    "S40601.F00", "B40601.F00", "B40401.F00", "S40401.F00", "S10061.F00", "B10061.F00",
    "S10201.F00", "B10201.F00", "S30091.F00", "B30091.F00", "B40241.F00", "S40241.F00",
    "B40651.F00", "S40651.F00", "B40481.F00", "S40481.F00", "B40371.F00", "S40371.F00"
  ),
  PairGroup = rep(paste0("Par", 1:12), each = 2)
)

# Junta com o dataframe principal
df_pcoa_wei <- df_pcoa_wei %>%
  left_join(pares, by = "SampleID")



library(ggplot2)
library(dplyr)

# Define os shapes para os sample plates
shapes <- c(
  "Sampleplate1" = 16,  # círculo
  "Sampleplate2" = 17,  # triângulo
  "Sampleplate3" = 18,  # losango
  "Sampleplate4" = 15   # quadrado
)

# Cria o gráfico
ggplot(df_pcoa_wei, aes(x = PCoA1, y = PCoA2)) +
  # Pontos cinza (sem PairGroup definido)
  geom_point(
    data = df_pcoa_wei %>% filter(is.na(PairGroup)),
    aes(shape = SamplePlate),
    color = "gray", size = 2, alpha = 0.5
  ) +
  # Pontos coloridos (com PairGroup)
  geom_point(
    data = df_pcoa_wei %>% filter(!is.na(PairGroup)),
    aes(color = PairGroup, shape = SamplePlate),
    size = 3
  ) +
  # Linhas conectando os pares
  geom_line(
    data = df_pcoa_wei %>%
      filter(!is.na(PairGroup)) %>%
      group_by(PairGroup) %>%
      filter(n() == 2),
    aes(group = PairGroup, color = PairGroup),
    linewidth = 1
  ) +
  # Nomes das amostras
  geom_text(
    data = df_pcoa_wei %>% filter(!is.na(PairGroup)),
    aes(label = SampleID, color = PairGroup),
    vjust = -0.8, size = 3.5
  ) +
  scale_shape_manual(values = shapes) +  # shapes tem que estar definido no seu ambiente
  labs(
    title = "PCoA - weighted UniFrac - BrasilxSuica",
    x = paste0("PCoA1 (", round(variance_explained.wei[1], 2), "%)"),
    y = paste0("PCoA2 (", round(variance_explained.wei[2], 2), "%)"),
    shape = "Sample Plate"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")


#======>  Weighted


# Para WEIGHTED
df_pcoa_wei <- df_pcoa_wei %>%
  mutate(SamplePlate = case_when(
    SampleID %in% c(  # <- AMOSTRAS ADICIONADAS MANUALMENTE
      "S20011.F00", "S20041.F00", "S20042.F00", "S20061.F00", "S20081.F00", "S20091.F00", "S20092.F00",
      "S20101.F00", "S30011.F00", "S30021.F00", "S30031.F00", "S30032.F00", "S30041.F00", "S30042.F00",
      "S30051.F00", "S30052.F00", "S30061.F00", "S30081.F00", "S30091.F00", "S30092.F00", "S30101.F00",
      "S30112.F00", "S30121.F00", "S30122.F00", "S30131.F00", "S30171.F00", "S30181.F00", "S30211.F00",
      "S30231.F00", "S30241.F00", "S30261.F00"
    ) ~ "Sampleplate2",
    
    grepl("^S1\\d{4}", SampleID) & SampleID <= "S10171.F00" ~ "Sampleplate1",
    grepl("^S1\\d{4}", SampleID) & SampleID > "S10171.F00" & SampleID <= "S40371.F00" ~ "Sampleplate2",
    grepl("^S4\\d{4}", SampleID) ~ "Sampleplate3",
    grepl("^B", SampleID) ~ "Sampleplate4",
    
    SampleID %in% c("SNEG55_Run5_L001", "SNEG5_Run5_DL001") ~ "Sampleplate3",
    SampleID %in% c("SNEG1_Run1_DL001", "SNEG2_Run2_DL001", "SNEG3_Run3_DL001", "SNEG4_Run4_DL001", "SPOSZymo_RuL001", "S20240319_PCRNEG") ~ "Sampleplate2",
    SampleID %in% c("SPooled.ExtL001", "SPooled_Extractioncontrol_POS_MOCK", "S20240125.PL001", "S20240131_PCRNEG") ~ "Sampleplate1",
    
    TRUE ~ "Desconhecido"
  ))


library(ggplot2)
library(dplyr)


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
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Pcoa_unweighted_BrasilxSuica_regionlabels.png",
       plot = uw1, width = 12, height = 8, dpi = 300)



#======= HEATMAP

# Lista com a ordem desejada (pares lado a lado)
ordem_amostras <- c(
  "B40371.F00", "S40371.F00",
  "B40481.F00", "S40481.F00",
  "B40651.F00", "S40651.F00",
  "B40241.F00", "S40241.F00",
  "B30091.F00", "S30091.F00",
  "B10201.F00", "S10201.F00",
  "B10061.F00", "S10061.F00",
  "B40401.F00", "S40401.F00",
  "B40601.F00", "S40601.F00",
  "B10111.F00", "S10111.F00",
  "B30021.F00", "S30021.F00",
  "B10011.F00", "S10011.F00"
)

# Verificar se todas as amostras estão na matriz
ordem_amostras <- ordem_amostras[ordem_amostras %in% colnames(SVs_rel)]

SVs_rel_ordenado <- SVs_rel[, ordem_amostras]

library(pheatmap)
library(RColorBrewer)

# Paleta com branco no centro (valores próximos de zero)
cores_personalizadas <- colorRampPalette(c("blue", "white", "red"))(100)

# Replotar o heatmap
pheatmap(
  mat = SVs_rel_filtrado,
  scale = "row",  # ou "none" se você quiser manter os valores originais
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  show_rownames = FALSE,
  main = "Heatmap - Abundância relativa (BrasilxSuiça)",
  color = cores_personalizadas
)


#========> BOXPLOT

# Transforma em long format
library(tidyverse)

# Garante que SVs_rel_ordenado está em formato data.frame
df <- as.data.frame(SVs_rel_ordenado)

# Adiciona a coluna com os nomes das ASVs
df$ASV <- rownames(df)

# Transforma em formato longo
df_long <- df %>%
  pivot_longer(
    cols = -ASV,
    names_to = "SampleID",
    values_to = "Abundancia"
  ) %>%
  mutate(
    Origem = ifelse(grepl("^B", SampleID), "Brasil", "Suíça"),
    Par = gsub("^[BS]", "", SampleID) # Ex: B40401.F00 e S40401.F00 → 40401.F00
  )

# Boxplot para comparar Brasil vs Suíça
ggplot(df_long, aes(x = Origem, y = Abundancia, group = Origem)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray") +
  geom_jitter(width = 0.2, alpha = 0.4, aes(color = Par)) +
  facet_wrap(~ASV, scales = "free_y") +
  labs(
    title = "Abundância relativa por ASV nas amostras pareadas",
    x = "", y = "Abundância relativa"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

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
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Pcoa_Weighted_BrasilxSuica_regionlabels.png",
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

