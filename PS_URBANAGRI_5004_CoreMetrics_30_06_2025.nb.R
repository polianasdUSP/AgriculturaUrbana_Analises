#Calculando core metrics direto no R, depois de aplicar Decontam nas 130 amostras

# Carregar pacotes necessários
library(phyloseq)
library(readr) # Importar metadados
library(tidyverse)


#===== Importar objetos

asv_table_final_decontam <- read.table(file = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/asv_table_final_decontam.tsv", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

metadados_130_completo <- read.csv(
  file = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/metadados_130.csv",
  header = TRUE,
  check.names = FALSE
)


#=========== Calcular ACE Antes da rarefação ============#




# Ler a tabela com os nomes corretos, incluindo a primeira coluna como "Feature.ID"
dados <- read.delim(
  "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/qiime_decontam_jun2025/feature-table.tsv",
  header = TRUE,
  skip = 1,
  row.names = 1,
  check.names = FALSE
)

# Transpor a matriz para colocar ASVs como colunas
dados_t <- t(dados)

dados_t_filtrado <- dados_t[, colSums(dados_t) > 0]

library(vegan)

resumo_ace <- estimateR(dados_t_filtrado)
str(resumo_ace)

ace_values <- resumo_ace["ACE", ]
head(ace_values)


# Criar dataframe com ACE por amostra
ace_df <- data.frame(
  IdVoluntario = colnames(resumo_ace),
  ACE = as.numeric(resumo_ace["S.ACE", ])
)

# Visualizar
head(ace_df)




