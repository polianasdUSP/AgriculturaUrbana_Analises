#install.packages("BiocManager")
BiocManager::install("decontam")

metadata <- read.delim(
  "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/qiime_novo_28_05_2025/metadata.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

metadata <- read.delim("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/qiime_novo_28_05_2025/metadata.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Remover sufixos tipo "_DNA_S42_L001"
#metadata$sample.id <- gsub("_DNA_S\\d+_L001$", "", metadata$sample.id)

# Substituir underscores por traços
#metadata$sample.id <- gsub("_", "-", metadata$sample.id)

# Remover linhas com "#q2:types"
#metadata <- metadata[metadata$sample.id != "#q2:types", ]

# Remover amostras que terminam com .F01
#metadata <- metadata[!grepl("\\.F01$", metadata$sample.id), ]

# Remover controles negativos (assumindo que contenham "NEG", "neg", "Control" ou "POSZymo" no nome)
#metadata <- metadata[!grepl("NEG|neg|Control|POSZymo", metadata$sample.id), ]

# Verificar resultado
#length(unique(metadata$sample.id))
#head(metadata$sample.id)


#C:\Users\polia\OneDrive\Desktop\EstatisticaR\AgrUrbana\16S_AgriUrbana\qiime_novo_28_05_2025

library(phyloseq)
library(decontam)
library(qiime2R)
library(microbiome)
library(vegan)


#criar objeto phyloseq

  library(phyloseq)

physeq <- qza_to_phyloseq(
  features="C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/qiime_novo_28_05_2025/table.qza",
  tree="C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/qiime_novo_28_05_2025/rooted-tree.qza",
  metadata = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/qiime_novo_28_05_2025/metadata.tsv"
)





#controles_negativos <- c(
 # "Pooled-Extractioncontrol-NEG_S56_L001",
  #"Pooled-Extractioncontrol-POS-MOCK_S55_L001",
  #"20240125-PCRNEG_S57_L001",
  #"20240131-PCRNEG_S96_L001",
  #"20240319_PCRNEG_S96_L001",
  #"NEG_Run1_DNA_S5_L001",
  #"NEG_Run2_DNA_S30_L001",
  #"NEG_Run3_DNA_S54_L001",
  #"NEG_Run4_DNA_S75_L001",
  #"POSZymo_Run1_DNA_S6_L001",
  #"NEG_Kit_Run5_DNA_S4_L001",
  #"NEG_Run5_DNA_S13_L001",
  #"20240325_Neg_Control_DNA_S70_L001",
  #"20240404_PCRneg_S95_L001"
#)


#========= Decontam para os 3 sample Sheets =========#


#====> s1

# Selecionar apenas amostras de SampleSheet1 que terminam com -F00
samples_s1 <- sample_names(physeq)[grepl("^10[01]\\d{2}-F00$", sample_names(physeq))]

# Controles negativos desse batch
neg_controls_s1 <- c(
  "Pooled-Extractioncontrol-NEG_S56_L001",
  "Pooled-Extractioncontrol-POS-MOCK_S55_L001",
  "20240125-PCRNEG_S57_L001",
  "20240131-PCRNEG_S96_L001"
)

# Subset do objeto phyloseq para esse batch
ps_s1 <- prune_samples(sample_names(physeq) %in% c(samples_s1, neg_controls_s1), physeq)

# Indicar os controles negativos
sample_data(ps_s1)$is.neg <- sample_names(ps_s1) %in% neg_controls_s1


#=====> s2

# SampleSheet2: manter apenas amostras terminando com _F00
samples_s2 <- sample_names(physeq)[grepl("^(10[0-4]|20[0-1]|30[0-2]|40[0-3])\\d{2}_F00(_DNA)?", sample_names(physeq))]

# Controles negativos para SampleSheet2
neg_controls_s2 <- c(
  "20240319_PCRNEG_S96_L001",
  "NEG_Run1_DNA_S5_L001",
  "NEG_Run2_DNA_S30_L001",
  "NEG_Run3_DNA_S54_L001",
  "NEG_Run4_DNA_S75_L001",
  "POSZymo_Run1_DNA_S6_L001"
)

# Subset do phyloseq para s2
ps_s2 <- prune_samples(sample_names(physeq) %in% c(samples_s2, neg_controls_s2), physeq)
sample_data(ps_s2)$is.neg <- sample_names(ps_s2) %in% neg_controls_s2


#====> s3

# SampleSheet3: manter apenas amostras terminando com _F00
samples_s3 <- sample_names(physeq)[grepl("^40[4-8]\\d{2}_F00(_DNA)?", sample_names(physeq))]

# Controles negativos para SampleSheet3
neg_controls_s3 <- c(
  "NEG_Kit_Run5_DNA_S4_L001",
  "NEG_Run5_DNA_S13_L001",
  "20240325_Neg_Control_DNA_S70_L001",
  "20240404_PCRneg_S95_L001"
)

# Subset do phyloseq para s3
ps_s3 <- prune_samples(sample_names(physeq) %in% c(samples_s3, neg_controls_s3), physeq)
sample_data(ps_s3)$is.neg <- sample_names(ps_s3) %in% neg_controls_s3


#=============================================================#
#                Decontam S1: threshold .01
#=============================================================#
#contam_s1 <- isContaminant(ps_s1, method = "prevalence", neg = "is.neg", threshold = 0.5)  


# 4. Criar coluna lógica com os negativos
sample_data(ps_s1)$is.neg <- sample_names(ps_s1) %in% neg_controls_s1

# 5. Rodar o decontam (método prevalência)
contam_s1 <- isContaminant(ps_s1, method = "prevalence", neg = "is.neg")  

# 6. Verificar quantos contaminantes foram encontrados
table(contam_s1$contaminant)

# 7. Filtrar ASVs não contaminantes
asvs_s1 <- taxa_names(ps_s1)[!contam_s1$contaminant]

# 8. Criar novo objeto phyloseq sem os contaminantes
ps_s1_clean <- prune_taxa(asvs_s1, ps_s1)

# 9. Verificação opcional
ps_s1_clean




#=======> CONTAMINANTES E SUA TAXONOMIA:

# Taxonomia dos contaminantes
contam_asvs <- rownames(contam_s1)[contam_s1$contaminant]
contam_tax <- tax_table(ps)[contam_asvs, ]

# Em formato data.frame
contam_tax_df <- as.data.frame(contam_tax)
print(contam_tax_df)



#write.csv(contam_tax_df, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/contaminantes_sampleSheet1.csv", row.names = TRUE)





#=> Quantos reads se perderam?

# 1. Total de reads por amostra antes e depois
reads_before <- sample_sums(ps_s1)
reads_after <- sample_sums(ps_s1_clean)

# 2. Calcular reads removidos e porcentagem
reads_removed <- reads_before - reads_after
percent_removed <- (reads_removed / reads_before) * 100

# 3. Criar dataframe com os resultados
removed_df <- data.frame(
  Sample = names(reads_before),
  Reads_Before = reads_before,
  Reads_After = reads_after,
  Reads_Removed = reads_removed,
  Percent_Removed = round(percent_removed, 2)
)

# 4. Visualizar os dados
removed_df <- removed_df[order(removed_df$Percent_Removed, decreasing = TRUE), ]
print(removed_df)


# 4. Remover os controles (NEG e POS)
removed_df_real <- removed_df %>%
  filter(!grepl("NEG|POS", Sample))

# 5. Reordenar Sample com base no % removido
removed_df_real$Sample <- factor(
  removed_df_real$Sample,
  levels = removed_df_real$Sample[order(removed_df_real$Percent_Removed, decreasing = TRUE)]
)

# 6. Criar gráfico de barras
ggplot(removed_df_real, aes(x = Sample, y = Percent_Removed)) +
  geom_bar(stat = "identity", fill = "#4682B4") +
  labs(
    title = "Percentual de Reads Removidos por Amostra (S1, sem controles). Threshold .01",
    x = "Amostra",
    y = "% de Reads Removidos"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )




#=============================================================#
#                 Decontam S2: threshold .01 
#=============================================================#

library(phyloseq)
library(decontam)


# 4. Criar a coluna lógica com os negativos
sample_data(ps_s2)$is.neg <- sample_names(ps_s2) %in% neg_controls_s2

# 5. Rodar o decontam (método prevalência)
contam_s2 <- isContaminant(ps_s2, method = "prevalence", neg = "is.neg")

# 6. Verificar quantos contaminantes foram identificados
table(contam_s2$contaminant)

# 7. Obter ASVs que não são contaminantes
asvs_s2 <- taxa_names(ps_s2)[!contam_s2$contaminant]

# 8. Criar novo objeto phyloseq sem os contaminantes
ps_s2_clean <- prune_taxa(asvs_s2, ps_s2)

# 9. Confirmar visualmente
ps_s2_clean

#VERIFICAR OS CONTAMINANTES

# 1. Identificar os ASVs considerados contaminantes
contam_asvs_s2 <- rownames(contam_s2)[contam_s2$contaminant]

# 2. Ver a taxonomia dos contaminantes
contam_tax_s2 <- tax_table(ps)[contam_asvs_s2, ]
contam_tax_s2_df <- as.data.frame(contam_tax_s2)

# 3. Visualizar
print(contam_tax_s2_df)



write.csv(contam_tax_s2_df, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/contaminantes_s2.csv", row.names = TRUE)

#generos dos contaminantes: 

table(contam_tax_s2_df$Genus)


#=> Quantos reads se perderam?

# 1. Total de reads por amostra antes e depois
reads_before2 <- sample_sums(ps_s2)
reads_after2 <- sample_sums(ps_s2_clean)

# 2. Calcular reads removidos e porcentagem
reads_removed2 <- reads_before2 - reads_after2
percent_removed2 <- (reads_removed2 / reads_before2) * 100

# 3. Criar dataframe com os resultados
removed_df2 <- data.frame(
  Sample = names(reads_before2),
  Reads_Before2 = reads_before2,
  Reads_After2 = reads_after2,
  Reads_Removed2 = reads_removed2,  # << CORRIGIDO AQUI
  Percent_Removed2 = round(percent_removed2, 2)
)

# 4. Visualizar os dados ordenados por % removida
removed_df2 <- removed_df2[order(removed_df2$Percent_Removed2, decreasing = TRUE), ]
print(removed_df2)



library(ggplot2)
library(dplyr)

# 1. Filtrar amostras reais (sem NEG ou POS)
removed_df2_real <- removed_df2 %>%
  filter(!grepl("NEG|POS", Sample))

# 2. Ordenar as amostras pelo percentual removido
removed_df2_real$Sample <- factor(
  removed_df2_real$Sample,
  levels = removed_df2_real$Sample[order(removed_df2_real$Percent_Removed2, decreasing = TRUE)]
)

# 3. Criar o gráfico com amostras ordenadas
ggplot(removed_df2_real, aes(x = Sample, y = Percent_Removed2)) +
  geom_bar(stat = "identity", fill = "#4682B4") +
  labs(
    title = "Percentual de Reads Removidos por Amostra (sem controles)",
    x = "Amostra",
    y = "% de Reads Removidos"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )



#=============================================================#
#                        Decontam S3
#=============================================================#

library(phyloseq)
library(decontam)



# 3. Subset do objeto phyloseq com amostras e controles
ps_s3 <- prune_samples(sample_names(ps) %in% c(samples_s3, neg_controls_s3), physeq)

# 4. Marcar os controles negativos
sample_data(ps_s3)$is.neg <- sample_names(ps_s3) %in% neg_controls_s3


#2. O limiar (threshold) do decontam pode estar conservador
#Por padrão, o threshold = 0.1, o que pode ser pouco sensível para detectar contaminantes que também aparecem em amostras reais.
#testar uma detecção mais rigorosa, pode tentar:

# 5. Rodar o método de prevalência do decontam
contam_s3 <- isContaminant(ps_s3, method = "prevalence", neg = "is.neg")

# 6. Verificar quantos contaminantes foram detectados
table(contam_s3$contaminant)

# 7. Obter os ASVs não contaminantes
asvs_s3 <- taxa_names(ps_s3)[!contam_s3$contaminant]

# 8. Criar objeto phyloseq limpo
ps_s3_clean <- prune_taxa(asvs_s3, ps_s3)

# 9. Verificação
ps_s3_clean

#TAXONOMIA DOS CONTAMINANTES

# 1. Identificar os ASVs considerados contaminantes
contam_asvs_s3 <- rownames(contam_s3)[contam_s3$contaminant]

# 2. Obter a taxonomia desses ASVs
contam_tax_s3 <- tax_table(ps)[contam_asvs_s3, ]
contam_tax_s3_df <- as.data.frame(contam_tax_s3)

# 3. Visualizar
print(contam_tax_s3_df)






write.csv(contam_tax_s3_df, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/contaminantes_s3.csv", row.names = TRUE)


#outras asvs a verificar: f__Mycoplasmataceae (2b9c8fa41ae2b255b4f3831f752010e5)  f__Succinivibrionaceae; g__Ruminobacter; s__gut_metagenome (3a0c5873e8683c4c3cec68faa75b3bb0) 

asvs_verificar <- c("2b9c8fa41ae2b255b4f3831f752010e5", "3a0c5873e8683c4c3cec68faa75b3bb0")




#=====> verificar os controles

controles <- sample_names(ps_s3)[grepl("NEG|Neg|PCRneg|Control", sample_names(ps_s3))]
ps_s3_controles <- subset_samples(ps_s3, sample_names(ps_s3) %in% controles)

#Ver composição (ex: por Genus ou Phylum):

# Transformar em proporção relativa
ps_s3_controles_rel <- transform_sample_counts(ps_s3_controles, function(x) x / sum(x))

library(phyloseq)
library(ggplot2)
library(dplyr)

# Derreter objeto phyloseq para data.frame
df_genus <- psmelt(genus_abund)

# Selecionar os 20 gêneros mais abundantes
top_genus <- df_genus %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 20) %>%
  pull(Genus)

# Filtrar apenas os 20 gêneros mais abundantes
df_top <- df_genus %>% filter(Genus %in% top_genus)

# Gráfico
ggplot(df_top, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundância Relativa por Gênero nos Controles",
       x = "Amostra", y = "Abundância Relativa") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Extrair abundância de Treponema nos controles
df_treponema <- psmelt(ps_s3_controles_rel) %>%
  filter(Genus == "Treponema")

# Visualizar
print(df_treponema)

# Se quiser um gráfico:
ggplot(df_treponema, aes(x = Sample, y = Abundance, fill = Sample)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundância de Treponema nos Controles",
       x = "Amostra", y = "Abundância Relativa") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = "none")



#=> Quantos reads se perderam?

# 1. Total de reads por amostra antes e depois
reads_before3 <- sample_sums(ps_s3)
reads_after3 <- sample_sums(ps_s3_clean)

# 3. Calcular reads removidos e porcentagem
reads_removed3 <- reads_before3 - reads_after3
percent_removed3 <- (reads_removed3 / reads_before3) * 100

# 3. Criar dataframe com os resultados
removed_df3 <- data.frame(
  Sample = names(reads_before3),
  Reads_Before3 = reads_before3,
  Reads_After3 = reads_after3,
  Reads_Removed3 = reads_removed3,  # << CORRIGIDO AQUI
  Percent_Removed3 = round(percent_removed3, 2)
)

# 4. Visualizar os dados ordenados por % removida
removed_df3 <- removed_df3[order(removed_df3$Percent_Removed3, decreasing = TRUE), ]
print(removed_df3)


# 1. Filtrar amostras reais (excluindo controles)
removed_df3_real <- removed_df3 %>%
  filter(
    !grepl("NEG|POS", Sample) & 
      !grepl("20240325_Neg_Control", Sample)
  )

# 2. Ordenar as amostras pelo percentual removido
removed_df3_real$Sample <- factor(
  removed_df3_real$Sample,
  levels = removed_df3_real$Sample[order(removed_df3_real$Percent_Removed3, decreasing = TRUE)]
)

# 3. Criar o gráfico com amostras ordenadas
ggplot(removed_df3_real, aes(x = Sample, y = Percent_Removed3)) +
  geom_bar(stat = "identity", fill = "#4682B4") +
  labs(
    title = "Percentual de Reads Removidos por Amostra (sem controles)",
    x = "Amostra",
    y = "% de Reads Removidos"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )





#=============================================================#
#                Juntar os três objetos limpos
#=============================================================#

# Juntar os objetos phyloseq limpos
physeq_descontaminado <- merge_phyloseq(ps_s1_clean, ps_s2_clean, ps_s3_clean)

# Verificar dimensões do objeto combinado
physeq_descontaminado


#=====> rarefação

#Verificar a distribuição de reads
#Antes de definir um ponto de corte para a rarefação:
  

hist(sample_sums(physeq_descontaminado),
     breaks = 50,
     main = "Distribuição de Reads por Amostra",
     xlab = "Número de Reads",
     col = "#66c2a5")


reads_por_amostra <- sample_sums(physeq_descontaminado)
reads_por_amostra[which.min(reads_por_amostra)]

#> reads_por_amostra[which.min(reads_por_amostra)]
#30081.F00 
#17693 

#tirei a amostra 30081.F00. Rarefação será com 25689 para maior profundidade e para manter amostra 10391

set.seed(123)
physeq_rarefied <- rarefy_even_depth(
  physeq_descontaminado,
  sample.size = 25689,
  rngseed = 123,
  verbose = TRUE
)

#verificar
sample_sums(physeq_rarefied)




#=================================================================================#
#       Criar beta diversidade
#=================================================================================#

#verificar se posso usar a tree sem podar



# ASVs que estão no physeq mas não estão na árvore
setdiff(taxa_names(physeq_rarefied), tree$tip.label) # => se retornar vazio, está ok!



# Adicionar a árvore ao objeto phyloseq
phy_tree(physeq_rarefied) <- tree



# Calcular UniFrac ponderado (weighted)
unifrac_weighted <- phyloseq::UniFrac(physeq_rarefied, weighted = TRUE)

# Calcular UniFrac não ponderado (unweighted)
unifrac_unweighted <- phyloseq::UniFrac(physeq_rarefied, weighted = FALSE, normalized = TRUE, parallel = TRUE, fast = TRUE)


#Realizar pcoa
ord_weighted <- ordinate(physeq_rarefied, method = "PCoA", distance = unifrac_weighted)

#Transformar em dataframe e adicionar nomes
df_pcoa_weighted <- as.data.frame(ord_weighted$vectors)
df_pcoa_weighted$SampleID <- rownames(df_pcoa_weighted)

#plotar

library(ggplot2)

ggplot(df_pcoa_weighted, aes(x = Axis.1, y = Axis.2)) +
  geom_point() +
  labs(
    title = "PCoA - Weighted UniFrac (sem outlier)",
    x = paste0("Axis 1 [", round(ord_weighted$values$Relative_eig[1] * 100, 1), "%]"),
    y = paste0("Axis 2 [", round(ord_weighted$values$Relative_eig[2] * 100, 1), "%]")
  ) +
  theme_minimal()



#nomes nos pontos

library(ggplot2)
library(ggrepel)

g <- ggplot(df_pcoa_weighted, aes(x = Axis.1, y = Axis.2)) +
  geom_point(size = 2) +
  geom_text(aes(label = Sample), size = 3, vjust = -0.5, hjust = 0.5) +
  labs(
    title = "PCoA - Weighted UniFrac (sem outlier)",
    x = "Axis 1",
    y = "Axis 2"
  ) +
  theme_minimal()

print(g)


df_pcoa_weighted$SampleID <- rownames(df_pcoa_weighted)

# Juntar com metadata (que tem SampleSheet)
df_pcoa_weighted <- merge(df_pcoa_weighted, metadata[, c("sample.id", "SampleSheet")],
                          by.x = "SampleID", by.y = "sample.id", all.x = TRUE)

metadata$SampleSheet <- NA  # Cria a coluna vazia

metadata$SampleSheet[1:19] <- "Batch1"
metadata$SampleSheet[20:109] <- "Batch2"
metadata$SampleSheet[110:130] <- "Batch3"



ggplot(df_pcoa_weighted, aes(x = Axis.1, y = Axis.2, color = SampleSheet)) +
  geom_point(size = 2) +
  labs(
    title = "PCoA - Weighted UniFrac (Decontam + Rarefied)",
    x = paste0("Axis 1 [", round(ord_weighted$values$Relative_eig[1] * 100, 1), "%]"),
    y = paste0("Axis 2 [", round(ord_weighted$values$Relative_eig[2] * 100, 1), "%]"),
    color = "Batch"
  ) +
  theme_minimal()


#===> unweighted

#Realizar pcoa
ord_unweighted <- ordinate(physeq_rarefied, method = "PCoA", distance = unifrac_unweighted)

#Transformar em dataframe e adicionar nomes
df_pcoa_unweighted <- as.data.frame(ord_unweighted$vectors)
df_pcoa_unweighted$SampleID <- rownames(df_pcoa_unweighted)



df_pcoa_unweighted$SampleID <- rownames(df_pcoa_unweighted)

# Juntar com metadata (que tem SampleSheet)
df_pcoa_unweighted <- merge(df_pcoa_unweighted, metadata[, c("sample.id", "SampleSheet")],
                          by.x = "SampleID", by.y = "sample.id", all.x = TRUE)





ggplot(df_pcoa_unweighted, aes(x = Axis.1, y = Axis.2, color = SampleSheet)) +
  geom_point(size = 2) +
  labs(
    title = "PCoA - unweighted UniFrac (Decontam + Rarefied)",
    x = paste0("Axis 1 [", round(ord_unweighted$values$Relative_eig[1] * 100, 1), "%]"),
    y = paste0("Axis 2 [", round(ord_unweighted$values$Relative_eig[2] * 100, 1), "%]"),
    color = "Batch"
  ) +
  theme_minimal()



#============================================================================#
#                    Threshold .05
#============================================================================#

#=============================================================#
#                Decontam S1: threshold .01
#=============================================================#
#contam_s1 <- isContaminant(ps_s1, method = "prevalence", neg = "is.neg", threshold = 0.5)  


# 4. Criar coluna lógica com os negativos
sample_data(ps_s1)$is.neg <- sample_names(ps_s1) %in% neg_controls_s1

# 5. Rodar o decontam (método prevalência)
contam_s1 <- isContaminant(ps_s1, method = "prevalence", neg = "is.neg")  

# 6. Verificar quantos contaminantes foram encontrados
table(contam_s1$contaminant)

# 7. Filtrar ASVs não contaminantes
asvs_s1 <- taxa_names(ps_s1)[!contam_s1$contaminant]

# 8. Criar novo objeto phyloseq sem os contaminantes
ps_s1_clean <- prune_taxa(asvs_s1, ps_s1)

# 9. Verificação opcional
ps_s1_clean




#=======> CONTAMINANTES E SUA TAXONOMIA:

# Taxonomia dos contaminantes
contam_asvs <- rownames(contam_s1)[contam_s1$contaminant]
contam_tax <- tax_table(ps)[contam_asvs, ]

# Em formato data.frame
contam_tax_df <- as.data.frame(contam_tax)
print(contam_tax_df)



#write.csv(contam_tax_df, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/contaminantes_sampleSheet1.csv", row.names = TRUE)





#=> Quantos reads se perderam?

# 1. Total de reads por amostra antes e depois
reads_before <- sample_sums(ps_s1)
reads_after <- sample_sums(ps_s1_clean)

# 2. Calcular reads removidos e porcentagem
reads_removed <- reads_before - reads_after
percent_removed <- (reads_removed / reads_before) * 100

# 3. Criar dataframe com os resultados
removed_df <- data.frame(
  Sample = names(reads_before),
  Reads_Before = reads_before,
  Reads_After = reads_after,
  Reads_Removed = reads_removed,
  Percent_Removed = round(percent_removed, 2)
)

# 4. Visualizar os dados
removed_df <- removed_df[order(removed_df$Percent_Removed, decreasing = TRUE), ]
print(removed_df)


# 4. Remover os controles (NEG e POS)
removed_df_real <- removed_df %>%
  filter(!grepl("NEG|POS", Sample))

# 5. Reordenar Sample com base no % removido
removed_df_real$Sample <- factor(
  removed_df_real$Sample,
  levels = removed_df_real$Sample[order(removed_df_real$Percent_Removed, decreasing = TRUE)]
)

# 6. Criar gráfico de barras
ggplot(removed_df_real, aes(x = Sample, y = Percent_Removed)) +
  geom_bar(stat = "identity", fill = "#4682B4") +
  labs(
    title = "Percentual de Reads Removidos por Amostra (S1, sem controles). Threshold .01",
    x = "Amostra",
    y = "% de Reads Removidos"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )





#=============================================================#
#                 Decontam S2: threshold .05 
#=============================================================#

library(phyloseq)
library(decontam)


# 4. Criar a coluna lógica com os negativos
sample_data(ps_s2)$is.neg <- sample_names(ps_s2) %in% neg_controls_s2

# 5. Rodar o decontam (método prevalência)
contam_s2 <- isContaminant(ps_s2, method = "prevalence", neg = "is.neg")

# 6. Verificar quantos contaminantes foram identificados
table(contam_s2$contaminant)

# 7. Obter ASVs que não são contaminantes
asvs_s2 <- taxa_names(ps_s2)[!contam_s2$contaminant]

# 8. Criar novo objeto phyloseq sem os contaminantes
ps_s2_clean <- prune_taxa(asvs_s2, ps_s2)

# 9. Confirmar visualmente
ps_s2_clean

#VERIFICAR OS CONTAMINANTES

# 1. Identificar os ASVs considerados contaminantes
contam_asvs_s2 <- rownames(contam_s2)[contam_s2$contaminant]

# 2. Ver a taxonomia dos contaminantes
contam_tax_s2 <- tax_table(ps)[contam_asvs_s2, ]
contam_tax_s2_df <- as.data.frame(contam_tax_s2)

# 3. Visualizar
print(contam_tax_s2_df)



write.csv(contam_tax_s2_df, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/contaminantes_s2.csv", row.names = TRUE)

#generos dos contaminantes: 

table(contam_tax_s2_df$Genus)


#=> Quantos reads se perderam?

# 1. Total de reads por amostra antes e depois
reads_before2 <- sample_sums(ps_s2)
reads_after2 <- sample_sums(ps_s2_clean)

# 2. Calcular reads removidos e porcentagem
reads_removed2 <- reads_before2 - reads_after2
percent_removed2 <- (reads_removed2 / reads_before2) * 100

# 3. Criar dataframe com os resultados
removed_df2 <- data.frame(
  Sample = names(reads_before2),
  Reads_Before2 = reads_before2,
  Reads_After2 = reads_after2,
  Reads_Removed2 = reads_removed2,  # << CORRIGIDO AQUI
  Percent_Removed2 = round(percent_removed2, 2)
)

# 4. Visualizar os dados ordenados por % removida
removed_df2 <- removed_df2[order(removed_df2$Percent_Removed2, decreasing = TRUE), ]
print(removed_df2)



library(ggplot2)
library(dplyr)

# 1. Filtrar amostras reais (sem NEG ou POS)
removed_df2_real <- removed_df2 %>%
  filter(!grepl("NEG|POS", Sample))

# 2. Ordenar as amostras pelo percentual removido
removed_df2_real$Sample <- factor(
  removed_df2_real$Sample,
  levels = removed_df2_real$Sample[order(removed_df2_real$Percent_Removed2, decreasing = TRUE)]
)

# 3. Criar o gráfico com amostras ordenadas
ggplot(removed_df2_real, aes(x = Sample, y = Percent_Removed2)) +
  geom_bar(stat = "identity", fill = "#4682B4") +
  labs(
    title = "Percentual de Reads Removidos por Amostra (sem controles)",
    x = "Amostra",
    y = "% de Reads Removidos"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )



#=============================================================#
#                        Decontam S3
#=============================================================#

library(phyloseq)
library(decontam)



# 3. Subset do objeto phyloseq com amostras e controles
ps_s3 <- prune_samples(sample_names(ps) %in% c(samples_s3, neg_controls_s3), physeq)

# 4. Marcar os controles negativos
sample_data(ps_s3)$is.neg <- sample_names(ps_s3) %in% neg_controls_s3


#2. O limiar (threshold) do decontam pode estar conservador
#Por padrão, o threshold = 0.1, o que pode ser pouco sensível para detectar contaminantes que também aparecem em amostras reais.
#testar uma detecção mais rigorosa, pode tentar:

# 5. Rodar o método de prevalência do decontam
contam_s3 <- isContaminant(ps_s3, method = "prevalence", neg = "is.neg")

# 6. Verificar quantos contaminantes foram detectados
table(contam_s3$contaminant)

# 7. Obter os ASVs não contaminantes
asvs_s3 <- taxa_names(ps_s3)[!contam_s3$contaminant]

# 8. Criar objeto phyloseq limpo
ps_s3_clean <- prune_taxa(asvs_s3, ps_s3)

# 9. Verificação
ps_s3_clean

#TAXONOMIA DOS CONTAMINANTES

# 1. Identificar os ASVs considerados contaminantes
contam_asvs_s3 <- rownames(contam_s3)[contam_s3$contaminant]

# 2. Obter a taxonomia desses ASVs
contam_tax_s3 <- tax_table(ps)[contam_asvs_s3, ]
contam_tax_s3_df <- as.data.frame(contam_tax_s3)

# 3. Visualizar
print(contam_tax_s3_df)






write.csv(contam_tax_s3_df, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/contaminantes_s3.csv", row.names = TRUE)


#outras asvs a verificar: f__Mycoplasmataceae (2b9c8fa41ae2b255b4f3831f752010e5)  f__Succinivibrionaceae; g__Ruminobacter; s__gut_metagenome (3a0c5873e8683c4c3cec68faa75b3bb0) 

asvs_verificar <- c("2b9c8fa41ae2b255b4f3831f752010e5", "3a0c5873e8683c4c3cec68faa75b3bb0")




#=====> verificar os controles

controles <- sample_names(ps_s3)[grepl("NEG|Neg|PCRneg|Control", sample_names(ps_s3))]
ps_s3_controles <- subset_samples(ps_s3, sample_names(ps_s3) %in% controles)

#Ver composição (ex: por Genus ou Phylum):

# Transformar em proporção relativa
ps_s3_controles_rel <- transform_sample_counts(ps_s3_controles, function(x) x / sum(x))

library(phyloseq)
library(ggplot2)
library(dplyr)

# Derreter objeto phyloseq para data.frame
df_genus <- psmelt(genus_abund)

# Selecionar os 20 gêneros mais abundantes
top_genus <- df_genus %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 20) %>%
  pull(Genus)

# Filtrar apenas os 20 gêneros mais abundantes
df_top <- df_genus %>% filter(Genus %in% top_genus)

# Gráfico
ggplot(df_top, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundância Relativa por Gênero nos Controles",
       x = "Amostra", y = "Abundância Relativa") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Extrair abundância de Treponema nos controles
df_treponema <- psmelt(ps_s3_controles_rel) %>%
  filter(Genus == "Treponema")

# Visualizar
print(df_treponema)

# Se quiser um gráfico:
ggplot(df_treponema, aes(x = Sample, y = Abundance, fill = Sample)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundância de Treponema nos Controles",
       x = "Amostra", y = "Abundância Relativa") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = "none")



#=> Quantos reads se perderam?

# 1. Total de reads por amostra antes e depois
reads_before3 <- sample_sums(ps_s3)
reads_after3 <- sample_sums(ps_s3_clean)

# 3. Calcular reads removidos e porcentagem
reads_removed3 <- reads_before3 - reads_after3
percent_removed3 <- (reads_removed3 / reads_before3) * 100

# 3. Criar dataframe com os resultados
removed_df3 <- data.frame(
  Sample = names(reads_before3),
  Reads_Before3 = reads_before3,
  Reads_After3 = reads_after3,
  Reads_Removed3 = reads_removed3,  # << CORRIGIDO AQUI
  Percent_Removed3 = round(percent_removed3, 2)
)

# 4. Visualizar os dados ordenados por % removida
removed_df3 <- removed_df3[order(removed_df3$Percent_Removed3, decreasing = TRUE), ]
print(removed_df3)


# 1. Filtrar amostras reais (excluindo controles)
removed_df3_real <- removed_df3 %>%
  filter(
    !grepl("NEG|POS", Sample) & 
      !grepl("20240325_Neg_Control", Sample)
  )

# 2. Ordenar as amostras pelo percentual removido
removed_df3_real$Sample <- factor(
  removed_df3_real$Sample,
  levels = removed_df3_real$Sample[order(removed_df3_real$Percent_Removed3, decreasing = TRUE)]
)

# 3. Criar o gráfico com amostras ordenadas
ggplot(removed_df3_real, aes(x = Sample, y = Percent_Removed3)) +
  geom_bar(stat = "identity", fill = "#4682B4") +
  labs(
    title = "Percentual de Reads Removidos por Amostra (sem controles)",
    x = "Amostra",
    y = "% de Reads Removidos"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )





#=============================================================#
#                Juntar os três objetos limpos
#=============================================================#

# Juntar os objetos phyloseq limpos
physeq_descontaminado <- merge_phyloseq(ps_s1_clean, ps_s2_clean, ps_s3_clean)

# Verificar dimensões do objeto combinado
physeq_descontaminado


#=====> rarefação

#Verificar a distribuição de reads
#Antes de definir um ponto de corte para a rarefação:


hist(sample_sums(physeq_descontaminado),
     breaks = 50,
     main = "Distribuição de Reads por Amostra",
     xlab = "Número de Reads",
     col = "#66c2a5")


reads_por_amostra <- sample_sums(physeq_descontaminado)
reads_por_amostra[which.min(reads_por_amostra)]

#> reads_por_amostra[which.min(reads_por_amostra)]
#30081.F00 
#17693 

#tirei a amostra 30081.F00. Rarefação será com 25689 para maior profundidade e para manter amostra 10391

set.seed(123)
physeq_rarefied <- rarefy_even_depth(
  physeq_descontaminado,
  sample.size = 25689,
  rngseed = 123,
  verbose = TRUE
)

#verificar
sample_sums(physeq_rarefied)




#=================================================================================#
#       Criar beta diversidade
#=================================================================================#

#verificar se posso usar a tree sem podar



# ASVs que estão no physeq mas não estão na árvore
setdiff(taxa_names(physeq_rarefied), tree$tip.label) # => se retornar vazio, está ok!



# Adicionar a árvore ao objeto phyloseq
phy_tree(physeq_rarefied) <- tree



# Calcular UniFrac ponderado (weighted)
unifrac_weighted <- phyloseq::UniFrac(physeq_rarefied, weighted = TRUE)

# Calcular UniFrac não ponderado (unweighted)
unifrac_unweighted <- phyloseq::UniFrac(physeq_rarefied, weighted = FALSE, normalized = TRUE, parallel = TRUE, fast = TRUE)


#Realizar pcoa
ord_weighted <- ordinate(physeq_rarefied, method = "PCoA", distance = unifrac_weighted)

#Transformar em dataframe e adicionar nomes
df_pcoa_weighted <- as.data.frame(ord_weighted$vectors)
df_pcoa_weighted$SampleID <- rownames(df_pcoa_weighted)

#plotar

library(ggplot2)

ggplot(df_pcoa_weighted, aes(x = Axis.1, y = Axis.2)) +
  geom_point() +
  labs(
    title = "PCoA - Weighted UniFrac (sem outlier)",
    x = paste0("Axis 1 [", round(ord_weighted$values$Relative_eig[1] * 100, 1), "%]"),
    y = paste0("Axis 2 [", round(ord_weighted$values$Relative_eig[2] * 100, 1), "%]")
  ) +
  theme_minimal()



#nomes nos pontos

library(ggplot2)
library(ggrepel)

g <- ggplot(df_pcoa_weighted, aes(x = Axis.1, y = Axis.2)) +
  geom_point(size = 2) +
  geom_text(aes(label = Sample), size = 3, vjust = -0.5, hjust = 0.5) +
  labs(
    title = "PCoA - Weighted UniFrac (sem outlier)",
    x = "Axis 1",
    y = "Axis 2"
  ) +
  theme_minimal()

print(g)


df_pcoa_weighted$SampleID <- rownames(df_pcoa_weighted)

# Juntar com metadata (que tem SampleSheet)
df_pcoa_weighted <- merge(df_pcoa_weighted, metadata[, c("sample.id", "SampleSheet")],
                          by.x = "SampleID", by.y = "sample.id", all.x = TRUE)

metadata$SampleSheet <- NA  # Cria a coluna vazia

metadata$SampleSheet[1:19] <- "Batch1"
metadata$SampleSheet[20:109] <- "Batch2"
metadata$SampleSheet[110:130] <- "Batch3"



ggplot(df_pcoa_weighted, aes(x = Axis.1, y = Axis.2, color = SampleSheet)) +
  geom_point(size = 2) +
  labs(
    title = "PCoA - Weighted UniFrac (Decontam + Rarefied)",
    x = paste0("Axis 1 [", round(ord_weighted$values$Relative_eig[1] * 100, 1), "%]"),
    y = paste0("Axis 2 [", round(ord_weighted$values$Relative_eig[2] * 100, 1), "%]"),
    color = "Batch"
  ) +
  theme_minimal()


#===> unweighted

#Realizar pcoa
ord_unweighted <- ordinate(physeq_rarefied, method = "PCoA", distance = unifrac_unweighted)

#Transformar em dataframe e adicionar nomes
df_pcoa_unweighted <- as.data.frame(ord_unweighted$vectors)
df_pcoa_unweighted$SampleID <- rownames(df_pcoa_unweighted)



df_pcoa_unweighted$SampleID <- rownames(df_pcoa_unweighted)

# Juntar com metadata (que tem SampleSheet)
df_pcoa_unweighted <- merge(df_pcoa_unweighted, metadata[, c("sample.id", "SampleSheet")],
                            by.x = "SampleID", by.y = "sample.id", all.x = TRUE)





ggplot(df_pcoa_unweighted, aes(x = Axis.1, y = Axis.2, color = SampleSheet)) +
  geom_point(size = 2) +
  labs(
    title = "PCoA - unweighted UniFrac (Decontam + Rarefied)",
    x = paste0("Axis 1 [", round(ord_unweighted$values$Relative_eig[1] * 100, 1), "%]"),
    y = paste0("Axis 2 [", round(ord_unweighted$values$Relative_eig[2] * 100, 1), "%]"),
    color = "Batch"
  ) +
  theme_minimal()

#....................................................................................#
#============================================================================#
#                    Threshold 1
#============================================================================#

#=============================================================#
#                Decontam S1: threshold .01
#=============================================================#
#contam_s1 <- isContaminant(ps_s1, method = "prevalence", neg = "is.neg", threshold = 0.5)  


# 4. Criar coluna lógica com os negativos
sample_data(ps_s1)$is.neg <- sample_names(ps_s1) %in% neg_controls_s1

# 5. Rodar o decontam (método prevalência)
contam_s1 <- isContaminant(ps_s1, method = "prevalence", neg = "is.neg")  

# 6. Verificar quantos contaminantes foram encontrados
table(contam_s1$contaminant)

# 7. Filtrar ASVs não contaminantes
asvs_s1 <- taxa_names(ps_s1)[!contam_s1$contaminant]

# 8. Criar novo objeto phyloseq sem os contaminantes
ps_s1_clean <- prune_taxa(asvs_s1, ps_s1)

# 9. Verificação opcional
ps_s1_clean




#=======> CONTAMINANTES E SUA TAXONOMIA:

# Taxonomia dos contaminantes
contam_asvs <- rownames(contam_s1)[contam_s1$contaminant]
contam_tax <- tax_table(ps)[contam_asvs, ]

# Em formato data.frame
contam_tax_df <- as.data.frame(contam_tax)
print(contam_tax_df)



#write.csv(contam_tax_df, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/contaminantes_sampleSheet1.csv", row.names = TRUE)





#=> Quantos reads se perderam?

# 1. Total de reads por amostra antes e depois
reads_before <- sample_sums(ps_s1)
reads_after <- sample_sums(ps_s1_clean)

# 2. Calcular reads removidos e porcentagem
reads_removed <- reads_before - reads_after
percent_removed <- (reads_removed / reads_before) * 100

# 3. Criar dataframe com os resultados
removed_df <- data.frame(
  Sample = names(reads_before),
  Reads_Before = reads_before,
  Reads_After = reads_after,
  Reads_Removed = reads_removed,
  Percent_Removed = round(percent_removed, 2)
)

# 4. Visualizar os dados
removed_df <- removed_df[order(removed_df$Percent_Removed, decreasing = TRUE), ]
print(removed_df)


# 4. Remover os controles (NEG e POS)
removed_df_real <- removed_df %>%
  filter(!grepl("NEG|POS", Sample))

# 5. Reordenar Sample com base no % removido
removed_df_real$Sample <- factor(
  removed_df_real$Sample,
  levels = removed_df_real$Sample[order(removed_df_real$Percent_Removed, decreasing = TRUE)]
)

# 6. Criar gráfico de barras
ggplot(removed_df_real, aes(x = Sample, y = Percent_Removed)) +
  geom_bar(stat = "identity", fill = "#4682B4") +
  labs(
    title = "Percentual de Reads Removidos por Amostra (S1, sem controles). Threshold .01",
    x = "Amostra",
    y = "% de Reads Removidos"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )





#=============================================================#
#                 Decontam S2: threshold .05 
#=============================================================#

library(phyloseq)
library(decontam)


# 4. Criar a coluna lógica com os negativos
sample_data(ps_s2)$is.neg <- sample_names(ps_s2) %in% neg_controls_s2

# 5. Rodar o decontam (método prevalência)
contam_s2 <- isContaminant(ps_s2, method = "prevalence", neg = "is.neg")

# 6. Verificar quantos contaminantes foram identificados
table(contam_s2$contaminant)

# 7. Obter ASVs que não são contaminantes
asvs_s2 <- taxa_names(ps_s2)[!contam_s2$contaminant]

# 8. Criar novo objeto phyloseq sem os contaminantes
ps_s2_clean <- prune_taxa(asvs_s2, ps_s2)

# 9. Confirmar visualmente
ps_s2_clean

#VERIFICAR OS CONTAMINANTES

# 1. Identificar os ASVs considerados contaminantes
contam_asvs_s2 <- rownames(contam_s2)[contam_s2$contaminant]

# 2. Ver a taxonomia dos contaminantes
contam_tax_s2 <- tax_table(ps)[contam_asvs_s2, ]
contam_tax_s2_df <- as.data.frame(contam_tax_s2)

# 3. Visualizar
print(contam_tax_s2_df)



write.csv(contam_tax_s2_df, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/contaminantes_s2.csv", row.names = TRUE)

#generos dos contaminantes: 

table(contam_tax_s2_df$Genus)


#=> Quantos reads se perderam?

# 1. Total de reads por amostra antes e depois
reads_before2 <- sample_sums(ps_s2)
reads_after2 <- sample_sums(ps_s2_clean)

# 2. Calcular reads removidos e porcentagem
reads_removed2 <- reads_before2 - reads_after2
percent_removed2 <- (reads_removed2 / reads_before2) * 100

# 3. Criar dataframe com os resultados
removed_df2 <- data.frame(
  Sample = names(reads_before2),
  Reads_Before2 = reads_before2,
  Reads_After2 = reads_after2,
  Reads_Removed2 = reads_removed2,  # << CORRIGIDO AQUI
  Percent_Removed2 = round(percent_removed2, 2)
)

# 4. Visualizar os dados ordenados por % removida
removed_df2 <- removed_df2[order(removed_df2$Percent_Removed2, decreasing = TRUE), ]
print(removed_df2)



library(ggplot2)
library(dplyr)

# 1. Filtrar amostras reais (sem NEG ou POS)
removed_df2_real <- removed_df2 %>%
  filter(!grepl("NEG|POS", Sample))

# 2. Ordenar as amostras pelo percentual removido
removed_df2_real$Sample <- factor(
  removed_df2_real$Sample,
  levels = removed_df2_real$Sample[order(removed_df2_real$Percent_Removed2, decreasing = TRUE)]
)

# 3. Criar o gráfico com amostras ordenadas
ggplot(removed_df2_real, aes(x = Sample, y = Percent_Removed2)) +
  geom_bar(stat = "identity", fill = "#4682B4") +
  labs(
    title = "Percentual de Reads Removidos por Amostra (sem controles)",
    x = "Amostra",
    y = "% de Reads Removidos"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )



#=============================================================#
#                        Decontam S3
#=============================================================#

library(phyloseq)
library(decontam)



# 3. Subset do objeto phyloseq com amostras e controles
ps_s3 <- prune_samples(sample_names(ps) %in% c(samples_s3, neg_controls_s3), physeq)

# 4. Marcar os controles negativos
sample_data(ps_s3)$is.neg <- sample_names(ps_s3) %in% neg_controls_s3


#2. O limiar (threshold) do decontam pode estar conservador
#Por padrão, o threshold = 0.1, o que pode ser pouco sensível para detectar contaminantes que também aparecem em amostras reais.
#testar uma detecção mais rigorosa, pode tentar:

# 5. Rodar o método de prevalência do decontam
contam_s3 <- isContaminant(ps_s3, method = "prevalence", neg = "is.neg")

# 6. Verificar quantos contaminantes foram detectados
table(contam_s3$contaminant)

# 7. Obter os ASVs não contaminantes
asvs_s3 <- taxa_names(ps_s3)[!contam_s3$contaminant]

# 8. Criar objeto phyloseq limpo
ps_s3_clean <- prune_taxa(asvs_s3, ps_s3)

# 9. Verificação
ps_s3_clean

#TAXONOMIA DOS CONTAMINANTES

# 1. Identificar os ASVs considerados contaminantes
contam_asvs_s3 <- rownames(contam_s3)[contam_s3$contaminant]

# 2. Obter a taxonomia desses ASVs
contam_tax_s3 <- tax_table(ps)[contam_asvs_s3, ]
contam_tax_s3_df <- as.data.frame(contam_tax_s3)

# 3. Visualizar
print(contam_tax_s3_df)






write.csv(contam_tax_s3_df, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/contaminantes_s3.csv", row.names = TRUE)


#outras asvs a verificar: f__Mycoplasmataceae (2b9c8fa41ae2b255b4f3831f752010e5)  f__Succinivibrionaceae; g__Ruminobacter; s__gut_metagenome (3a0c5873e8683c4c3cec68faa75b3bb0) 

asvs_verificar <- c("2b9c8fa41ae2b255b4f3831f752010e5", "3a0c5873e8683c4c3cec68faa75b3bb0")




#=====> verificar os controles

controles <- sample_names(ps_s3)[grepl("NEG|Neg|PCRneg|Control", sample_names(ps_s3))]
ps_s3_controles <- subset_samples(ps_s3, sample_names(ps_s3) %in% controles)

#Ver composição (ex: por Genus ou Phylum):

# Transformar em proporção relativa
ps_s3_controles_rel <- transform_sample_counts(ps_s3_controles, function(x) x / sum(x))

library(phyloseq)
library(ggplot2)
library(dplyr)

# Derreter objeto phyloseq para data.frame
df_genus <- psmelt(genus_abund)

# Selecionar os 20 gêneros mais abundantes
top_genus <- df_genus %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 20) %>%
  pull(Genus)

# Filtrar apenas os 20 gêneros mais abundantes
df_top <- df_genus %>% filter(Genus %in% top_genus)

# Gráfico
ggplot(df_top, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundância Relativa por Gênero nos Controles",
       x = "Amostra", y = "Abundância Relativa") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Extrair abundância de Treponema nos controles
df_treponema <- psmelt(ps_s3_controles_rel) %>%
  filter(Genus == "Treponema")

# Visualizar
print(df_treponema)

# Se quiser um gráfico:
ggplot(df_treponema, aes(x = Sample, y = Abundance, fill = Sample)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundância de Treponema nos Controles",
       x = "Amostra", y = "Abundância Relativa") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = "none")



#=> Quantos reads se perderam?

# 1. Total de reads por amostra antes e depois
reads_before3 <- sample_sums(ps_s3)
reads_after3 <- sample_sums(ps_s3_clean)

# 3. Calcular reads removidos e porcentagem
reads_removed3 <- reads_before3 - reads_after3
percent_removed3 <- (reads_removed3 / reads_before3) * 100

# 3. Criar dataframe com os resultados
removed_df3 <- data.frame(
  Sample = names(reads_before3),
  Reads_Before3 = reads_before3,
  Reads_After3 = reads_after3,
  Reads_Removed3 = reads_removed3,  # << CORRIGIDO AQUI
  Percent_Removed3 = round(percent_removed3, 2)
)

# 4. Visualizar os dados ordenados por % removida
removed_df3 <- removed_df3[order(removed_df3$Percent_Removed3, decreasing = TRUE), ]
print(removed_df3)


# 1. Filtrar amostras reais (excluindo controles)
removed_df3_real <- removed_df3 %>%
  filter(
    !grepl("NEG|POS", Sample) & 
      !grepl("20240325_Neg_Control", Sample)
  )

# 2. Ordenar as amostras pelo percentual removido
removed_df3_real$Sample <- factor(
  removed_df3_real$Sample,
  levels = removed_df3_real$Sample[order(removed_df3_real$Percent_Removed3, decreasing = TRUE)]
)

# 3. Criar o gráfico com amostras ordenadas
ggplot(removed_df3_real, aes(x = Sample, y = Percent_Removed3)) +
  geom_bar(stat = "identity", fill = "#4682B4") +
  labs(
    title = "Percentual de Reads Removidos por Amostra (sem controles)",
    x = "Amostra",
    y = "% de Reads Removidos"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )





#=============================================================#
#                Juntar os três objetos limpos
#=============================================================#

# Juntar os objetos phyloseq limpos
physeq_descontaminado <- merge_phyloseq(ps_s1_clean, ps_s2_clean, ps_s3_clean)

# Verificar dimensões do objeto combinado
physeq_descontaminado


#=====> rarefação

#Verificar a distribuição de reads
#Antes de definir um ponto de corte para a rarefação:


hist(sample_sums(physeq_descontaminado),
     breaks = 50,
     main = "Distribuição de Reads por Amostra",
     xlab = "Número de Reads",
     col = "#66c2a5")


reads_por_amostra <- sample_sums(physeq_descontaminado)
reads_por_amostra[which.min(reads_por_amostra)]

#> reads_por_amostra[which.min(reads_por_amostra)]
#30081.F00 
#17693 

#tirei a amostra 30081.F00. Rarefação será com 25689 para maior profundidade e para manter amostra 10391

set.seed(123)
physeq_rarefied <- rarefy_even_depth(
  physeq_descontaminado,
  sample.size = 25689,
  rngseed = 123,
  verbose = TRUE
)

#verificar
sample_sums(physeq_rarefied)




#=================================================================================#
#       Criar beta diversidade
#=================================================================================#

#verificar se posso usar a tree sem podar



# ASVs que estão no physeq mas não estão na árvore
setdiff(taxa_names(physeq_rarefied), tree$tip.label) # => se retornar vazio, está ok!



# Adicionar a árvore ao objeto phyloseq
phy_tree(physeq_rarefied) <- tree



# Calcular UniFrac ponderado (weighted)
unifrac_weighted <- phyloseq::UniFrac(physeq_rarefied, weighted = TRUE)

# Calcular UniFrac não ponderado (unweighted)
unifrac_unweighted <- phyloseq::UniFrac(physeq_rarefied, weighted = FALSE, normalized = TRUE, parallel = TRUE, fast = TRUE)


#Realizar pcoa
ord_weighted <- ordinate(physeq_rarefied, method = "PCoA", distance = unifrac_weighted)

#Transformar em dataframe e adicionar nomes
df_pcoa_weighted <- as.data.frame(ord_weighted$vectors)
df_pcoa_weighted$SampleID <- rownames(df_pcoa_weighted)

#plotar

library(ggplot2)

ggplot(df_pcoa_weighted, aes(x = Axis.1, y = Axis.2)) +
  geom_point() +
  labs(
    title = "PCoA - Weighted UniFrac (sem outlier)",
    x = paste0("Axis 1 [", round(ord_weighted$values$Relative_eig[1] * 100, 1), "%]"),
    y = paste0("Axis 2 [", round(ord_weighted$values$Relative_eig[2] * 100, 1), "%]")
  ) +
  theme_minimal()



#nomes nos pontos

library(ggplot2)
library(ggrepel)

g <- ggplot(df_pcoa_weighted, aes(x = Axis.1, y = Axis.2)) +
  geom_point(size = 2) +
  geom_text(aes(label = Sample), size = 3, vjust = -0.5, hjust = 0.5) +
  labs(
    title = "PCoA - Weighted UniFrac (sem outlier)",
    x = "Axis 1",
    y = "Axis 2"
  ) +
  theme_minimal()

print(g)


df_pcoa_weighted$SampleID <- rownames(df_pcoa_weighted)

# Juntar com metadata (que tem SampleSheet)
df_pcoa_weighted <- merge(df_pcoa_weighted, metadata[, c("sample.id", "SampleSheet")],
                          by.x = "SampleID", by.y = "sample.id", all.x = TRUE)

metadata$SampleSheet <- NA  # Cria a coluna vazia

metadata$SampleSheet[1:19] <- "Batch1"
metadata$SampleSheet[20:109] <- "Batch2"
metadata$SampleSheet[110:130] <- "Batch3"



ggplot(df_pcoa_weighted, aes(x = Axis.1, y = Axis.2, color = SampleSheet)) +
  geom_point(size = 2) +
  labs(
    title = "PCoA - Weighted UniFrac (Decontam + Rarefied)",
    x = paste0("Axis 1 [", round(ord_weighted$values$Relative_eig[1] * 100, 1), "%]"),
    y = paste0("Axis 2 [", round(ord_weighted$values$Relative_eig[2] * 100, 1), "%]"),
    color = "Batch"
  ) +
  theme_minimal()


#===> unweighted

#Realizar pcoa
ord_unweighted <- ordinate(physeq_rarefied, method = "PCoA", distance = unifrac_unweighted)

#Transformar em dataframe e adicionar nomes
df_pcoa_unweighted <- as.data.frame(ord_unweighted$vectors)
df_pcoa_unweighted$SampleID <- rownames(df_pcoa_unweighted)



df_pcoa_unweighted$SampleID <- rownames(df_pcoa_unweighted)

# Juntar com metadata (que tem SampleSheet)
df_pcoa_unweighted <- merge(df_pcoa_unweighted, metadata[, c("sample.id", "SampleSheet")],
                            by.x = "SampleID", by.y = "sample.id", all.x = TRUE)





ggplot(df_pcoa_unweighted, aes(x = Axis.1, y = Axis.2, color = SampleSheet)) +
  geom_point(size = 2) +
  labs(
    title = "PCoA - unweighted UniFrac (Decontam + Rarefied)",
    x = paste0("Axis 1 [", round(ord_unweighted$values$Relative_eig[1] * 100, 1), "%]"),
    y = paste0("Axis 2 [", round(ord_unweighted$values$Relative_eig[2] * 100, 1), "%]"),
    color = "Batch"
  ) +
  theme_minimal()
