#=============================================================#
#                        Decontam S1: threshold .05
#=============================================================#


# 4. Criar coluna lógica com os negativos
sample_data(ps_s1)$is.neg <- sample_names(ps_s1) %in% neg_controls_s1

# 5. Rodar o decontam (método prevalência)
contam_s1_ <- isContaminant(ps_s1, method = "prevalence", neg = "is.neg", threshold= 0.5)  

# 6. Verificar quantos contaminantes foram encontrados
table(contam_s1_$contaminant)

# 7. Filtrar ASVs não contaminantes
asvs_s1_ <- taxa_names(ps_s1)[!contam_s1_$contaminant]

# 8. Criar novo objeto phyloseq sem os contaminantes
ps_s1_clean_ <- prune_taxa(asvs_s1_, ps_s1)

# 9. Verificação opcional
ps_s1_clean_




#=======> CONTAMINANTES E SUA TAXONOMIA:

# Taxonomia dos contaminantes
contam_asvs_ <- rownames(contam_s1_)[contam_s1_$contaminant]
contam_tax_ <- tax_table(ps)[contam_asvs_, ]

# Em formato data.frame
contam_tax_df_ <- as.data.frame(contam_tax_)
print(contam_tax_df_)



#write.csv(contam_tax_df, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/contaminantes_sampleSheet1.csv", row.names = TRUE)





#=> Quantos reads se perderam?

# 1. Total de reads por amostra antes e depois
reads_before_ <- sample_sums(ps_s1)
reads_after_ <- sample_sums(ps_s1_clean_)

# 2. Calcular reads removidos e porcentagem
reads_removed_ <- reads_before_ - reads_after_
percent_removed_ <- (reads_removed_ / reads_before_) * 100

# 3. Criar dataframe com os resultados
removed_df_ <- data.frame(
  Sample = names(reads_before_),
  Reads_Before_ = reads_before_,
  Reads_After_ = reads_after_,
  Reads_Removed_ = reads_removed_,
  Percent_Removed_ = round(percent_removed_, 2)
)

# 4. Visualizar os dados
removed_df_ <- removed_df_[order(removed_df_$Percent_Removed_, decreasing = TRUE), ]
print(removed_df_)


# 4. Remover os controles (NEG e POS)
removed_df_real_ <- removed_df_ %>%
  filter(!grepl("NEG|POS", Sample))

# 5. Reordenar Sample com base no % removido
removed_df_real_$Sample <- factor(
  removed_df_real_$Sample,
  levels = removed_df_real_$Sample[order(removed_df_real_$Percent_Removed_, decreasing = TRUE)]
)

# 6. Criar gráfico de barras
# 6. Criar gráfico de barras
ggplot(removed_df_real_, aes(x = Sample, y = Percent_Removed_)) +
  geom_bar(stat = "identity", fill = "#4682B4") +
  labs(
    title = "Percentual de Reads Removidos por Amostra (S1, sem controles)Threshold .05",
    x = "Amostra",
    y = "% de Reads Removidos"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )


#=============================================================#
# S1: comparar visualmente o percentual de reads removidos por 
#amostra entre os dois thresholds (0.01 e 0.05
#=============================================================#



# 1. Adicionar coluna indicando o threshold
removed_df_real$Threshold <- "0.01"
removed_df_real_$Threshold <- "0.05"

# 2. Renomear coluna Percent_Removed_ para Percent_Removed para unir
colnames(removed_df_real_)[colnames(removed_df_real_) == "Percent_Removed_"] <- "Percent_Removed"
colnames(removed_df_real_)[colnames(removed_df_real_) == "Reads_Before_"] <- "Reads_Before"
colnames(removed_df_real_)[colnames(removed_df_real_) == "Reads_After_"] <- "Reads_After"
colnames(removed_df_real_)[colnames(removed_df_real_) == "Reads_Removed_"] <- "Reads_Removed"

# 3. Unir os dois dataframes
df_comparacao <- bind_rows(removed_df_real, removed_df_real_)

# 4. Reordenar amostras por média de remoção (opcional)
df_comparacao$Sample <- factor(
  df_comparacao$Sample,
  levels = df_comparacao %>%
    group_by(Sample) %>%
    summarise(media = mean(Percent_Removed)) %>%
    arrange(desc(media)) %>%
    pull(Sample)
)

# 5. Plot
ggplot(df_comparacao, aes(x = Sample, y = Percent_Removed, fill = Threshold)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Comparação do % de Reads Removidos por Amostra",
    subtitle = "Thresholds: 0.01 vs 0.05",
    x = "Amostra",
    y = "% de Reads Removidos",
    fill = "Threshold"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )




