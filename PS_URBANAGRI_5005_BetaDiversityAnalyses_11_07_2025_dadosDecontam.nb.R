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

library(dplyr)

metadados_grupos_saude_binario <- metadados_grupos_saude_binario %>%
  mutate(Health_Status = case_when(
    Health_Status == "Saudável" ~ "Healthy",
    Health_Status == "Transição" ~ "Transition",
    Health_Status == "Doente" ~ "Unhealthy",
    TRUE ~ Health_Status  # mantém o valor original se não for nenhum dos acima
  ))


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
permanova_result_unwei <- adonis2(as.dist(unweighted_unifrac) ~ Health_Status, data = metadados_grupos_saude_binario)

# Ver resultado
print(permanova_result_unwei)

#Extrair valor-p e R² para colocar no gráfico
p_value <- permanova_result_unwei$`Pr(>F)`[1]
r2_value <- permanova_result_unwei$R2[1]

# Criar label formatado
permanova_label_unwei <- paste0("PERMANOVA: R² = ", round(r2_value, 3), 
                          ", p = ", format.pval(p_value, digits = 3, eps = .001))




uw1 <- ggplot(df_pcoa_unw, aes(x = PCoA1, y = PCoA2, color = Health_Status)) +
  geom_point(size = 8, alpha = 0.7) +  # aumentei o tamanho dos pontos
  labs(
    title = "PCoA - Unweighted UniFrac",
    x = paste0("PCoA1 (", round(variance_explained.unw[1], 2), "%)"),
    y = paste0("PCoA2 (", round(variance_explained.unw[2], 2), "%)")
  ) +
  annotate("text", 
           x = 0.3, y = 0.25,  # ajuste se necessário
           label = permanova_label_unwei, 
           hjust = 1, size = 10) +  # aumentei o tamanho da fonte da anotação
  theme_minimal(base_size = 22) +  # aumentei o tamanho base da fonte do gráfico
  theme(legend.position = "right")

print(uw1)

# Salvar com alta resolução
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Pcoa_unweighted_Decontam_Health_Status_labels.png",
       plot = uw1, width = 12, height = 8, dpi = 600)  # dpi 600 para publicação


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
cat("PCoA1:", variance_explained.wei[1], "%\n")
cat("PCoA2:", variance_explained.wei[2], "%\n")

# Rodar PERMANOVA
permanova_result_wei <- adonis2(as.dist(weighted_unifrac) ~ Health_Status, data = metadados_grupos_saude_binario)

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
  geom_point(size = 8, alpha = 0.7) +  # pontos maiores
  labs(
    title = "PCoA - Weighted UniFrac",
    x = paste0("PCoA1 (", round(variance_explained.wei[1], 2), "%)"),
    y = paste0("PCoA2 (", round(variance_explained.wei[2], 2), "%)")
  ) +
  annotate("text", 
           x = 0.3, y = 0.25,
           label = permanova_label, 
           hjust = 1, size = 10) +  # fonte maior da anotação
  theme_minimal(base_size = 22) +  # fonte geral maior
  theme(legend.position = "right")


# 5. Salvar o gráfico
# Salvar com alta resolução
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Pcoa_weighted_Decontam_Health_Status_labels.png",
       plot = w1, width = 12, height = 8, dpi = 600)  # dpi 600 para publicação

library(patchwork)



# Combine os gráficos lado a lado
p_combined <- uw1 + w1 + plot_layout(ncol = 2, widths = c(1, 1))

# Salve com proporções mais realistas (ex: A4 em cm: 18x10)
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/PCoA_unweighted_weighted_side_by_side_corrigido.png",
       plot = p_combined, width = 18, height = 10, units = "cm", dpi = 600)

#====== Beta Diversity com Bray curtis =========#

bray_curtis_raw <- as.matrix(read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/qiime_decontam_jun2025/core-metrics-results/bray_curtis_distance_matrix.qza")$data)

# Limpar bray_curtis
bray_curtis <- bray_curtis_raw[!(rownames(bray_curtis_raw) %in% ids_remover),
                                             !(colnames(bray_curtis_raw) %in% ids_remover)]

# Passo 1: Garantir que os metadados tenham os mesmos IDs da matriz
# (apenas os que sobraram após limpeza)
ids_validos <- rownames(bray_curtis)
metadados_filtrados <- metadados_binarios_puros[metadados_binarios_puros$Sample.id %in% ids_validos, ]

# Passo 2: Garantir que a ordem das linhas do metadata bata com a matriz
metadados_filtrados <- metadados_filtrados[match(ids_validos, metadados_filtrados$Sample.id), ]

# Verificar se está tudo certo
stopifnot(identical(rownames(bray_curtis), metadados_filtrados$Sample.id))

# Passo 3: Rodar PERMANOVA
set.seed(123)  # Para reprodutibilidade
resultado_permanova_bray <- adonis2(bray_curtis ~ Health_Status, data = metadados_grupos_saude_binario)

# Mostrar resultado
print(resultado_permanova_bray)

#=======> Jaccard com Health_Status

jaccard_raw <- as.matrix(read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/qiime_decontam_jun2025/core-metrics-results/jaccard_distance_matrix.qza")$data)

# Unweighted UniFrac
jaccard <- jaccard_raw[!(rownames(jaccard_raw) %in% ids_remover),
                               !(colnames(jaccard_raw) %in% ids_remover)]




set.seed(123)  # Para reprodutibilidade
resultado_permanova_jaccard <- adonis2(jaccard ~ Health_Status, data = metadados_grupos_saude_binario)

# Mostrar resultado
print(resultado_permanova_jaccard)


#========================================#
#       Beta com Parasitologico
#========================================#


metadados.all.filtrado$Parasitas <- ifelse(metadados.all.filtrado$Parasitological == "Sim", "Positivo", "Negativo")


library(vegan)

# PERMANOVA com Unweighted UniFrac
adonis_unw <- adonis2(unweighted_unifrac ~ Parasitas, data = metadados.all.filtrado)
print(adonis_unw)

# PERMANOVA com Weighted UniFrac
adonis_w <- adonis2(weighted_unifrac ~ Parasitas, data = metadados.all.filtrado)
print(adonis_w)

library(dplyr)

df_pcoa_wei_parasitas <- pcoa_result.wei$vectors %>%
  as.data.frame() %>%
  tibble::rownames_to_column("SampleID") %>%
  left_join(metadados.all.filtrado[, c("Sample.id", "Parasitas")], by = c("SampleID" = "Sample.id"))

df_pcoa_unw_parasitas <- pcoa_result.unw$vectors %>%
  as.data.frame() %>%
  tibble::rownames_to_column("SampleID") %>%
  left_join(metadados.all.filtrado[, c("Sample.id", "Parasitas")], by = c("SampleID" = "Sample.id"))



library(ggplot2)

# Etiqueta da PERMANOVA
permanova_label_w <- "PERMANOVA\nR² = 0.013 | p = 0.099"

pcoa_w <- ggplot(df_pcoa_wei_parasitas, aes(x = Axis.1, y = Axis.2, color = Parasitas)) +
  geom_point(size = 5, alpha = 0.7) +  # pontos grandes e translúcidos
  labs(
    title = "PCoA - Weighted UniFrac",
    x = paste0("PCoA1 (", round(variance_explained.wei[1], 2), "%)"),
    y = paste0("PCoA2 (", round(variance_explained.wei[2], 2), "%)")
  ) +
  annotate("text", x = 0.3, y = 0.25, label = permanova_label_w, size = 6, hjust = 1) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "right")

print (pcoa_w)


permanova_label_unw <- "PERMANOVA\nR² = 0.017 | p = 0.005"  # Substitua após rodar o teste

pcoa_unw <- ggplot(df_pcoa_unw_parasitas, aes(x = Axis.1, y = Axis.2, color = Parasitas)) +
  geom_point(size = 5, alpha = 0.7) +
  labs(
    title = "PCoA - Unweighted UniFrac",
    x = paste0("PCoA1 (", round(variance_explained.unw[1], 2), "%)"),
    y = paste0("PCoA2 (", round(variance_explained.unw[2], 2), "%)")
  ) +
  annotate("text", x = 0.3, y = 0.25, label = permanova_label_unw, size = 6, hjust = 1) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "right")


#======= Clusters ============#

library(factoextra)

# Supondo que você já tenha isso assim:
coords_unw <- pcoa_result.unw$vectors[, 1:2]
coords_w <- pcoa_result.wei$vectors[, 1:2]

#2. K-means clustering com k = 2
set.seed(123)
kmeans_unw <- kmeans(coords_unw, centers = 2, nstart = 25)
kmeans_w <- kmeans(coords_w, centers = 2, nstart = 25)




# Unweighted
fviz_cluster(kmeans_unw, data = coords_unw,
             geom = "point", main = "K-means Clustering (Unweighted)")

# Weighted
fviz_cluster(kmeans_w, data = coords_w,
             geom = "point", main = "K-means Clustering (Weighted)")


library(cluster)

# Unweighted
sil_unw <- silhouette(kmeans_unw$cluster, dist(coords_unw))
plot(sil_unw, main = "Silhouette Plot (Unweighted)")

# Weighted
sil_w <- silhouette(kmeans_w$cluster, dist(coords_w))
plot(sil_w, main = "Silhouette Plot (Weighted)")

# Supondo que você tenha o objeto metadata
metadados_grupos_saude_binario$Cluster_unw <- factor(kmeans_unw$cluster)


# Garantir que os nomes estejam sincronizados
clusters_filtrados <- kmeans_w$cluster[names(kmeans_w$cluster) %in% metadados_grupos_saude_binario$Sample.id]

# Ordenar os clusters para que correspondam à ordem de Sample.id
clusters_ordenados <- clusters_filtrados[match(metadados_grupos_saude_binario$Sample.id, names(clusters_filtrados))]

# Atribuir ao metadados
metadados_grupos_saude_binario$Cluster_w <- factor(clusters_ordenados)


metadados_grupos_saude_binario$Cluster_w <- factor(kmeans_w$cluster)

# Use a matriz de distância original (por exemplo, dist_unw ou dist_w)
adonis2(dist_unw ~ Cluster_unw, data = metadados.all.filtrado)
adonis2(dist_w ~ Cluster_w, data = metadados.all.filtrado)


#Cross-tabulação (tabela de contingência): quantos indivíduos caíram no mesmo grupo nos dois métodos. Se os clusters forem bem diferentes, os números nas diagonais serão baixos.

table(metadados.all.filtrado$Cluster_unw, metadados.all.filtrado$Cluster_w)

# Índice de concordância (Adjusted Rand Index): Usa o pacote mclust para calcular o ARI (quanto mais perto de 1, mais semelhantes os agrupamentos).

library(mclust)
adjustedRandIndex(metadados.all.filtrado$Cluster_unw, metadados.all.filtrado$Cluster_w)

#============================================================================#
#                               Resultado com Cluster                       #      
#============================================================================#

library(vegan)

# Unweighted UniFrac + Clusters gerados via K-means
adonis2(unweighted_unifrac ~ Cluster_unw, data = metadados_grupos_saude_binario)

# Weighted UniFrac + Clusters gerados via K-means
adonis2(weighted_unifrac ~ Cluster_w, data = metadados_grupos_saude_binario)

#============================================================================#
#                        Weighted e Variaveis                                #      
#============================================================================#
#Variáveis de saúde metabólica e inflamatória:

vars_saude <- c(
  "BMI", "Waist_Hip_Ratio", "Systolic", "Diastolic",
  "HbA1c", "Glucose", "Insulin", "CRP",
  "Triglycerides", "Cholesterol", "LDL", "HDL", "VLDL",
  "TGO", "TGP", "GGT",
  "IL6", "IL10", "IL17A", "TNF", "IFNGamma", "IL2", "IL4", "Metabolic_Score"
)

#Código para comparar as variáveis entre os clusters (Wilcoxon):
# Loop para testar associação entre cada variável e o Cluster_w

results_saude <- lapply(vars_saude, function(var) {
  formula <- as.formula(paste(var, "~ Cluster_w"))
  test <- tryCatch(wilcox.test(formula, data = metadados_grupos_saude_binario),
                   error = function(e) NA)
  data.frame(Variable = var,
             p_value = ifelse(is.na(test), NA, round(test$p.value, 4)))
})

# Resultado final em tabela
results_saude_df <- do.call(rbind, results_saude)
results_saude_df <- results_saude_df[order(results_saude_df$p_value), ]
results_saude_df

#========> Boxplot de Insulin por Cluster_w
# Instale se necessário:
# install.packages("ggpubr")

library(ggplot2)
library(dplyr)
library(ggpubr)

# Garantir que Cluster_w seja fator
metadados.all.filtrado$Cluster_w <- as.factor(metadados.all.filtrado$Cluster_w)

# Calcular n por grupo
insulin_n <- metadados.all.filtrado %>%
  group_by(Cluster_w) %>%
  summarise(n = sum(!is.na(Insulin)))

# Criar o gráfico com valor de p (teste de Wilcoxon)
ggplot(metadados.all.filtrado, aes(x = Cluster_w, y = Insulin, fill = Cluster_w)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label.y = max(metadados.all.filtrado$Insulin, na.rm = TRUE) + 5) +
  scale_x_discrete(labels = paste0("Cluster ", insulin_n$Cluster_w, "\n(n=", insulin_n$n, ")")) +
  labs(title = "Insulin by Microbiota Cluster (Weighted UniFrac)",
       x = "Cluster",
       y = "Insulin (μU/mL)") +
  theme_minimal()


#========> Boxplot de TNF por cluster_w

library(ggplot2)
library(dplyr)

# Garantir que Cluster_w seja fator
metadados.all.filtrado$Cluster_w <- as.factor(metadados.all.filtrado$Cluster_w)

# Calcular n por grupo para TNF
tnf_n <- metadados.all.filtrado %>%
  group_by(Cluster_w) %>%
  summarise(n = sum(!is.na(TNF)))





# Garantir que Cluster_w seja fator
metadados.all.filtrado$Cluster_w <- as.factor(metadados.all.filtrado$Cluster_w)

# Criar gráfico com valor de p do teste de Wilcoxon (ou t.test)
ggplot(metadados.all.filtrado, aes(x = Cluster_w, y = TNF, fill = Cluster_w)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label.y = max(metadados.all.filtrado$TNF, na.rm = TRUE) + 2) +
  scale_x_discrete(labels = c("Cluster 1", "Cluster 2")) +
  labs(title = "TNF by Microbiota Cluster (Weighted UniFrac)",
       x = "Cluster",
       y = "TNF (pg/mL)") +
  theme_minimal()


#========> Boxplot de Metabolic_Score por cluster_w
# Garantir que Cluster_w seja fator
metadados_grupos_saude_binario$Cluster_w <- as.factor(metadados_grupos_saude_binario$Cluster_w)

# Calcular n por grupo
insulin_n <- metadados_grupos_saude_binario %>%
  group_by(Cluster_w) %>%
  summarise(n = sum(!is.na(Metabolic_Score)))

# Criar o gráfico com valor de p (teste de Wilcoxon)
ggplot(metadados_grupos_saude_binario, aes(x = Cluster_w, y = Metabolic_Score, fill = Cluster_w)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label.y = max(metadados_grupos_saude_binario$Metabolic_Score, na.rm = TRUE) + 5) +
  scale_x_discrete(labels = paste0("Cluster ", insulin_n$Cluster_w, "\n(n=", insulin_n$n, ")")) +
  labs(title = "Metabolic_Score by Microbiota Cluster (Weighted UniFrac)",
       x = "Cluster",
       y = "Insulin (μU/mL)") +
  theme_minimal()


#====== Clusters vs Dieta


# Selecionar apenas variáveis numéricas de dieta (colunas 2 a 29)
vars_dieta <- colnames(metadados.dieta.all)[c(2:29, 38:41)]


# Rodar o teste de Wilcoxon para cada variável vs Cluster_w
results_dieta <- lapply(vars_dieta, function(var) {
  formula <- as.formula(paste(var, "~ Cluster_w"))
  test <- tryCatch(wilcox.test(formula, data = metadados.dieta.all),
                   error = function(e) NA)
  data.frame(
    Variable = var,
    p_value = ifelse(is.na(test), NA, round(test$p.value, 4))
  )
})

# Juntar os resultados e ordenar por p-valor
results_dieta_df <- do.call(rbind, results_dieta) %>%
  arrange(p_value)

# Visualizar
results_dieta_df


#==============> Boxplot

library(ggplot2)
library(dplyr)
library(ggpubr)

# Calcular n por grupo
n_calcio <- metadados.dieta.all %>%
  group_by(Cluster_w) %>%
  summarise(n = sum(!is.na(calcio_mg)))

n_fibra <- metadados.dieta.all %>%
  group_by(Cluster_w) %>%
  summarise(n = sum(!is.na(fibra_alimentar_g)))

# Boxplot para Cálcio
plot_calcio <- ggplot(metadados.dieta.all, aes(x = factor(Cluster_w), y = calcio_mg, fill = factor(Cluster_w))) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_x_discrete(labels = paste0("Cluster ", n_calcio$Cluster_w, "\n(n=", n_calcio$n, ")")) +
  labs(title = "Calcium (mg) by Microbiota Cluster (Weighted Unifrac)",
       x = "Cluster", y = "Residual Calcium (mg)") +
  theme_minimal()

# Boxplot para Fibra
plot_fibra <- ggplot(metadados.dieta.all, aes(x = factor(Cluster_w), y = fibra_alimentar_g, fill = factor(Cluster_w))) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_x_discrete(labels = paste0("Cluster ", n_fibra$Cluster_w, "\n(n=", n_fibra$n, ")")) +
  labs(title = "Dietary Fiber (g) by Microbiota Cluster (Weighted Unifrac)",
       x = "Cluster", y = "Residual Fiber (g)") +
  theme_minimal()

# Mostrar os plots
plot_calcio
plot_fibra


#====================================================#
#         BHEI e NOVA
#=====================================================#

# BHEI
bhei <- read.csv2("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/Dieta/BHEI_R_scores.csv", stringsAsFactors = FALSE)


# NOVA
library(readxl)
library(readxl)

nova <- read_excel("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas/nova_wide_percentagens.xlsx")

metadados_grupos_saude_binario2 <- read_excel("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados_grupos_saude_binario.xlsx") 



#selecionar colunas nevessarias:
# Selecionar colunas relevantes de cada tabela
bhei_sel <- bhei[, c("IdVoluntario", "BHEI_R_Score_Total")]
nova_sel <- nova[, c("IdVoluntario", "Percentual_NOVA_group_1", "Percentual_NOVA_group_2", "Percentual_NOVA_group_3")]

#Renomear a coluna de ID em nova_sel para coincidir com metadados.dieta.all:
colnames(nova_sel)[1] <- "Sample.id"
colnames(bhei_sel)[1] <- "Sample.id"

#Mesclar tudo com metadados.dieta.all:
metadados.dieta.all <- metadados.dieta.all %>%
  left_join(bhei_sel, by = "Sample.id") %>%
  left_join(nova_sel, by = "Sample.id")


#=======================================================#
#     Colorir as PCoA com o que deu significante
#=======================================================#

# PCoA weighted
df_pcoa_wei <- left_join(df_pcoa_wei, metadados_grupos_saude_binario, by = c("SampleID" = "Sample.id"))


# PCoA unweighted
df_pcoa_unw <- left_join(df_pcoa_unw, metadados_grupos_saude_binario, by = c("SampleID" = "Sample.id") )


library(ggplot2)

ggplot(df_pcoa_wei, aes(x = PCoA1, y = PCoA2, color = TNF)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "D", direction = -1) +
  labs(title = "PCoA - Weighted UniFrac (cor: TNF)",
       x = "PCoA1", y = "PCoA2", color = "TNF") +
  theme_minimal()


ggplot(df_pcoa_wei, aes(x = PCoA1, y = PCoA2, color = Insulin)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "D", direction = -1) +
  labs(title = "PCoA - Weighted UniFrac (cor: Insulin)",
       x = "PCoA1", y = "PCoA2", color = "Insulin") +
  theme_minimal()

ggplot(df_pcoa_wei, aes(x = PCoA1, y = PCoA2, color = Metabolic_Score)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "D", direction = -1) +
  labs(title = "PCoA - Weighted UniFrac (cor: Metabolic Score)",
       x = "PCoA1", y = "PCoA2", color = "Metabolic Score") +
  theme_minimal()

# PCoA weighted
df_pcoa_wei <- left_join(df_pcoa_wei, metadados.dieta.all, by = c("SampleID" = "Sample.id"))


# PCoA unweighted
df_pcoa_unw <- left_join(df_pcoa_unw, metadados.dieta.all, by = c("SampleID" = "Sample.id") )


ggplot(df_pcoa_wei, aes(x = PCoA1, y = PCoA2, color = calcio_mg)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "D", direction = -1) +
  labs(title = "PCoA - Weighted UniFrac (calcium (mg))",
       x = "PCoA1", y = "PCoA2", color = "Calcium") +
  theme_minimal()


ggplot(df_pcoa_wei, aes(x = PCoA1, y = PCoA2, color = fibra_alimentar_g)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "D", direction = -1) +
  labs(title = "PCoA - Weighted UniFrac (Dietary Fiber)",
       x = "PCoA1", y = "PCoA2", color = "Dietary Fiber") +
  theme_minimal()
