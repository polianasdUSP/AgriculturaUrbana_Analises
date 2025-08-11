
  
# Carregar as bibliotecas necessárias
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggcorrplot)
library(pheatmap)
library(grid)
library(patchwork)
library(ggpubr)
library(ggpmisc)
library(qiime2R)
library(tibble)


#========= Importar dados Alpha Diversity =====#

#=====> Shannon


# Importar o arquivo .qza
shannon <- read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/qiime_decontam_jun2025/core-metrics-results/shannon_vector.qza")

# Transformar em data frame com Sample.id e valor de Shannon
shannon_df <- shannon$data %>%
  rownames_to_column(var = "Sample.id") 




#====> faith
faith_pd <- read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/qiime_decontam_jun2025/core-metrics-results/faith_pd_vector.qza")


#verificar
faith_pd$data

# Transformar em data.frame com nomes corretos
# Criar data.frame com as colunas corretas
faith_pd <- faith_pd$data %>%
  select(V1, V2) %>%
  rename(Sample.id = V1, faith_pd = V2)


#===> Evenness


# Importar o arquivo .qza
evenness <- read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/qiime_decontam_jun2025/core-metrics-results/evenness_vector.qza")

# Transformar os rownames em coluna e renomear a coluna de valor
evenness_df <- evenness$data %>%
  rownames_to_column(var = "Sample.id") 

#=====> Observed



# Importar o arquivo .qza
observed_features <- read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/qiime_decontam_jun2025/core-metrics-results/observed_features_vector.qza")

# Transformar rownames em coluna e renomear a coluna de valor
observed_df <- observed_features$data %>%
  rownames_to_column(var = "Sample.id")

#======> Chao1


# Importar o arquivo .qza
chao1 <- read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/qiime_decontam_jun2025/core-metrics-results/chao1_vector.qza")

# Transformar rownames em coluna e renomear
chao1_df <- chao1$data %>%
  rownames_to_column(var = "Sample.id")


#======> Simpson


# Importar o arquivo .qza
simpson <- read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/qiime_decontam_jun2025/core-metrics-results/simpson_vector.qza")

# Transformar rownames em coluna e renomear
simpson_df <- simpson$data %>%
  rownames_to_column(var = "Sample.id")

#======> Alpha.all

library(dplyr)
library(purrr)

# Lista com todos os dataframes
alpha_list <- list(simpson_df, shannon_df, chao1_df, observed_df, evenness_df, faith_pd)

# Juntar todos pelo Sample.id
alpha.all <- reduce(alpha_list, left_join, by = "Sample.id")


#=======================================#
#    Juntas metadados e alpha
#=======================================#



metadados.saude.alpha <- metadados.saude.alpha %>%
  left_join(alpha.all, by = "Sample.id")


#write.csv(metadados.alpha.all, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados.alpha.all130.csv", row.names = FALSE)


#metadados_grupos_saude_binario <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados_grupos_saude_binario.csv", stringsAsFactors = FALSE)




# Juntar os dados com base em Sample.id
#metadados_grupos_saude_binario <- merge(
#  metadados.saude.alpha,
#  metadados_grupos_saude_binario[, c("Sample.id", "Health_Status")],
#  by = "Sample.id",
#  all.x = TRUE
#)

library(ggplot2)
library(ggpubr)

#Shannon X Health_Status

metadados_grupos_saude_binario$Health_Status <- factor(
  metadados_grupos_saude_binario$Health_Status,
  levels = c("Low Metabolic Risk", "Intermediate Metabolic Risk", "High Metabolic Risk")
)



library(ggplot2)
library(ggpubr)

ggplot(metadados_grupos_saude_binario, aes(x = Health_Status, y = shannon_entropy, fill = Health_Status)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "black", size = 1.5) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("#6BAED6", "#74C476", "#FB6A4A")) +
  labs(
    title = "Shannon Diversity by Metabolic Health Status",
    x = "Health Status",
    y = "Shannon Entropy"
  ) +
  coord_cartesian(ylim = c(3.5, 7.5)) +  # ajuste conforme sua distribuição
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(method = "kruskal.test", label.y = 7.4)

ggsave("shannon_vs_healthstatus.png", dpi = 600, width = 7, height = 5, units = "in")



#Faith_pd X Health Status

ggplot(metadados_grupos_saude_binario, aes(x = Health_Status, y = faith_pd, fill = Health_Status)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("#6BAED6", "#74C476", "#FB6A4A")) +
  labs(
    title = "Faith's Phylogenetic Diversity by Metabolic Health Status",
    x = "Health Status",
    y = "Faith PD"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(method = "kruskal.test", label.y = max(metadados_grupos_saude_binario$faith_pd, na.rm = TRUE) + 1)



#Chao1 X Health Status

ggplot(metadados_grupos_saude_binario, aes(x = Health_Status, y = chao1, fill = Health_Status)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("#6BAED6", "#74C476", "#FB6A4A")) +
  labs(
    title = "Chao1 Richness by Metabolic Health Status",
    x = "Health Status",
    y = "Chao1 Index"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(method = "kruskal.test", label.y = max(metadados_grupos_saude_binario$chao1, na.rm = TRUE) + 10)

#Simpson X Health Score
ggplot(metadados_grupos_saude_binario, aes(x = Health_Status, y = simpson, fill = Health_Status)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +  # oculta outliers para destacar o box
  geom_jitter(width = 0.2, alpha = 0.6, color = "black", size = 1.5) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("#6BAED6", "#74C476", "#FB6A4A")) +
  labs(
    title = "Simpson Index by Metabolic Health Status",
    x = "Health Status",
    y = "Simpson Index"
  ) +
  coord_cartesian(ylim = c(0.85, 1)) +  # foca na faixa útil
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(method = "kruskal.test", label.y = 0.995)


library(ggplot2)
library(ggpubr)

#pielou_evenness x health status

library(ggplot2)
library(ggpubr)

ggplot(metadados_grupos_saude_binario, aes(x = Health_Status, y = pielou_evenness, fill = Health_Status)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "black", size = 1.5) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("#6BAED6", "#74C476", "#FB6A4A")) +
  labs(
    title = "Pielou's Evenness by Metabolic Health Status",
    x = "Health Status",
    y = "Pielou's Evenness"
  ) +
  coord_cartesian(ylim = c(0.6, 1)) +  # ajuste conforme sua distribuição
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(method = "kruskal.test", label.y = 0.98)




#observed x health status

library(ggplot2)
library(ggpubr)

ggplot(metadados_grupos_saude_binario, aes(x = Health_Status, y = observed_features, fill = Health_Status)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "black", size = 1.5) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("#6BAED6", "#74C476", "#FB6A4A")) +
  labs(
    title = "Observed Features by Metabolic Health Status",
    x = "Health Status",
    y = "Observed Richness"
  ) +
  coord_cartesian(ylim = c(min(metadados_grupos_saude_binario$observed_features, na.rm = TRUE) - 10,
                           max(metadados_grupos_saude_binario$observed_features, na.rm = TRUE) + 10)) +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(method = "kruskal.test", label.y = max(metadados_grupos_saude_binario$observed_features, na.rm = TRUE) + 5)



#=======   Loop pra salvar todos de uma vez ======#

library(ggplot2)
library(ggpubr)

# Lista dos índices de diversidade alfa
indices_alpha <- c("shannon_entropy", "chao1", "faith_pd", "simpson", "pielou_evenness", "observed_features")

# Cores para os grupos
cores_grupos <- c("#6BAED6", "#74C476", "#FB6A4A")

# Caminho de destino
caminho <- "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/"

# Criar pasta se não existir
if (!dir.exists(caminho)) dir.create(caminho, recursive = TRUE)

# Loop para gerar e salvar os gráficos
for (indice in indices_alpha) {
  
  # Extraindo os valores para o eixo y
  y_min <- min(metadados_grupos_saude_binario[[indice]], na.rm = TRUE)
  y_max <- max(metadados_grupos_saude_binario[[indice]], na.rm = TRUE)
  margem <- (y_max - y_min) * 0.15
  label_y <- y_max + margem * 0.5
  
  p <- ggplot(metadados_grupos_saude_binario, aes_string(x = "Health_Status", y = indice, fill = "Health_Status")) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, color = "black", size = 1.5) +
    stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") +
    scale_fill_manual(values = cores_grupos) +
    labs(
      title = paste0(indice, " by Metabolic Health Status"),
      x = "Health Status",
      y = indice
    ) +
    coord_cartesian(ylim = c(y_min - margem * 0.3, y_max + margem)) +
    theme_minimal() +
    theme(legend.position = "none") +
    stat_compare_means(method = "kruskal.test", label.y = label_y)
  
  # Salvar como PNG
  ggsave(
    filename = paste0(caminho, indice, "_boxplot.png"),
    plot = p,
    width = 7, height = 5, dpi = 300
  )
}


#========> Health Status com parametros de saude

library(ggplot2)
library(ggpubr)

# Caminho para salvar os gráficos
caminho <- "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Box_plots/"

# Criar pasta se não existir
if (!dir.exists(caminho)) dir.create(caminho, recursive = TRUE)

# Variáveis de saúde (excluindo alpha diversidade + variáveis não clínicas)
variaveis_saude <- setdiff(colnames(metadados_grupos_saude_binario), 
                           c("Sample.id", "Sex", "geometry", "Health_Status",
                             "shannon_entropy", "pielou_evenness", "faith_pd",
                             "observed_features", "simpson", "chao1"))

# Loop para criar e salvar os boxplots
for (var in variaveis_saude) {
  
  y_min <- min(metadados_grupos_saude_binario[[var]], na.rm = TRUE)
  y_max <- max(metadados_grupos_saude_binario[[var]], na.rm = TRUE)
  margem <- (y_max - y_min) * 0.15
  label_y <- y_max + margem * 0.5
  
  p <- ggplot(metadados_grupos_saude_binario, aes_string(x = "Health_Status", y = var, fill = "Health_Status")) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, color = "black", size = 1.5) +
    stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") +
    scale_fill_manual(values = c("#6BAED6", "#74C476", "#FB6A4A")) +
    labs(
      title = paste0(var, " by Metabolic Health Status"),
      x = "Health Status",
      y = var
    ) +
    coord_cartesian(ylim = c(y_min - margem * 0.3, y_max + margem)) +
    theme_minimal() +
    theme(legend.position = "none") +
    stat_compare_means(method = "kruskal.test", label.y = label_y)
  
  ggsave(
    filename = paste0(caminho, var, "_boxplot.png"),
    plot = p,
    width = 7, height = 5, dpi = 300
  )
}



#=================================================#
#            Health Status x Diet                #
#=================================================#


#Juntar os dados com base em Sample.id

metadados_grupos_saude_binario_diet <- merge(
  metadados.dieta.all,
  metadados_grupos_saude_binario[, c("Sample.id", "Health_Status")],
  by = "Sample.id",
  all.x = TRUE
)


library(ggplot2)
library(ggpubr)

# Caminho para salvar os gráficos
caminho <- "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Box_plots/Boxplot_dieta/"

# Criar pasta se não existir
if (!dir.exists(caminho)) dir.create(caminho, recursive = TRUE)

# Variáveis dietéticas (exclui ID e Health_Status)
variaveis_dieta <- setdiff(colnames(metadados_grupos_saude_binario_diet), c("Sample.id", "Health_Status"))

# Loop para criar e salvar os boxplots
for (var in variaveis_dieta) {
  
  y_min <- min(metadados_grupos_saude_binario_diet[[var]], na.rm = TRUE)
  y_max <- max(metadados_grupos_saude_binario_diet[[var]], na.rm = TRUE)
  margem <- (y_max - y_min) * 0.15
  label_y <- y_max + margem * 0.5
  
  p <- ggplot(metadados_grupos_saude_binario_diet, aes_string(x = "Health_Status", y = var, fill = "Health_Status")) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, color = "black", size = 1.5) +
    stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") +
    scale_fill_manual(values = c("#6BAED6", "#74C476", "#FB6A4A")) +
    labs(
      title = paste0(var, " by Metabolic Health Status"),
      x = "Health Status",
      y = var
    ) +
    coord_cartesian(ylim = c(y_min - margem * 0.3, y_max + margem)) +
    theme_minimal() +
    theme(legend.position = "none") +
    stat_compare_means(method = "kruskal.test", label.y = label_y)
  
  # Salvar como PNG
  ggsave(
    filename = paste0(caminho, var, "_boxplot.png"),
    plot = p,
    width = 7, height = 5, dpi = 300
  )
}




#============================================#
#    r    Alpha e Parasitologico
#============================================#


indices_alpha <- c("shannon_entropy", "chao1", "observed_features", 
                   "pielou_evenness", "faith_pd", "simpson")


# Criar uma lista para armazenar os resultados
resultado_tests <- lapply(indices_alpha, function(indice) {
  teste <- wilcox.test(as.formula(paste(indice, "~ Parasitas")), data = metadados.all.filtrado)
  data.frame(
    Index = indice,
    W = teste$statistic,
    p_value = teste$p.value,
    row.names = NULL
  )
})

# Combinar em um único data.frame
resultado_tests_df <- do.call(rbind, resultado_tests)
resultado_tests_df$p_value_rounded <- round(resultado_tests_df$p_value, 4)

resultado_tests_df


library(ggplot2)
library(rlang)

plots_alpha <- list()

for (i in 1:nrow(resultado_tests_df)) {
  indice <- resultado_tests_df$Index[i]
  p_val <- resultado_tests_df$p_value_rounded[i]
  
  p <- ggplot(metadados.all.filtrado, aes(x = Parasitas, y = !!sym(indice), fill = Parasitas)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
    labs(
      title = paste0(indice, " (p = ", p_val, ")"),
      x = "", y = indice
    ) +
    theme_minimal(base_size = 16) +
    theme(legend.position = "none")
  
  plots_alpha[[indice]] <- p
}


print(plots_alpha[["shannon_entropy"]])
print(plots_alpha[["chao1"]])
print(plots_alpha[["observed_features"]])
print(plots_alpha[["pielou_evenness"]])
print(plots_alpha[["faith_pd"]])
print(plots_alpha[["simpson"]])





# Caminho de destino
output_dir <- "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/Box_plots"

# Salvar cada plot como PNG em alta resolução
for (indice in names(plots_alpha)) {
  ggsave(
    filename = paste0(output_dir, "/Boxplot_", indice, "_Parasitas.png"),
    plot = plots_alpha[[indice]],
    width = 12, height = 10, units = "cm", dpi = 600
  )
}

