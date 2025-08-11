
library(ape)

#=================





# Selecionar variáveis de interesse
vars_metabolicas2 <- c("Sample.id", "Glucose", "Triglycerides", "Cholesterol",
                       "LDL", "HDL", "HbA1c", "Systolic", "Diastolic", 
                       "CRP", "TGO", "TGP", "GGT")

# Separar a coluna de ID
sample_ids2 <- metadados.all.filtrado$Sample.id

# Padronizar apenas as variáveis numéricas / Remover NAs e padronizar (z-score)
df_scaled2 <- metadados.all.filtrado %>%
  select(-Sample.id) %>%
  select(all_of(vars_metabolicas2[-1])) %>%  # remove "Sample.id" do vetor
  scale(center = TRUE, scale = TRUE) %>%
  as.data.frame()

# Adicionar o Sample.id de volta
df_scaled2$Sample.id <- sample_ids2

#Código para reorganizar as colunas (Sample.id primeiro)
df_scaled2 <- df_scaled2 %>%
  relocate(Sample.id, .before = 1)


#Remover NAs
df_scaled_clean2 <- df_scaled2 %>% drop_na()

#Ou, se quiser ver quantas linhas foram removidas:
#n_inicial <- nrow(df_scaled)
#df_scaled_clean <- df_scaled %>% drop_na()
#n_final <- nrow(df_scaled_clean)
#cat("Linhas removidas:", n_inicial - n_final, "\n")


#================> Criar 3 grupos separados (super saudavel, intermediario e )

# 1. Rodar o k-means com 3 grupos

set.seed(123)  # Para reprodutibilidade

# Rodar k-means (removendo Sample.id antes)
kmeans_result2 <- kmeans(df_scaled_clean2 %>% select(-Sample.id), centers = 3)

# Ver quantos em cada grupo
table(kmeans_result2$cluster)

#Adicionar a coluna de grupo ao dataframe
# Adicionar o cluster ao df_scaled
df_scaled_clean2$HealthGroup <- as.factor(kmeans_result2$cluster)

## Juntar com metadados.all.filtrado pelo Sample.id
metadados_clusterizado2 <- metadados.all.filtrado %>%
  inner_join(df_scaled_clean2 %>% select(Sample.id, HealthGroup), by = "Sample.id")


# Observar médias por grupo
df_scaled_clean2 %>%
  group_by(HealthGroup) %>%
  summarise(across(where(is.numeric), mean))

# Depois de avaliar, renomear:
#metadados_clusterizado$HealthGroup_named <- recode_factor(
#  metadados_clusterizado$HealthGroup,
#  `1` = "Optimal Health",
#  `2` = "High Risk Metabolic Profile"

)

#============> PCoA

#fazer uma PCoA com seus dados de saúde:

# Remover Sample.id e a coluna de grupo antes de calcular distância
dados_pcoa2 <- df_scaled_clean2 %>%
  select(-Sample.id, -HealthGroup)

# Calcular matriz de distância
dist_matrix2 <- dist(dados_pcoa2, method = "euclidean")

# Rodar PCoA com cmdscale
pcoa_result2 <- cmdscale(dist_matrix2, k = 2, eig = TRUE)

# Criar dataframe com coordenadas
pcoa_df2 <- as.data.frame(pcoa_result2$points)
colnames(pcoa_df2) <- c("PCoA1", "PCoA2")




# Adicionar os grupos
pcoa_df2$Sample.id <- df_scaled_clean2$Sample.id
pcoa_df2$HealthGroup <- metadados_clusterizado2$HealthGroup

library(ggplot2)

p4 <- ggplot(pcoa_df2, aes(x = PCoA1, y = PCoA2, color = HealthGroup)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCoA Based on Metabolic Markers",
    x = "PCoA1",
    y = "PCoA2",
    color = "Health Group"
  ) +
  theme_minimal() +
  scale_color_manual(values = c(
    "1" = "#1b9e77",
    "2" = "#d95f02"
    
  ))



print(p4)


#=====> Permanova

library(vegan)
permanova_health_group2 <- adonis2(dist_matrix2 ~ metadados_clusterizado2$HealthGroup)

r2 <- round(permanova_health_group2$R2[1], 3)
pval <- format.pval(permanova_health_group2$`Pr(>F)`[1], digits = 3, eps = .001)


#!!!! ====== Este grafico abaixo ======= !!!!#

p5 <- ggplot(pcoa_df2, aes(x = PCoA1, y = PCoA2, color = HealthGroup)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCoA Based on Metabolic Markers (n=66)",
    x = "PCoA1",
    y = "PCoA2",
    color = "Health Group",
    caption = paste0("PERMANOVA: R² = ", r2, ", p = ", pval)
  ) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c(
    "1" = "#1b9e77",
    "2" = "#d95f02"
  ))

print(p5)

ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/PCoA_HealthGroups.png", plot = p5, width = 10, height = 6, dpi = 300)


#========> Colorir PCoA com todos os parametros de saúde, pra entender quais tem mais peso

library(dplyr)


#Joins da pcoa com os parametros
pcoa_merged2 <- pcoa_df2 %>%
  left_join(metadados.all.filtrado %>% 
              select(Sample.id, Age, HbA1c, INSULINA, Glucose, pulse, Cholesterol,
                     LDL, HDL, Triglycerides, TGO, TGP, GGT, BMI, PCR, BMI), 
            by = "Sample.id")


library(ggplot2)



#HbA1c

ggplot(pcoa_merged2, aes(x = PCoA1, y = PCoA2, color = HbA1c)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_gradient(low = "#d9f0a3", high = "#005a32") +
  labs(title = "PCoA colored by HbA1c",
       x = "PCoA1", y = "PCoA2", color = "HbA1c") +
  theme_minimal(base_size = 16)




#loop para todos os parametros

library(ggplot2)
library(scales)

# Lista de variáveis que você quer colorir (cada uma de cada vez)
variaveis2 <- c("Age", "HbA1c", "Glucose", "INSULINA", "pulse", "Cholesterol", 
                "LDL", "HDL", "Triglycerides", "TGO", "TGP", "GGT", 
                "BMI", "PCR", "BMI")



# Loop para gerar gráfico por variável
for (var in variaveis2) {
  p <- ggplot(pcoa_merged2, aes(x = PCoA1, y = PCoA2)) +
    geom_point(aes_string(color = var), size = 5, alpha = 0.9) +
    scale_color_viridis_c(option = "D", direction = -1) +  # cor sólida com viridis
    labs(
      title = paste("PCoA colored by", var),
      x = "PCoA1", y = "PCoA2", color = var
    ) +
    theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", size = 22),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16)
    )
  
  # Salvar o gráfico
  ggsave(
    filename = paste0("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/PCoA_plots/PCoA_plots_cluster/PCoA_", var, ".png"),
    plot = p,
    width = 10, height = 7, dpi = 300
  )
}


#====== Cluster binarios (sim/nao) ========#

#Criar metadados binários para fazer kmeans (distância euclideana) e pcoa

metadados_upset <- metadados.all.filtrado %>%
  mutate(
    High_Triglycerides = ifelse(!is.na(Triglycerides) & Triglycerides > 150, 1, 0),
    High_Fasting_Glucose = ifelse(!is.na(Glucose) & Glucose > 100, 1, 0),
    Obesity = ifelse(BMI >= 30, 1, 0),  # BMI tem 0 NA, sem ifelse extra
    Low_HDL = ifelse(!is.na(HDL) & HDL < 40, 1, 0),
    High_Blood_Pressure = ifelse(!is.na(Systolic) & !is.na(Diastolic) &
                                   (Systolic >= 130 | Diastolic >= 85), 1, 0)
  )



# índice Jaccard, por padrão no vegdist(), ignora co-ausências (linhas com todos os 0s). Usar a distância Gower, que trata zeros como ausência legítima

# Carregar pacotes necessários
library(dplyr)
library(cluster)   # para daisy()
library(ggplot2)

# 1. Selecionar variáveis de interesse
vars_upset <- c("Sample.id", "High_Triglycerides", "High_Fasting_Glucose", 
                "Obesity", "Low_HDL", "High_Blood_Pressure")

# 2. Criar dataframe com essas variáveis
df_upset <- metadados_upset %>%
  select(all_of(vars_upset))

# 3. Remover NAs
df_upset_clean <- df_upset %>% drop_na()

# 4. Rodar k-means clustering com 3 grupos (sem incluir Sample.id)
set.seed(123)
kmeans_result <- kmeans(df_upset_clean %>% select(-Sample.id), centers = 3)

# 5. Adicionar coluna com grupos ao dataframe
df_upset_clean$HealthGroup <- as.factor(kmeans_result$cluster)

# 6. Calcular distância Gower (preserva amostras com tudo 0)
dados_binarios <- df_upset_clean %>% select(-Sample.id, -HealthGroup)
dist_matrix <- daisy(dados_binarios, metric = "euclidean")

# 7. Rodar PCoA
pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)

# 8. Criar dataframe com coordenadas PCoA
pcoa_df <- as.data.frame(pcoa_result$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$Sample.id <- df_upset_clean$Sample.id
pcoa_df$HealthGroup <- df_upset_clean$HealthGroup

# 9. Calcular variância explicada por eixo
eig_vals <- pcoa_result$eig
var_exp <- round(100 * eig_vals[1:2] / sum(eig_vals), 1)

# 10. Plotar gráfico PCoA
e <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = HealthGroup)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCoA Groups by Clinical Reference Ranges",
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    color = "Health Group"
  ) +
  theme_minimal(base_size = 16) +
  scale_color_manual(values = c("1" = "#1b9e77", 
                                "2" = "#d95f02", 
                                "3" = "#7570b3"))

print(e)



# Salvar o gráfico

ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/PCoA_cluster", plot = e, width = 10, height = 6, dpi = 300)




#media de risco por grupo

df_upset_clean %>%
  group_by(HealthGroup) %>%
  summarise(across(where(is.numeric), mean))


#Ver quantos pontos estao em coordenadas iguais

pcoa_df %>%
  group_by(PCoA1, PCoA2) %>%
  summarise(count = n()) %>%
  filter(count > 1)

#Pontos ficaram sobrepostos, entao usar jitter pra visualizar melhor:

e1 <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = HealthGroup)) +
  geom_jitter(width = 0.02, height = 0.02, size = 3, alpha = 0.7) +
  labs(
    title = "PCoA dos Grupos por Critérios Clínicos",
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    color = "Grupo de Saúde"
  ) +
  theme_minimal(base_size = 16) +
  scale_color_manual(values = c(
    "1" = "#1b9e77",
    "2" = "#d95f02",
    "3" = "#7570b3"
  ))


print(e1)
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/PCoA_cluster_jitter", plot = e1, width = 10, height = 6, dpi = 300)


#====>  Variacao explicada 

pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)
eig_vals <- pcoa_result$eig
var_exp <- round(100 * eig_vals / sum(eig_vals), 1)

var_exp[1]


cat("Variação explicada por PCoA1:", var_exp[1], "%\n")
cat("Variação explicada por PCoA2:", var_exp[2], "%\n")



#Joins da pcoa com os parametros
pcoa_merged3 <- pcoa_df %>%
  left_join(metadados.all.filtrado %>% 
              select(Sample.id, Age, HbA1c, INSULINA, Glucose, pulse, HDL, Triglycerides, TGO, TGP, GGT, BMI, PCR, BMI), 
            by = "Sample.id")


library(ggplot2)


# Lista de variáveis
variaveis3 <- c("Age", "Glucose", "pulse", "HDL",
                "Triglycerides", "BMI")

# Loop para gerar gráfico por variável
for (var in variaveis3) {
  p <- ggplot(pcoa_merged3, aes(x = PCoA1, y = PCoA2)) +
    geom_jitter(aes_string(color = var), size = 5, alpha = 0.9, width = 0.1, height = 0.1) +
    scale_color_viridis_c(option = "D", direction = -1) +
    labs(
      title = paste("PCoA colored by", var),
      x = "PCoA1", y = "PCoA2", color = var
    ) +
    theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", size = 22),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16)
    )
  
  ggsave(
    filename = paste0("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/PCoA_plots/PCoA_plots_cluster/PCoA_Upset_", var, ".png"),
    plot = p,
    width = 10, height = 7, dpi = 300
  )
}


#==================== Clusters binários = Completo ===========#

#Criar metadados binários para fazer kmeans (distância euclideana) e pcoa

metadados_grupos_saude_binario <- metadados.all.filtrado %>%
  mutate(
    # Pressão arterial
    High_Blood_Pressure = ifelse(!is.na(Systolic) & !is.na(Diastolic) &
                                   (Systolic >= 130 | Diastolic >= 85), 1, 0),
    
    # Glucose de jejum
    High_Fasting_Glucose = case_when(
      is.na(Glucose) ~ NA_integer_,
      Glucose <= 100 ~ 0,                 # normal
      Glucose > 100 & Glucose < 126 ~ 1,  # alterada (pré-diabetes)
      Glucose >= 126 ~ 2                  # muito alta (diabetes)
    ),
    
    # HbA1c em porcentagem
    HbA1c_cat = case_when(
      is.na(HbA1c) ~ NA_integer_,
      HbA1c < 6.0 ~ 0,
      HbA1c >= 6.0 & HbA1c < 6.5 ~ 1,
      HbA1c >= 6.5 ~ 2
    ),
    
    
    # BMI categorizado por faixa etária
    BMI_cat = case_when(
      is.na(BMI) | is.na(Age) ~ NA_integer_,
      Age >= 60 & BMI < 22 ~ 1,
      Age >= 60 & BMI >= 22 & BMI < 27 ~ 0,
      Age >= 60 & BMI >= 27 & BMI < 30 ~ 1,
      Age >= 60 & BMI >= 30 & BMI < 35 ~ 2,
      Age >= 60 & BMI >= 35 ~ 3,
      Age < 60 & BMI < 18.5 ~ 1,
      Age < 60 & BMI >= 18.5 & BMI < 25 ~ 0,
      Age < 60 & BMI >= 25 & BMI < 30 ~ 1,
      Age < 60 & BMI >= 30 & BMI < 35 ~ 2,
      Age < 60 & BMI >= 35 & BMI < 40 ~ 3,
      Age < 60 & BMI >= 40 ~ 4
    ),
    
    # Triglicerídeos
    TRIG_cat = case_when(
      is.na(Triglycerides) ~ NA_integer_,
      Triglycerides < 150 ~ 0,
      Triglycerides >= 150 & Triglycerides < 200 ~ 1,
      Triglycerides >= 200 & Triglycerides < 500 ~ 2,
      Triglycerides >= 500 ~ 3
    ),
    
    # Cholesterol total
    COL_cat = case_when(
      is.na(Cholesterol) ~ NA_integer_,
      Cholesterol <= 190 ~ 0,
      Cholesterol > 190 & Cholesterol <= 240 ~ 1,
      Cholesterol > 240 ~ 2
    ),
    
    
    # HDL
    HDL_cat = case_when(
      is.na(HDL) ~ NA_integer_,
      HDL < 40 ~ 1,        # risco
      HDL >= 40 ~ 0 ),       # normal ou bom
    
    
    
    # LDL
    LDL_cat = case_when(
      is.na(LDL) ~ NA_integer_,
      LDL < 100 ~ 0,
      LDL >= 100 & LDL < 130 ~ 0,
      LDL >= 130 & LDL < 160 ~ 1,
      LDL >= 160 & LDL < 190 ~ 2,
      LDL >= 190 ~ 3
    ),
    
    # GGT por sexo
    GGT_cat = case_when(
      is.na(GGT) | is.na(Sex) ~ NA_integer_,
      Sex == "Masculino" & GGT <= 48 ~ 0,
      Sex == "Masculino" & GGT > 48 ~ 1,
      Sex == "Feminino" & GGT <= 32 ~ 0,
      Sex == "Feminino" & GGT > 32 ~ 1
    )
  )





# Vamos usar distancia euclideana


# Carregar pacotes necessários
library(dplyr)
library(cluster)   # para daisy()
library(ggplot2)

# 1. Selecionar variáveis de interesse
vars_grupo_bin <- c("Sample.id", "High_Blood_Pressure", "High_Fasting_Glucose",               
                    "HbA1c_cat",  "BMI_cat","TRIG_cat", "COL_cat","HDL_cat","LDL_cat",                           
                    "GGT_cat")

# 2. Criar dataframe com essas variáveis
df_grupo_bin <- metadados_grupos_saude_binario %>%
  select(all_of(vars_grupo_bin))

# 3. Remover NAs
df_grupo_bin_clean <- df_grupo_bin %>% drop_na()

# 4. Rodar k-means clustering com 3 grupos (sem incluir Sample.id)
set.seed(123)
kmeans_result_grupo_bin <- kmeans(df_grupo_bin_clean %>% select(-Sample.id), centers = 3)

# 5. Adicionar coluna com grupos ao dataframe
df_grupo_bin_clean$HealthGroup <- as.factor(kmeans_result_grupo_bin$cluster)

# 6. Calcular distância Gower (preserva amostras com tudo 0)
dados_binarios_grupo_bin <- df_grupo_bin_clean %>% select(-Sample.id, -HealthGroup)
dist_matrix_grupo_bin <- daisy(dados_binarios_grupo_bin, metric = "euclidean")

# 7. Rodar PCoA
pcoa_result_grupo_bin <- cmdscale(dist_matrix_grupo_bin, k = 2, eig = TRUE)

# 8. Criar dataframe com coordenadas PCoA
pcoa_df_grupo_bin <- as.data.frame(pcoa_result_grupo_bin$points)
colnames(pcoa_df_grupo_bin) <- c("PCoA1", "PCoA2")
pcoa_df_grupo_bin$Sample.id <- df_grupo_bin_clean$Sample.id
pcoa_df_grupo_bin$HealthGroup <- df_grupo_bin_clean$HealthGroup

# 9. Calcular variância explicada por eixo
eig_vals_grupo_bin <- pcoa_result_grupo_bin$eig
var_exp_grupo_bin <- round(100 * eig_vals_grupo_bin[1:2] / sum(eig_vals_grupo_bin), 1)

# 10. Plotar gráfico PCoA
g.bin <- ggplot(pcoa_df_grupo_bin, aes(x = PCoA1, y = PCoA2, color = HealthGroup)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCoA Health Groups by Clinical Reference Ranges",
    x = paste0("PCoA1 (", var_exp_grupo_bin[1], "%)"),
    y = paste0("PCoA2 (", var_exp_grupo_bin[2], "%)"),
    color = "Health Group"
  ) +
  theme_minimal(base_size = 16) +
  scale_color_manual(values = c("1" = "#1b9e77", 
                                "2" = "#d95f02", 
                                "3" = "#7570b3"))

print(g.bin)




# Salvar o gráfico

ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/PCoA_cluster", plot = g.bin width = 10, height = 6, dpi = 300)


# Ver qual grupo cada voluntário pertence
# Visualizar os primeiros voluntários e seus grupos
head(df_grupo_bin_clean %>% select(Sample.id, HealthGroup))

# Contagem de voluntários por grupo
table(df_grupo_bin_clean$HealthGroup)

library(dplyr)

#estatistica descritiva de cada grupo
df_grupo_bin_clean %>%
  group_by(HealthGroup) %>%
  summarise(across(-Sample.id, ~ mean(.x, na.rm = TRUE)))

# voluntarios do grupo 1: 

# Ver apenas os voluntários no grupo 1
grupo_1 <- df_grupo_bin_clean %>%
  filter(HealthGroup == 1)

# Visualizar
View(grupo_1)

grupo_1_ids <- grupo_1$Sample.id
print(grupo_1_ids)


#Código para gerar silhuetas para vários k
library(cluster)

# Preparar os dados sem Sample.id e sem HealthGroup (já removido antes)
dados_kmeans <- df_grupo_bin_clean %>% 
  select(-Sample.id, -HealthGroup)

# Calcular matriz de distância euclidiana
dist_matrix <- dist(dados_kmeans, method = "euclidean")

library(factoextra)
library(patchwork)

# Gerar plots de silhueta usando factoextra::fviz_silhouette()
plots_silhueta <- list()

for (k in 2:5) {
  set.seed(123)
  km <- kmeans(dados_kmeans, centers = k, nstart = 25)
  sil <- silhouette(km$cluster, dist_matrix)
  
  p <- fviz_silhouette(sil) + 
    ggtitle(paste("Silhouette - k =", k)) +
    theme_minimal(base_size = 14)
  
  plots_silhueta[[k]] <- p
}

# Juntar todos os gráficos com patchwork
(plots_silhueta[[2]] | plots_silhueta[[3]]) / (plots_silhueta[[4]] | plots_silhueta[[5]])


#gerar grafico do average silhouette score por n of cluster


library(cluster)
library(ggplot2)

# Calcular matriz de distância (euclidiana)
dist_matrix <- dist(dados_kmeans, method = "euclidean")

# Criar vetor para armazenar a silhueta média
sil_means <- data.frame(k = 2:10, sil_mean = NA)

# Calcular silhueta média para k = 2 a 10
for (i in 2:10) {
  set.seed(123)
  km <- kmeans(dados_kmeans, centers = i, nstart = 25)
  sil <- silhouette(km$cluster, dist_matrix)
  sil_means$sil_mean[i - 1] <- mean(sil[, "sil_width"])
}

# Plotar o gráfico
ggplot(sil_means, aes(x = k, y = sil_mean)) +
  geom_line(size = 1.2, color = "#1f78b4") +
  geom_point(size = 3, color = "#e31a1c") +
  labs(
    title = "Average Silhouette Score by Number of Clusters (k)",
    x = "Number of Clusters (k)",
    y = "Average Silhouette Score"
  ) +
  theme_minimal(base_size = 16)





#colorir por parametros

library(ggplot2)
library(dplyr)

# 1. Selecionar variáveis contínuas
vars_continuas <- c("Systolic", "Diastolic", "Weight", "BMI", 
                    "HbA1c", "Cholesterol", "LDL", "HDL", "Triglycerides", 
                    "GGT", "Glucose")

# 2. Juntar os metadados com as coordenadas da PCoA
pcoa_dados <- pcoa_df_grupo_bin %>%
  left_join(metadados_grupos_saude_binario %>% 
              select(Sample.id, all_of(vars_continuas)), 
            by = "Sample.id")

library(ggplot2)
library(dplyr)

# Criar pasta de saída (ajuste o caminho conforme quiser)
dir_out <- "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/PCoA_plots/PCoA_plots/PCoA_por_variavel"
dir.create(dir_out, showWarnings = FALSE)

# Loop para gerar e salvar os gráficos
for (var in vars_continuas) {
  p <- ggplot(pcoa_dados, aes(x = PCoA1, y = PCoA2, color = .data[[var]])) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_viridis_c(option = "D", direction = -1) +
    labs(
      title = paste("PCoA colored by", var),
      x = paste0("PCoA1 (", var_exp_grupo_bin[1], "%)"),
      y = paste0("PCoA2 (", var_exp_grupo_bin[2], "%)"),
      color = var
    ) +
    theme_minimal(base_size = 14)
  
  # Salvar como PNG
  ggsave(
    filename = paste0(dir_out, "/", "PCoA_", var, ".png"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}


#====== PCOA só com sim e nao (1 e 0)
#Código para PCoA com distância Euclidiana (apenas variáveis binárias puras)




#=========================================================#
#         Análise de Correspondencia multipla (MCA)
#=========================================================#

# Instalar caso ainda não tenha
#install.packages("FactoMineR")
#install.packages("factoextra")

# Carregar pacotes
library(FactoMineR)
library(factoextra)


# Preparar os dados (removendo Sample.id e HealthGroup, se ainda não fez)
dados_mca <- df_grupo_bin_clean %>%
  select(-Sample.id, -HealthGroup)

# Garantir que tudo seja fator
dados_mca <- dados_mca %>% mutate(across(everything(), as.factor))

# Rodar MCA
res.mca <- MCA(dados_mca, graph = FALSE)

#Visualizar com seus clusters (HealthGroup):
fviz_mca_ind(res.mca,
             habillage = df_grupo_bin_clean$HealthGroup, # cor por grupo
             addEllipses = TRUE, ellipse.level = 0.95,
             repel = TRUE,
             label = "none") +
  labs(title = "Multiple Correspondence Analysis (MCA) - Health Groups") +
  theme_minimal(base_size = 14)


# Ver quantas pessoas foram incluídas na MCA
nrow(res.mca$ind$coord)

# Ver quantas linhas foram retiradas
nrow(dados_mca) - nrow(res.mca$ind$coord)


#===> Incluir saudaveis

# 1. Remover apenas NAs, mas manter quem tem todos 0s
df_binarios_completos <- df_grupo_bin %>%
  drop_na()

# 2. Criar grupo 0 para quem tem todos os valores 0 nas variáveis (exceto Sample.id)
df_binarios_completos$HealthGroup <- ifelse(
  rowSums(df_binarios_completos[ , -1]) == 0, 0, NA  # -1 remove Sample.id
)

# 3. Preencher os grupos para os que já tinham cluster definido
df_binarios_completos$HealthGroup[is.na(df_binarios_completos$HealthGroup)] <- kmeans_result_grupo_bin$cluster

# 4. Agora, rodar a MCA novamente com todos os grupos
dados_mca <- df_binarios_completos %>%
  select(-Sample.id, -HealthGroup) %>%
  mutate(across(everything(), as.factor))

res.mca <- MCA(dados_mca, graph = FALSE)

# 5. Plotando com o novo grupo
fviz_mca_ind(res.mca,
             habillage = as.factor(df_binarios_completos$HealthGroup),
             addEllipses = TRUE,
             ellipse.level = 0.95,
             label = "none", repel = TRUE) +
  labs(title = "MCA com Grupo 0 (Completamente Saudáveis)") +
  theme_minimal(base_size = 14)

#======> Silhouete

# Carregar pacotes
library(cluster)
library(factoextra)

# 1. Coordenadas dos indivíduos na MCA
coords_mca <- as.data.frame(res.mca$ind$coord)

# 2. Agrupamentos (fator)
grupos_mca <- as.factor(df_binarios_completos$HealthGroup)

# 3. Matriz de distância entre indivíduos
dist_mca <- dist(coords_mca)

# 4. Calcular silhueta
silhueta_mca <- silhouette(as.numeric(grupos_mca), dist_mca)

# 5. Plotar em uma única figura
fviz_silhouette(silhueta_mca,
                palette = "Dark2",
                ggtheme = theme_minimal(base_size = 14),
                print.summary = TRUE) +
  labs(title = "Silhueta dos Grupos no Espaço MCA")


#calcular e plotar o average silhouette score no espaço MCA:

# Carregar pacotes
library(cluster)
library(ggplot2)
library(dplyr)

# Coordenadas dos indivíduos no espaço MCA
coords_mca <- as.data.frame(res.mca$ind$coord)

# Avaliar silhueta média para diferentes valores de k
silhouette_scores <- data.frame()

for (k in 2:10) {
  set.seed(123)
  km <- kmeans(coords_mca, centers = k, nstart = 25)
  sil <- silhouette(km$cluster, dist(coords_mca))
  avg_sil <- mean(sil[, 3])  # terceira coluna = largura da silhueta
  
  silhouette_scores <- rbind(silhouette_scores,
                             data.frame(k = k, avg_silhouette = avg_sil))
}

# Plotar
ggplot(silhouette_scores, aes(x = k, y = avg_silhouette)) +
  geom_line(group = 1, color = "steelblue", linewidth = 1) +
  geom_point(size = 3, color = "red") +
  labs(
    title = "Average Silhouette Score por Número de Clusters (k)",
    x = "Número de Clusters (k)",
    y = "Average Silhouette Score"
  ) +
  theme_minimal(base_size = 14)


#==============================================================#
#              Analise de Correspondencia
#==============================================================#

# Carregar pacotes
library(dplyr)
library(FactoMineR)
library(factoextra)

# ========================
# 1. Criar planilha binária (0 = saudável, 1 = alterado)
# ========================
# 1. Criar planilha binária (0 = saudável, 1 = alterado)
metadados_binarios_puros <- metadados.all.filtrado %>%
  mutate(
    High_Blood_Pressure = ifelse(Systolic >= 130 | Diastolic >= 85, 1, 0),
    High_Fasting_Glucose = ifelse(Glucose > 126, 1, 0),
    HbA1c_bin = ifelse(HbA1c >= 6.5, 1, 0),
    Obese = ifelse(BMI >= 30, 1, 0),
    High_Trig = ifelse(Triglycerides >= 200, 1, 0),
    High_Col = ifelse(Cholesterol > 200, 1, 0),
    Low_HDL = ifelse(HDL < 40, 1, 0),
    High_LDL = ifelse(LDL >= 160, 1, 0),
    High_GGT = case_when(
      Sex == "Masculino" & GGT > 48 ~ 1,
      Sex == "Feminino" & GGT > 32 ~ 1,
      TRUE ~ 0
    )
  ) %>%
  select(Sample.id, High_Blood_Pressure, High_Fasting_Glucose, HbA1c_bin,
         Obese, High_Trig, High_Col, Low_HDL, High_LDL, High_GGT)



  
  # ========================
# 2. Garantir que todas as colunas (exceto Sample.id) sejam fatores com níveis 0 e 1
# ========================

# 2. Garantir que todas as colunas (exceto Sample.id) sejam fatores com níveis 0 e 1
dados_ca_bin_puros <- metadados_binarios_puros %>%
  select(-Sample.id) %>%
  mutate(across(everything(), ~ factor(.x, levels = c(0, 1))))

# ========================
# 3. Rodar a MCA
# ========================
res.mca.binarios <- MCA(dados_ca_bin_puros, graph = FALSE)

# ========================
# 4. Identificar voluntários completamente saudáveis (tudo 0)
# ========================
saudaveis_ids <- metadados_binarios_puros %>%
  filter(rowSums(select(., -Sample.id)) == 0) %>%
  pull(Sample.id)

# ========================
# 5. Construir vetor de grupos (0 para saudáveis, 1-3 do kmeans)
# ========================

# IDs usados na MCA
ids_mca <- metadados_binarios_puros$Sample.id

# Obter os grupos do kmeans para os IDs com dados
healthgroup_bin <- df_grupo_bin_clean %>%
  filter(Sample.id %in% ids_mca) %>%
  select(Sample.id, HealthGroup)

# Criar vetor final com grupo 0 para saudáveis
healthgroup_completo <- data.frame(Sample.id = ids_mca) %>%
  left_join(healthgroup_bin, by = "Sample.id") %>%
  mutate(HealthGroup = ifelse(Sample.id %in% saudaveis_ids, "0", as.character(HealthGroup))) %>%
  mutate(HealthGroup = factor(HealthGroup, levels = c("0", "1", "2", "3")))

# ========================
# 6. Plotar a MCA colorida por grupo (com elipses)
# ========================
fviz_mca_ind(res.mca.binarios,
             habillage = healthgroup_completo$HealthGroup,
             palette = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"),
             addEllipses = TRUE,
             ellipse.level = 0.95,
             label = "none",
             repel = TRUE) +
  labs(title = "MCA - Grupos de Saúde - 
              Zero = completamente saudáveis") +
  theme_minimal(base_size = 14)


ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/PCoA_MCA_grupos_saude_contraste.png", width = 10, height = 6, dpi = 300, bg = "white")


#====> colorir com parametros continuos


# Carregar pacotes
library(dplyr)
library(ggplot2)
library(viridis)

# 1. Coordenadas fixas da MCA
coords <- as.data.frame(res.mca.binarios$ind$coord)
colnames(coords)[1:2] <- c("Dim.1", "Dim.2")  # <- ESSA LINHA RESOLVE
coords$Sample.id <- metadados_binarios_puros$Sample.id

# 2. Variáveis contínuas para colorir
variaveis_continuas <- c("Glucose", "HbA1c", "Cholesterol", "LDL", "HDL",
                         "Triglycerides", "TGO", "TGP", "GGT", "BMI")

# 3. Loop para gerar gráficos
for (var in variaveis_continuas) {
  
  dados_plot <- coords %>%
    left_join(metadados.all.filtrado %>% select(Sample.id, all_of(var)), by = "Sample.id")
  
  p <- ggplot(dados_plot, aes(x = Dim.1, y = Dim.2, color = .data[[var]])) +
    geom_jitter(width = 0.03, height = 0.03, size = 3, alpha = 0.9) +
    scale_color_viridis_c(option = "D", direction = -1, na.value = "gray90") +
    labs(
      title = paste("MCA - Colorido por", var),
      x = paste0("Dim1 (", round(res.mca.binarios$eig[1,2], 1), "%)"),
      y = paste0("Dim2 (", round(res.mca.binarios$eig[2,2], 1), "%)"),
      color = var
    ) +
    theme_minimal(base_size = 14)
  
  print(p)
  
  ggsave(filename = paste0("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/PCoA_plots/PCoA_plots/PCoA_MCA_colorido_por_", var, ".png"),
         plot = p,
         width = 9, height = 6, dpi = 300, bg = "white")
}





#=====================================================#

#       Score metabolico 

#=====================================================#


#Criar metadados binários para fazer kmeans (distância euclideana) e pcoa

metadados_grupos_saude_binario <- metadados.all.filtrado %>%
  mutate(
    # Pressão arterial
    High_Blood_Pressure = ifelse(!is.na(Systolic) & !is.na(Diastolic) &
                                   (Systolic >= 130 | Diastolic >= 85), 1, 0),
    
    # Glucose de jejum
    High_Fasting_Glucose = case_when(
      is.na(Glucose) ~ NA_integer_,
      Glucose <= 100 ~ 0,                 # normal
      Glucose > 100 & Glucose < 126 ~ 1,  # alterada (pré-diabetes)
      Glucose >= 126 ~ 2                  # muito alta (diabetes)
    ),
    
    # HbA1c em porcentagem
    HbA1c_cat = case_when(
      is.na(HbA1c) ~ NA_integer_,
      HbA1c < 6.0 ~ 0,
      HbA1c >= 6.0 & HbA1c < 6.5 ~ 1,
      HbA1c >= 6.5 ~ 2
    ),
    
    
    # BMI categorizado por faixa etária
    BMI_cat = case_when(
      is.na(BMI) | is.na(Age) ~ NA_integer_,
      Age >= 60 & BMI < 22 ~ 1,
      Age >= 60 & BMI >= 22 & BMI < 27 ~ 0,
      Age >= 60 & BMI >= 27 & BMI < 30 ~ 1,
      Age >= 60 & BMI >= 30 & BMI < 35 ~ 2,
      Age >= 60 & BMI >= 35 ~ 3,
      Age < 60 & BMI < 18.5 ~ 1,
      Age < 60 & BMI >= 18.5 & BMI < 25 ~ 0,
      Age < 60 & BMI >= 25 & BMI < 30 ~ 1,
      Age < 60 & BMI >= 30 & BMI < 35 ~ 2,
      Age < 60 & BMI >= 35 & BMI < 40 ~ 3,
      Age < 60 & BMI >= 40 ~ 4
    ),
    
    # Triglicerídeos
    TRIG_cat = case_when(
      is.na(Triglycerides) ~ NA_integer_,
      Triglycerides < 150 ~ 0,
      Triglycerides >= 150 & Triglycerides < 200 ~ 1,
      Triglycerides >= 200 & Triglycerides < 500 ~ 2,
      Triglycerides >= 500 ~ 3
    ),
    
    # Cholesterol total
    COL_cat = case_when(
      is.na(Cholesterol) ~ NA_integer_,
      Cholesterol <= 190 ~ 0,
      Cholesterol > 190 & Cholesterol <= 240 ~ 1,
      Cholesterol > 240 ~ 2
    ),
    
    
    # HDL
    HDL_cat = case_when(
      is.na(HDL) ~ NA_integer_,
      HDL < 40 ~ 1,        # risco
      HDL >= 40 ~ 0 ),       # normal ou bom
    
    
    
    # LDL
    LDL_cat = case_when(
      is.na(LDL) ~ NA_integer_,
      LDL < 100 ~ 0,
      LDL >= 100 & LDL < 130 ~ 0,
      LDL >= 130 & LDL < 160 ~ 1,
      LDL >= 160 & LDL < 190 ~ 2,
      LDL >= 190 ~ 3
    ),
    
    # GGT por sexo
    GGT_cat = case_when(
      is.na(GGT) | is.na(Sex) ~ NA_integer_,
      Sex == "Masculino" & GGT <= 48 ~ 0,
      Sex == "Masculino" & GGT > 48 ~ 1,
      Sex == "Feminino" & GGT <= 32 ~ 0,
      Sex == "Feminino" & GGT > 32 ~ 1
    )
  )


#Metabolic Score

metadados_grupos_saude_binario <- metadados_grupos_saude_binario %>%
  mutate(
    Metabolic_Score = rowSums(across(c(
      High_Blood_Pressure, High_Fasting_Glucose, HbA1c_cat,
      BMI_cat, TRIG_cat, COL_cat, HDL_cat, LDL_cat, GGT_cat
    ), ~replace_na(., 0)))
  )


#nomear colunas
metadados_grupos_saude_binario <- metadados_grupos_saude_binario %>%
  mutate(
    Health_Status = case_when(
      Metabolic_Score <= 1 ~ "Saudável",
      Metabolic_Score >= 2 & Metabolic_Score <= 4 ~ "Transição",
      Metabolic_Score >= 5 ~ "Doente"
    )
  )


#plotar

library(ggplot2)

ggplot(metadados_grupos_saude_binario, aes(x = Health_Status, y = Metabolic_Score, fill = Health_Status)) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  scale_fill_manual(values = c("Saudável" = "#1b9e77", "Transição" = "#d95f02", "Doente" = "#7570b3")) +
  labs(title = "Distribuição do Metabolic Score por Grupo de Saúde",
       x = "Grupo de Saúde", y = "Metabolic Score") +
  theme_minimal(base_size = 14)


#comparar parametros entre os grupos:

ggplot(metadados_grupos_saude_binario, aes(x = Health_Status, y = Glucose, fill = Health_Status)) +
  geom_boxplot() +
  labs(title = "Glicemia em jejum por Grupo de Saúde", x = "", y = "Glucose (mg/dL)") +
  theme_minimal(base_size = 14)

#Ver se é significante
kruskal.test(Glucose ~ Health_Status, data = metadados_grupos_saude_binario)

#install.packages("FSA")  # se ainda não tiver
library(FSA)

# Teste de Dunn com correção de Bonferroni
dunnTest(Glucose ~ Health_Status, data = metadados_grupos_saude_binario, method = "bonferroni")



library(ggplot2)
library(ggpubr)


#Grafico com valor de kruskall-wallis

ggplot(metadados_grupos_saude_binario, aes(x = Health_Status, y = Glucose, fill = Health_Status)) +
  geom_boxplot() +
  labs(title = "Glicemia em jejum por Grupo de Saúde", x = "", y = "Glucose (mg/dL)") +
  theme_minimal(base_size = 14) +
  stat_compare_means(method = "kruskal.test", label.y = max(metadados_grupos_saude_binario$Glucose, na.rm = TRUE) + 5)


# Instale se necessário:
# install.packages("ggpubr")

library(ggpubr)

# Boxplot com estatísticas do Kruskal-Wallis e post-hoc de Dunn
ggplot(metadados_grupos_saude_binario, aes(x = Health_Status, y = Glucose, fill = Health_Status)) +
  geom_boxplot() +
  stat_compare_means(method = "kruskal.test", label.y = max(metadados_grupos_saude_binario$Glucose, na.rm = TRUE) + 10) +  # Teste global
  stat_compare_means(comparisons = list(
    c("Doente", "Saudável"),
    c("Doente", "Transição"),
    c("Saudável", "Transição")
  ),
  method = "wilcox.test", label = "p.signif", hide.ns = TRUE,
  label.y = c(
    max(metadados_grupos_saude_binario$Glucose, na.rm = TRUE) + 30,
    max(metadados_grupos_saude_binario$Glucose, na.rm = TRUE) + 20,
    max(metadados_grupos_saude_binario$Glucose, na.rm = TRUE) + 10)
  ) +
  labs(title = "Glicemia em jejum por Grupo de Saúde", x = "", y = "Glucose (mg/dL)") +
  theme_minimal(base_size = 14)


#HbA1c

#Ver se é significante
kruskal.test(HbA1c ~ Health_Status, data = metadados_grupos_saude_binario)


# Teste de Dunn com correção de Bonferroni
dunnTest(HbA1c ~ Health_Status, data = metadados_grupos_saude_binario, method = "bonferroni")




#Grafico com valor de kruskall-wallis

ggplot(metadados_grupos_saude_binario, aes(x = Health_Status, y = HbA1c, fill = Health_Status)) +
  geom_boxplot() +
  labs(title = "HbA1c por Grupo de Saúde", x = "", y = "HbA1c") +
  theme_minimal(base_size = 14) +
  stat_compare_means(method = "kruskal.test", label.y = max(metadados_grupos_saude_binario$Glucose, na.rm = TRUE) + 5)


# Instale se necessário:
# install.packages("ggpubr")

library(ggpubr)

# Boxplot com estatísticas do Kruskal-Wallis e post-hoc de Dunn
ggplot(metadados_grupos_saude_binario, aes(x = Health_Status, y = Glucose, fill = Health_Status)) +
  geom_boxplot() +
  stat_compare_means(method = "kruskal.test", label.y = max(metadados_grupos_saude_binario$Glucose, na.rm = TRUE) + 10) +  # Teste global
  stat_compare_means(comparisons = list(
    c("Doente", "Saudável"),
    c("Doente", "Transição"),
    c("Saudável", "Transição")
  ),
  method = "wilcox.test", label = "p.signif", hide.ns = TRUE,
  label.y = c(
    max(metadados_grupos_saude_binario$Glucose, na.rm = TRUE) + 30,
    max(metadados_grupos_saude_binario$Glucose, na.rm = TRUE) + 20,
    max(metadados_grupos_saude_binario$Glucose, na.rm = TRUE) + 10)
  ) +
  labs(title = "Glicemia em jejum por Grupo de Saúde", x = "", y = "Glucose (mg/dL)") +
  theme_minimal(base_size = 14)


#====>  Grafico de pizza representando os grupos

library(ggplot2)
library(dplyr)



#metadados_grupos_saude_binario <- metadados_grupos_saude_binario %>%
#  mutate(Health_Status = recode(Health_Status,
#                                "Saudável" = "Low Metabolic Risk",
#                                "Transição" = "Intermediate Metabolic Risk",
#                                "Doente" = "High Metabolic Risk"))


table(metadados_grupos_saude_binario$Health_Status)



library(ggplot2)
library(dplyr)

# Criar a tabela de frequência manualmente
grupo_tab <- as.data.frame(table(metadados_grupos_saude_binario$Health_Status))
colnames(grupo_tab) <- c("Health_Status", "n")

# Calcular porcentagens
grupo_tab$percent <- grupo_tab$n / sum(grupo_tab$n) * 100
grupo_tab$label <- paste0("\n", grupo_tab$n, " (", round(grupo_tab$percent, 1), "%)")


#ESSE GRAFICO!!!

# Gráfico com letras maiores
g_pizza <- ggplot(grupo_tab, aes(x = "", y = n, fill = Health_Status)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void(base_size = 18) +  # Aumenta tamanho da fonte
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 6) +  # Texto maior
  labs(title = "Metabolic Risk Distribution (n=128)") +
  scale_fill_brewer(palette = "Set2")

# Salvar o gráfico em alta resolução
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/pizza_grupo_metabolico.png",
       plot = g_pizza, width = 10, height = 8, dpi = 300)


#===== UPSET =====#
library(dplyr)
library(tidyr)

metadados_binarios_puros <- metadados.all.filtrado %>%
  mutate(
    # Pressão arterial
    High_Blood_Pressure = ifelse(Systolic >= 130 | Diastolic >= 85, 1, 0),
    
    # Glucose
    High_Fasting_Glucose = case_when(
      Glucose <= 100 ~ 0,
      Glucose > 100 & Glucose < 126 ~ 1,
      Glucose >= 126 ~ 1,
      TRUE ~ NA_integer_
    ),
    
    # HbA1c
    HbA1c_cat = case_when(
      HbA1c < 6.0 ~ 0,
      HbA1c >= 6.0 & HbA1c < 6.5 ~ 1,
      HbA1c >= 6.5 ~ 1,
      TRUE ~ NA_integer_
    ),
    
    # BMI categorizado por idade
    BMI_cat = case_when(
      Age >= 60 & BMI < 22 ~ 0,
      Age >= 60 & BMI >= 22 & BMI < 27 ~ 0,
      Age >= 60 & BMI >= 27 & BMI < 30 ~ 0,
      Age >= 60 & BMI >= 30 & BMI < 35 ~ 1,
      Age >= 60 & BMI >= 35 ~ 1,
      Age < 60 & BMI < 18.5 ~ 1,
      Age < 60 & BMI >= 18.5 & BMI < 25 ~ 0,
      Age < 60 & BMI >= 25 & BMI < 30 ~ 0,
      Age < 60 & BMI >= 30 & BMI < 35 ~ 1,
      Age < 60 & BMI >= 35 & BMI < 40 ~ 1,
      Age < 60 & BMI >= 40 ~ 1,
      TRUE ~ NA_integer_
    ),
    
    # Triglicerídeos
    TRIG_cat = case_when(
      Triglycerides < 150 ~ 0,
      Triglycerides >= 150 & Triglycerides < 200 ~ 1,
      Triglycerides >= 200 & Triglycerides < 500 ~ 1,
      Triglycerides >= 500 ~ 1,
      TRUE ~ NA_integer_
    ),
    
    # Cholesterol total
    COL_cat = case_when(
      Cholesterol <= 200 ~ 0,
      Cholesterol > 200 & Cholesterol <= 240 ~ 1,
      Cholesterol > 240 ~ 1,
      TRUE ~ NA_integer_
    ),
    
    # HDL
    HDL_cat = case_when(
      HDL < 40 ~ 1,
      HDL >= 40 ~ 0,
      TRUE ~ NA_integer_
    ),
    
    # LDL
    LDL_cat = case_when(
      LDL < 130 ~ 0,
      LDL >= 130 & LDL < 160 ~ 1,
      LDL >= 160 & LDL < 190 ~ 1,
      LDL >= 190 ~ 1,
      TRUE ~ NA_integer_
    ),
    
    # GGT por sexo
    GGT_cat = case_when(
      Sex == "Masculino" & GGT > 48 ~ 1,
      Sex == "Masculino" & GGT <= 48 ~ 0,
      Sex == "Feminino" & GGT > 32 ~ 1,
      Sex == "Feminino" & GGT <= 32 ~ 0,
      TRUE ~ NA_integer_
    )
  ) %>%
  # Substituir todos os NA das variáveis binárias/categóricas por 0
  mutate(across(
    c(High_Blood_Pressure, High_Fasting_Glucose, HbA1c_cat, BMI_cat, TRIG_cat,
      COL_cat, HDL_cat, LDL_cat, GGT_cat),
    ~ replace_na(., 0)
  ))




#install.packages("ComplexUpset")
library(ComplexUpset)
library(ggplot2)


# Selecionar colunas binárias (0 ou 1)
colunas_binarias <- c("High_Blood_Pressure",              
                      "High_Fasting_Glucose",              "HbA1c_cat" ,                         "BMI_cat" ,                          
                      "TRIG_cat",                           "COL_cat"  ,                          "HDL_cat" ,                          
                      "LDL_cat" ,                           "GGT_cat"                           )

# Subset do dataframe com essas variáveis
dados_upset <- metadados_binarios_puros[, colunas_binarias]

# Garantir que todas sejam numéricas (necessário para ComplexUpset)
dados_upset[] <- lapply(dados_upset, as.numeric)

library(ComplexUpset)
library(ggplot2)


# Gráfico funcional sem mapeamentos problemáticos


# Gráfico final = idosos
p <-upset(
  dados_upset,
  intersect = colunas_binarias,
  base_annotations = list(
    'Intersection size' = intersection_size(
      mapping = aes(fill = "azul"),
      text = list(size = 4)
    ) + 
      scale_fill_manual(values = c("azul" = "#003399"), guide = "none")
  ),
  set_sizes = upset_set_size(
    mapping = aes(fill = "vermelho")
  ) + 
    scale_fill_manual(values = c("vermelho" = "#8B0000"), guide = "none"),
  width_ratio = 0.2
) +
  labs(title = "Metabolic Health Profile (n=128)")


# Exibir
print(p)


ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/upset_metabolic_healthprofile_completen130.png", plot = p, width = 18, height = 10, dpi = 300)

#======> Citocinas e Metabolic Health Score

#Boxplots com metadados.metabolic.score e tratamento de NA

library(ggplot2)
library(ggpubr)

# Lista das citocinas
citocinas <- c("IL17A", "IFNGamma", "TNF", "IL10", "IL6", "IL4", "IL2")

# Loop para gerar boxplots + Kruskal
for (cit in citocinas) {
  dados <-metadados.all.filtrado[, c("Health_Status", cit)]
  dados <- dados[!is.na(dados[[cit]]), ]  # remove NAs da citocina
  
  # Teste de Kruskal
  p.kruskal <- kruskal.test(reformulate("Health_Status", cit), data = dados)
  print(paste(cit, "- p =", round(p.kruskal$p.value, 4)))
  
  # Gráfico
  p <- ggplot(dados, aes(x = Health_Status, y = .data[[cit]], fill = Health_Status)) +
    geom_boxplot() +
    labs(title = paste0(cit, " por grupo metabólico (p = ", round(p.kruskal$p.value, 3), ")"),
         x = "Health Status", y = paste(cit, "(pg/mL)")) +
    theme_minimal(base_size = 14) +
    scale_fill_brewer(palette = "Spectral") +
    theme(legend.position = "none")
  
  print(p)
}


#====== Criar um índice inflamatório e outro antiinflamatório (z-score)====#

# Selecionar subconjuntos
inflamatorias <- c("IL6", "TNF", "IFNGamma", "IL17A")
antiinflamatorias <- c("IL10", "IL4", "IL2")

# Escalar (z-score) as citocinas
z_inflam <- scale(metadados_grupos_saude_binario[, inflamatorias])
z_antiinflam <- scale(metadados_grupos_saude_binario[, antiinflamatorias])

# Calcular o escore médio por indivíduo
metadados_grupos_saude_binario$Inflam_Score <- rowMeans(z_inflam, na.rm = TRUE)
metadados_grupos_saude_binario$AntiInflam_Score <- rowMeans(z_antiinflam, na.rm = TRUE)

#=====> Comparar os escores entre os grupos


# Kruskal-Wallis para cada escore
kruskal.test(Inflam_Score ~ Health_Status, data = metadados_grupos_saude_binario)
kruskal.test(AntiInflam_Score ~ Health_Status, data = metadados_grupos_saude_binario)

#=======> Boxplots com interpretação

library(ggplot2)

ggplot(metadados_grupos_saude_binario, aes(x = Health_Status, y = Inflam_Score, fill = Health_Status)) +
  geom_boxplot() +
  labs(title = "Escore inflamatório (z-score médio)", y = "Inflam_Score") +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Spectral")

ggplot(metadados_grupos_saude_binario, aes(x = Health_Status, y = AntiInflam_Score, fill = Health_Status)) +
  geom_boxplot() +
  labs(title = "Escore anti-inflamatório (z-score médio)", y = "AntiInflam_Score") +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Spectral")


#=======>


library(pheatmap)

# Corrigir lista de variáveis (sem duplicata de BMI)
metabolicos <- unique(c("HbA1c", "Glucose", "INSULINA", "Cholesterol", 
                        "LDL", "HDL", "Triglycerides", "TGO", "TGP", "GGT", 
                        "BMI", "PCR"))

# Citocinas + escores
citocinas <- c("IL17A", "IFNGamma", "TNF", "IL10", "IL6", "IL4", "IL2", "Inflam_Score", "AntiInflam_Score")

# Inicializar matrizes
cor_matrix <- matrix(NA, nrow = length(citocinas), ncol = length(metabolicos),
                     dimnames = list(citocinas, metabolicos))
p_matrix <- cor_matrix

# Calcular correlação de Spearman com verificação
for (i in citocinas) {
  for (j in metabolicos) {
    x <- metadados_grupos_saude_binario[[i]]
    y <- metadados_grupos_saude_binario[[j]]
    
    if (is.numeric(x) && is.numeric(y)) {
      dados_validos <- complete.cases(x, y)
      if (sum(dados_validos) > 5) {
        teste <- suppressWarnings(cor.test(x[dados_validos], y[dados_validos], method = "spearman"))
        cor_matrix[i, j] <- round(teste$estimate, 2)
        p_matrix[i, j] <- teste$p.value
      }
    }
  }
}

# Remover colunas/linhas com NA em excesso
cor_matrix_clean <- cor_matrix[rowSums(is.na(cor_matrix)) < ncol(cor_matrix), 
                               colSums(is.na(cor_matrix)) < nrow(cor_matrix)]

p_matrix_clean <- p_matrix[rownames(cor_matrix_clean), colnames(cor_matrix_clean)]

# Criar matriz de asteriscos e rótulos
sig_labels <- ifelse(p_matrix_clean < 0.05, "*", "")
labels <- matrix(paste0(cor_matrix_clean, sig_labels), nrow = nrow(cor_matrix_clean),
                 dimnames = dimnames(cor_matrix_clean))

#====> Preparar os dados em formato longo

library(dplyr)
library(tidyr)
library(ggplot2)

# Converter a matriz de correlação para long format
cor_df <- as.data.frame(cor_matrix_clean) %>%
  mutate(Citocina = rownames(.)) %>%
  pivot_longer(-Citocina, names_to = "Parametro", values_to = "Correlacao")

# Adicionar os p-valores e os asteriscos de significância
p_df <- as.data.frame(p_matrix_clean) %>%
  mutate(Citocina = rownames(.)) %>%
  pivot_longer(-Citocina, names_to = "Parametro", values_to = "p_value") %>%
  mutate(asterisco = ifelse(p_value < 0.05, "*", ""))

# Juntar as duas informações
heatmap_df <- left_join(cor_df, p_df, by = c("Citocina", "Parametro")) %>%
  mutate(label = paste0(round(Correlacao, 2), asterisco))

#====> Criar heatmap com ggplot2

heatmap_plot <- ggplot(heatmap_df, aes(x = Parametro, y = Citocina, fill = Correlacao)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 4) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick3",
                       midpoint = 0, limit = c(-1, 1), name = "Spearman\nrho") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) +
  labs(title = "Correlação Spearman: Citocinas vs. Parâmetros Metabólicos",
       x = "Parâmetros Metabólicos", y = "Citocinas")

# Visualizar
print(heatmap_plot)

# Salvar
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/heatmap_citocinas_ggplot.png",
       plot = heatmap_plot, width = 12, height = 8, dpi = 300)



#=====> testar citocinas com alpha
library(ggplot2)

# 1. Listas de variáveis
citocinas <- c("IL17A", "IFNGamma", "TNF", "IL10", "IL6", "IL4", "IL2", "Inflam_Score", "AntiInflam_Score")
alphas <- c("shannon_entropy", "simpson", "pielou_evenness", "chao1", "faith_pd", "observed_features")



# 2. Loop para cada par (citocina x índice de diversidade)
for (cytokine in citocinas) {
  for (alpha in alphas) {
    
    # Remover NA dos dois vetores
    dados_validos <- complete.cases(metadados.metabolic.score[[cytokine]], metadados.metabolic.score[[alpha]])
    
    if (sum(dados_validos) > 5) {
      # Calcular correlação Spearman
      resultado <- cor.test(metadados.metabolic.score[[cytokine]][dados_validos],
                            metadados.metabolic.score[[alpha]][dados_validos],
                            method = "spearman")
      
      # Criar o gráfico
      p <- ggplot(metadados.metabolic.score[dados_validos, ], 
                  aes_string(x = alpha, y = cytokine)) +
        geom_point(color = "#1f78b4", size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", se = FALSE, color = "darkred") +
        theme_minimal(base_size = 14) +
        labs(
          title = paste("R =", round(resultado$estimate, 2), "- p =", signif(resultado$p.value, 3)),
          subtitle = paste(cytokine, "vs", alpha),
          x = alpha,
          y = cytokine
        )
      
      # Salvar em PNG
      nome_arquivo <- paste0("scatter_", cytokine, "_vs_", alpha, ".png")
      ggsave(nome_arquivo, p, width = 6, height = 5, dpi = 300)
    }
  }
}


