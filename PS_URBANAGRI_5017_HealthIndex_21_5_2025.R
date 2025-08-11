
#Codigo para GMW2

library(tibble)
library(purrr)

arquivos <- list.files(
  path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/Metagenomica/gmwi2",
  pattern = "_GMWI2\\.txt$",
  full.names = TRUE
)

ler_gmwi2_simples <- function(arquivo) {
  valor <- scan(arquivo, what = numeric(), quiet = TRUE)
  amostra <- gsub("_GMWI2\\.txt", "", basename(arquivo))
  tibble(SampleID = amostra, GMWI2 = valor)
}

resultados_gmwi2 <- map_dfr(arquivos, ler_gmwi2_simples)


#Juntar em metadados.all.filtrado e colocar a classe de cada pessoa.


metadados.all.filtrado$GMWI2_class <- case_when(
  is.na(metadados.all.filtrado$GMWI2) ~ NA_character_,
  metadados.all.filtrado$GMWI2 >= 1.0 ~ "Very Healthy",
  metadados.all.filtrado$GMWI2 >  0.5 ~ "Healthy",
  metadados.all.filtrado$GMWI2 >  0.0 ~ "Slightly Healthy",
  metadados.all.filtrado$GMWI2 == 0.0 ~ "Indeterminate",
  metadados.all.filtrado$GMWI2 > -0.5 ~ "Slightly Unhealthy",
  metadados.all.filtrado$GMWI2 > -1.0 ~ "Unhealthy",
  metadados.all.filtrado$GMWI2 <= -1.0 ~ "Very Unhealthy"
)


#install.packages("writexl")  # só precisa rodar uma vez
library(writexl)

write_xlsx(metadados.all.filtrado, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados_all_filtrado.xlsx")



#=================

#Criar em metadados.all.filtrado uma coluna chamada Pulse que seja Systolic - Diastolic
metadados.all.filtrado <- metadados.all.filtrado %>%
  mutate(pulse = Systolic - Diastolic)




# Selecionar variáveis de interesse
vars_metabolicas <- c("Sample.id", "GLICOSE", "TRIGLICERIDES", "COLESTEROL",
                      "LDL", "HDL", "HbA1c", "pulse",
                      "PCR", "TGO", "TGP", "GGT")

# Separar a coluna de ID
sample_ids <- metadados.all.filtrado$Sample.id

# Padronizar apenas as variáveis numéricas / Remover NAs e padronizar (z-score)
df_scaled <- metadados.all.filtrado %>%
  select(-Sample.id) %>%
  select(all_of(vars_metabolicas[-1])) %>%  # remove "Sample.id" do vetor
  scale(center = TRUE, scale = TRUE) %>%
  as.data.frame()

# Adicionar o Sample.id de volta
df_scaled$Sample.id <- sample_ids

#Código para reorganizar as colunas (Sample.id primeiro)
df_scaled <- df_scaled %>%
  relocate(Sample.id, .before = 1)


#Remover NAs
df_scaled_clean <- df_scaled %>% drop_na()

#Ou, se quiser ver quantas linhas foram removidas:
#n_inicial <- nrow(df_scaled)
#df_scaled_clean <- df_scaled %>% drop_na()
#n_final <- nrow(df_scaled_clean)
#cat("Linhas removidas:", n_inicial - n_final, "\n")


#================> Criar 3 grupos separados (super saudavel, intermediario e )

# 1. Rodar o k-means com 3 grupos

set.seed(123)  # Para reprodutibilidade

# Rodar k-means (removendo Sample.id antes)
kmeans_result <- kmeans(df_scaled_clean %>% select(-Sample.id), centers = 3)

# Ver quantos em cada grupo
table(kmeans_result$cluster)

#Adicionar a coluna de grupo ao dataframe
# Adicionar o cluster ao df_scaled
df_scaled_clean$HealthGroup <- as.factor(kmeans_result$cluster)

## Juntar com metadados.all.filtrado pelo Sample.id
metadados_clusterizado <- metadados.all.filtrado %>%
  inner_join(df_scaled_clean %>% select(Sample.id, HealthGroup), by = "Sample.id")


# Observar médias por grupo
df_scaled_clean %>%
  group_by(HealthGroup) %>%
  summarise(across(where(is.numeric), mean))

# Depois de avaliar, renomear:
metadados_clusterizado$HealthGroup_named <- recode_factor(
  metadados_clusterizado$HealthGroup,
  `1` = "Optimal Health",
  `2` = "High Risk Metabolic Profile",
  `3` = "In Transition"
)

#============> PCoA

#fazer uma PCoA com seus dados de saúde:

# Remover Sample.id e a coluna de grupo antes de calcular distância
dados_pcoa <- df_scaled_clean %>%
  select(-Sample.id, -HealthGroup)

# Calcular matriz de distância
dist_matrix <- dist(dados_pcoa, method = "euclidean")

# Rodar PCoA
pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)

#pcoa_result <- pcoa(dist_matrix)
#pcoa_result$values$Relative_eig  # proporção explicada por cada eixo

#round(pcoa_result$values$Relative_eig[1] * 100, 2)


# Criar dataframe com coordenadas
pcoa_df <- as.data.frame(pcoa_result$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")

# Adicionar os grupos
pcoa_df$Sample.id <- df_scaled_clean$Sample.id
pcoa_df$HealthGroup <- metadados_clusterizado$HealthGroup

library(ggplot2)

p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = HealthGroup)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCoA Based on Metabolic and Inflammatory Markers",
    x = "PCoA1",
    y = "PCoA2",
    color = "Health Group"
  ) +
  theme_minimal() +
  scale_color_manual(values = c(
    "1" = "#1b9e77",
    "2" = "#d95f02",
    "3" = "#7570b3"
  ))



print(p)


#=====> Permanova

library(vegan)
permanova_health_group <- adonis2(dist_matrix ~ metadados_clusterizado$HealthGroup)

r2 <- round(permanova_health_group$R2[1], 3)
pval <- format.pval(permanova_health_group$`Pr(>F)`[1], digits = 3, eps = .001)


#!!!! ====== Este grafico abaixo ======= !!!!#

p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = HealthGroup)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCoA Based on Metabolic and Inflammatory Markers (n=66)",
    x = "PCoA1",
    y = "PCoA2",
    color = "Health Group",
    caption = paste0("PERMANOVA: R² = ", r2, ", p = ", pval)
  ) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c(
    "1" = "#1b9e77",
    "2" = "#d95f02",
    "3" = "#7570b3"
  ))

ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/PCoA_HealthGroups.png", plot = p, width = 10, height = 6, dpi = 300)


#========> Colorir PCoA com todos os parametros de saúde, pra entender quais tem mais peso

library(dplyr)


#Joins da pcoa com os parametros
pcoa_merged <- pcoa_df %>%
  left_join(metadados.all.filtrado %>% 
              select(Sample.id, Age, HbA1c, INSULINA, GLICOSE, pulse, COLESTEROL,
                     LDL, HDL, VLDL, TRIGLICERIDES, TGO, TGP, GGT, IMC, PCR, IFNGamma, IL2, IL4, IL6, IL10, IL17A, TNF), 
            by = "Sample.id")
library(ggplot2)

#HbA1c

ggplot(pcoa_merged, aes(x = PCoA1, y = PCoA2, color = HbA1c)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_gradient(low = "#d9f0a3", high = "#005a32") +
  labs(title = "PCoA colored by HbA1c",
       x = "PCoA1", y = "PCoA2", color = "HbA1c") +
  theme_minimal(base_size = 16)


#IL6

ggplot(pcoa_merged, aes(x = PCoA1, y = PCoA2, color = IL6)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_gradient(low = "#fee0d2", high = "#de2d26") +
  labs(title = "PCoA colored by IL-6",
       x = "PCoA1", y = "PCoA2", color = "IL-6") +
  theme_minimal(base_size = 16)


#loop para todos os parametros

library(ggplot2)
library(scales)

# Lista de variáveis que você quer colorir (cada uma de cada vez)
variaveis <- c("Age", "HbA1c", "GLICOSE", "INSULINA", "pulse", "COLESTEROL", 
               "LDL", "HDL", "VLDL", "TRIGLICERIDES", "TGO", "TGP", "GGT", 
               "IMC", "IFNGamma", "IL2", "IL4", "IL6", 
               "IL10", "IL17A", "TNF", "PCR")



# Loop para gerar gráfico por variável
for (var in variaveis) {
  p <- ggplot(pcoa_merged, aes(x = PCoA1, y = PCoA2)) +
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
    filename = paste0("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/PCoA_plots/PCoA_plots/PCoA_", var, ".png"),
    plot = p,
    width = 10, height = 7, dpi = 300
  )
}



#====================Nova classificacao

metadados_grupos_saude_binario <- metadados.all.filtrado %>%
  mutate(
    # Pressão arterial sistólica e diastólica em 3 categorias
    High_Blood_Pressure = case_when(
      is.na(Systolic) | is.na(Diastolic) ~ NA_integer_,
      Systolic >= 140 | Diastolic >= 90 ~ 2,  # muito alta
      Systolic >= 130 | Diastolic >= 85 ~ 1,  # alta
      TRUE ~ 0  # normal
    ),
    
    # Glicose de jejum em 3 categorias
    High_Fasting_Glucose = case_when(
      is.na(GLICOSE) ~ NA_integer_,
      GLICOSE >= 127 ~ 2,  # muito alta
      GLICOSE >= 100 ~ 1,  # alta
      TRUE ~ 0  # normal
    ),
    
    # HbA1c em porcentagem
    HbA1c_cat = case_when(
      is.na(HbA1c) ~ NA_integer_,
      HbA1c < 6.0 ~ 0,
      HbA1c >= 6.0 & HbA1c < 6.5 ~ 1,
      HbA1c >= 6.5 ~ 2
    ),
    
    
    # IMC categorizado por faixa etária
    IMC_cat = case_when(
      is.na(IMC) | is.na(Age) ~ NA_integer_,
      Age >= 60 & IMC < 22 ~ 0,
      Age >= 60 & IMC >= 22 & IMC < 27 ~ 0,
      Age >= 60 & IMC >= 27 & IMC < 30 ~ 1,
      Age >= 60 & IMC >= 30 & IMC < 35 ~ 1,
      Age >= 60 & IMC >= 35 ~ 2,
      Age < 60 & IMC < 18.5 ~ 1,
      Age < 60 & IMC >= 18.5 & IMC < 25 ~ 0,
      Age < 60 & IMC >= 25 & IMC < 30 ~ 1,
      Age < 60 & IMC >= 30 & IMC < 35 ~ 1,
      Age < 60 & IMC >= 35 & IMC < 40 ~ 1,
      Age < 60 & IMC >= 40 ~ 2
    ),
    
    # Triglicerídeos
    TRIG_cat = case_when(
      is.na(TRIGLICERIDES) ~ NA_integer_,
      TRIGLICERIDES < 150 ~ 0,
      TRIGLICERIDES >= 150 & TRIGLICERIDES < 200 ~ 1,
      TRIGLICERIDES >= 200 & TRIGLICERIDES < 500 ~ 1,
      TRIGLICERIDES >= 500 ~ 2
    ),
    
    # Colesterol total em 3 categorias
    COL_cat = case_when(
      is.na(COLESTEROL) ~ NA_integer_,
      COLESTEROL >= 240 ~ 2,  # muito alto
      COLESTEROL > 200 ~ 1,   # alto
      TRUE ~ 0                # normal
    ),
    
    
    # HDL
    HDL_cat = case_when(
      is.na(HDL) ~ NA_integer_,
      HDL < 40 ~ 1,
      HDL >= 40 & HDL <= 60 ~ 0,
      HDL > 60 ~ 0
    ),
    
    # LDL
    LDL_cat = case_when(
      is.na(LDL) ~ NA_integer_,
      LDL < 100 ~ 0,
      LDL >= 100 & LDL < 130 ~ 0,
      LDL >= 130 & LDL < 160 ~ 1,
      LDL >= 160 & LDL < 190 ~ 1,
      LDL >= 190 ~ 2
    ),
    
    # TGP (ALT)
    TGP_cat = case_when(
      is.na(TGP) ~ NA_integer_,
      TGP <= 50 ~ 0,
      TGP > 50 & TGP <= 150 ~ 1,
      TGP > 150 ~ 2
    ),
    
    # TGO (AST)
    TGO_cat = case_when(
      is.na(TGO) ~ NA_integer_,
      TGO <= 50 ~ 0,
      TGO > 50 & TGO <= 150 ~ 1,
      TGO > 150 ~ 2
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


