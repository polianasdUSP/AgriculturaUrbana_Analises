#---
#  title: "BHEI-REVISED"
#output: html_notebook

  
#write.csv(metadados.all.filtrado,
#          file = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana#/Planilhas_UrbanAgri/Para_Dissertação/metadados_n105_completo.csv",
#           row.names = FALSE,
#           fileEncoding = "UTF-8")




  
  

  # Carregar pacotes
library(dplyr)
library(gt)
library(webshot2)
library(ggplot2)
library(dplyr)



# Criar o gráfico e salvar em um objeto
p_bhei <- ggplot(BHEI_R_scores130, aes(x = BHEI_R_Score_Total)) +
  geom_histogram(aes(y = ..density..), bins = 20, fill = "#619CFF", color = "black") +
  geom_density(color = "darkblue", size = 1) +
  labs(title = "Distribution of BHEI-R Scores",
       x = "BHEI-R Score (0–100)",
       y = "Density") +
  theme_minimal()

# Salvar em alta resolução
ggsave(
  filename = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/BHEI_histograma_alta.png",
  plot = p_bhei,
  width = 10, height = 7, dpi = 300
)




cores_nova <- c("1" = "#F8766D",  # vermelho
                "2" = "#00BA38",  # verde
                "3" = "#619CFF")  # azul





# Calcular estatísticas descritivas
bhei_summary <- BHEI_R_scores130 %>%
  summarise(
    Min = min(BHEI_R_Score_Total, na.rm = TRUE),
    Max = max(BHEI_R_Score_Total, na.rm = TRUE),
    Mean = mean(BHEI_R_Score_Total, na.rm = TRUE),
    SD = sd(BHEI_R_Score_Total, na.rm = TRUE)
  )

# Criar variável de tercis
BHEI_R_scores130 <- BHEI_R_scores130 %>%
  mutate(
    BHEI_Tercile = ntile(BHEI_R_Score_Total, 3),
    Tercile_Label = case_when(
      BHEI_Tercile == 1 ~ "Low",
      BHEI_Tercile == 2 ~ "Intermediate",
      BHEI_Tercile == 3 ~ "High"
    )
  )

# Obter os cortes e n por grupo
tercil_info <- BHEI_R_scores130 %>%
  group_by(Tercile_Label) %>%
  summarise(
    Min = round(min(BHEI_R_Score_Total, na.rm = TRUE), 2),
    Max = round(max(BHEI_R_Score_Total, na.rm = TRUE), 2),
    N = n()
  ) %>%
  arrange(factor(Tercile_Label, levels = c("Low", "Intermediate", "High")))

# Criar tabela com gt
tabela_gt <- gt(tercil_info) %>%
  tab_header(title = "BHEI-R Score Tertile Summary") %>%
  cols_label(
    Tercile_Label = "Diet Quality Group",
    Min = "Min",
    Max = "Max",
    N = "n"
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.align = "center"
  )

# Salvar a tabela como imagem
gtsave(
  tabela_gt,
  filename = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/BHEI_Tertile_Summary.png"
)





# Calcular média e desvio padrão
resumo_bhei <- BHEI_R_scores130 %>%
  summarise(
    Mean = round(mean(BHEI_R_Score_Total, na.rm = TRUE), 2),
    SD = round(sd(BHEI_R_Score_Total, na.rm = TRUE), 2),
    Min = round(min(BHEI_R_Score_Total, na.rm = TRUE), 2),
    Max = round(max(BHEI_R_Score_Total, na.rm = TRUE), 2),
    N = sum(!is.na(BHEI_R_Score_Total))
  )

# Transformar para gt
tabela_gt <- resumo_bhei %>%
  gt() %>%
  tab_header(
    title = "Summary of BHEI-R Total Score"
  ) %>%
  cols_label(
    Mean = "Mean",
    SD = "Standard Deviation",
    Min = "Minimum",
    Max = "Maximum",
    N = "n"
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.align = "center"
  )

# Salvar como imagem
gtsave(
  tabela_gt,
  filename = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/BHEI_Total_Summary.png"
)

library(ggplot2)
library(dplyr)

# Criar variável de tercis se ainda não tiver
BHEI_R_scores130 <- BHEI_R_scores130 %>%
  mutate(BHEI_group = ntile(BHEI_R_Score_Total, 3)) %>%
  mutate(BHEI_group = factor(BHEI_group, 
                             levels = c(1, 2, 3),
                             labels = c("Low", "Intermediate", "High")))

# Contagem por grupo
bhei_counts <- BHEI_R_scores130 %>%
  count(BHEI_group)

# Gráfico de barras
ggplot(bhei_counts, aes(x = BHEI_group, y = n, fill = BHEI_group)) +
  geom_bar(stat = "identity") +
  labs(title = "Participants by Diet Quality Terciles",
       x = "Diet Quality",
       y = "Number of Participants") +
  scale_fill_manual(values = c("Low" = "#F8766D", "Intermediate" = "#00BA38", "High" = "#619CFF")) +
  theme_minimal(base_size = 14)



#✅ Classificação sugerida para o BHEI-R (0 a 100):
#Escore BHEI-R	Classificação da Qualidade da Dieta
#≥ 60	Boa qualidade (alta)
#45 – 59,9	Qualidade intermediária
#< 45	Baixa qualidade (precisa melhorar)

#Previdelli et al. (2011) – Validação do BHEI-R

#Andrade et al. (2016) – Qualidade da dieta de adultos brasileiros

#IBGE (PNS) – Relatórios com escores de qualidade da alimentação




# 1. Classificar participantes conforme a escala real do BHEI-R
BHEI_R_scores130 <- BHEI_R_scores130 %>%
  mutate(BHEI_class = case_when(
    BHEI_R_Score_Total < 45 ~ "Low",
    BHEI_R_Score_Total >= 45 & BHEI_R_Score_Total < 60 ~ "Intermediate",
    BHEI_R_Score_Total >= 60 ~ "High"
  ))

# 2. Contar número de participantes por grupo
tabela_bhei_real <- BHEI_R_scores130 %>%
  count(BHEI_class)

# 3. Gráfico de barras
plot_bhei_class <- ggplot(BHEI_R_scores130, aes(x = BHEI_class, fill = BHEI_class)) +
  geom_bar() +
  scale_fill_manual(values = c("Low" = "#F8766D", "Intermediate" = "#00BA38", "High" = "#619CFF")) +
  labs(title = "Participants by Diet Quality (BHEI-R categories)",
       x = "Diet Quality",
       y = "Number of Participants") +
  theme_minimal(base_size = 14)

# 4. Mostrar o gráfico
print(plot_bhei_class)

# 5. Salvar em alta resolução
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/bhei_class_real.png",
       plot = plot_bhei_class, width = 8, height = 6, dpi = 300)

table(BHEI_R_scores130$BHEI_class)
