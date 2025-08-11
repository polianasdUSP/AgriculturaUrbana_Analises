PS_URBANAGRI_5013_Dieta_NOVA_07_01_2025

# Instala e carrega os pacotes necessários
library(dplyr)
library(readxl)
library(writexl)
library(fmsb)
library(ggplot2)
library(dplyr)
library(reshape2)
library(readxl)
library(tidyr)




recall <- read.delim("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/AgriculturaUrbana_Analises/urbagri_recall_nova - urbagri_recall_nova_ParaR.tsv", 
                     header = TRUE)

#tbca <- read.table("~/USP/projetos/rCodingClub/tbca_completa2024-10-04.txt", sep = "\t", header = T, row.names = 1)

# Ler a primeira aba do arquivo
tbca <- read_excel("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/Dieta/Dados_composição_TBCA_Hoffmann11022022v6.xlsx")


#Compare os valores únicos entre as colunas TBCA_code e cod_alimento para garantir que há correspondências:
unique(recall$TBCA_code)
unique(tbca$cod_alimento)

#Corrigir Erros de Formatação. Erros comuns, como espaços extras, podem ser corrigidos removendo-os:
recall$TBCA_code <- trimws(recall$TBCA_code)
tbca$cod_alimento <- trimws(tbca$cod_alimento)


# Merge dos dataframes
recall.withEnergy <- merge(
  x = recall[c("TBCA_code", "IdVoluntario", "R24H", "Quantidade", "NOVA_group", "NOVA_subgroup", "NOVA_subgroup_description_1")],
  y = tbca[c("cod_alimento", "energia2_kcal")],  # Inclua "cod_alimento" explicitamente
  by.x = "TBCA_code",
  by.y = "cod_alimento",
  all.x = TRUE  # Mantém todas as linhas de recall
)

# Verificar os primeiros dados
head(recall.withEnergy)
colnames(recall)



#convert energy based on "quantidade consumida"
recall.withEnergy$consumedCals <- recall.withEnergy$energia2_kcal*recall.withEnergy$Quantidade/100
# calculate energy per NOVA subcategory



library(dplyr)

recall.withEnergy.NOVA <- recall.withEnergy %>%
  group_by(IdVoluntario) %>% 
  mutate(nRcord = length(unique(R24H))) %>% 
  ungroup() %>% 
  group_by(IdVoluntario, R24H, NOVA_group, NOVA_subgroup, NOVA_subgroup_description_1) %>%
  summarise(consumoTotalDiario = sum(consumedCals, na.rm = TRUE), 
            nRcord = mean(nRcord, na.rm = TRUE), .groups = "drop") %>% 
  group_by(IdVoluntario, NOVA_group, NOVA_subgroup, NOVA_subgroup_description_1) %>% 
  summarise(somaConsumoDiario = sum(consumoTotalDiario, na.rm = TRUE), 
            nRcord = mean(nRcord, na.rm = TRUE), .groups = "drop") %>% 
  mutate(mediaConsumoDiario = somaConsumoDiario / nRcord) %>% 
  select(-somaConsumoDiario)  # Removendo a coluna intermediária que não é mais necessária



# calc percentages for each person 
recall.withEnergy.NOVA.percentage <- recall.withEnergy.NOVA %>%
  group_by(IdVoluntario) %>%
  mutate(consumoTotal = sum(mediaConsumoDiario)) %>%
  group_by(NOVA_subgroup, .add = TRUE) %>%  # Substituí `add` por `.add`
  mutate(prop = 100 * mediaConsumoDiario / consumoTotal)  # Calcula a proporção


#alternative
recall.withEnergy.NOVA %>%
  group_by(IdVoluntario) %>%
  mutate(per =  mediaConsumoDiario/sum(mediaConsumoDiario)) %>% 
  ungroup


#write.csv(recall.withEnergy.NOVA, 
          "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas/recall_withEnergy_NOVA.csv", 
          row.names = FALSE)

#write.table(recall.withEnergy.NOVA, 
            "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas/recall_withEnergy_NOVA.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)



#write_xlsx(recall.withEnergy.NOVA, 
           "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas/recall_withEnergy_NOVA.xlsx")


#write.csv(recall.withEnergy.NOVA.percentage, 
          "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas/recall_withEnergy_NOVA.percentage.csv", 
          row.names = FALSE)

#write.table(recall.withEnergy.NOVA.percentage, 
            "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas/recall_withEnergy_NOVA.percentage.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

#library(writexl)

#write_xlsx(recall.withEnergy.NOVA.percentage, 
#           "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas/recall_withEnergy_NOVA.percentage.xlsx")



#other transformations
recall.withEnergy.NOVA %>%
  group_by(IdVoluntario) %>%
  summarise(ConsumoDiarioTotal = sum(mediaConsumoDiario, na.rm = T))




# Calcular a porcentagem de cada NOVA_group por voluntário
nova_percentages <- recall.withEnergy.NOVA %>%
  group_by(IdVoluntario, NOVA_group) %>%
  summarise(ConsumoGrupo = sum(mediaConsumoDiario, na.rm = TRUE)) %>%
  mutate(Percentual = ConsumoGrupo / sum(ConsumoGrupo) * 100) %>%
  ungroup()

# Visualizar os primeiros registros
head(nova_percentages)

# Salvar o resultado em um arquivo Excel
#library(writexl)
#write_xlsx(nova_percentages, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas/nova_percentages_voluntario_agrUrbana.xlsx")


#====> Limpar "Nova Percentages" para ter apenas os Ids de metadados.all.filtrado
# Extrair apenas os números dos Sample.id
ids_metadados_limpos <- gsub("S(\\d{5})\\.F00", "\\1", metadados.all.filtrado$Sample.id)

# Filtrar com base nos IDs numéricos
consumo_nova_filtrado <- nova_percentages %>%
  filter(IdVoluntario %in% ids_metadados_limpos)


# verificar se estão batendo
length(unique(nova_percentages$IdVoluntario))          # Total original
length(unique(consumo_nova_filtrado$IdVoluntario)) # Total mantido

nova_percentages <- consumo_nova_filtrado



# Resumo descritivo por NOVA_group
resumo_consumo <- nova_percentages %>%
  group_by(NOVA_group) %>%
  summarise(
    MediaConsumo = mean(Percentual),
    MedianaConsumo = median(Percentual),
    DesvioPadrao = sd(Percentual)
  )
print(resumo_consumo)


# Instale se ainda não tiver
# install.packages("gt")
# webshot2::install_webshot()  # só precisa rodar uma vez

library(dplyr)
library(gt)
library(webshot2)

# Calcular estatísticas por grupo NOVA
resumo_consumo <- nova_percentages %>%
  group_by(NOVA_group) %>%
  summarise(
    Mean = mean(Percentual, na.rm = TRUE),
    Median = median(Percentual, na.rm = TRUE),
    SD = sd(Percentual, na.rm = TRUE),
    Min = min(Percentual, na.rm = TRUE),
    Max = max(Percentual, na.rm = TRUE)
  )

# Transformar NOVA_group em texto para evitar ".00"
resumo_consumo$NOVA_group <- as.character(resumo_consumo$NOVA_group)

# Criar a tabela formatada com `gt`
tabela_gt <- resumo_consumo %>%
  gt() %>%
  tab_header(
    title = "Consumption Summary by NOVA Group (n=130)"
  ) %>%
  fmt_number(
    columns = c(Mean, Median, SD, Min, Max),
    decimals = 2
  ) %>%
  cols_label(
    NOVA_group = "NOVA Group",
    Mean = "Mean (%)",
    Median = "Median (%)",
    SD = "Standard Deviation (%)",
    Min = "Minimum (%)",
    Max = "Maximum (%)"
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.align = "center"
  )

# Salvar como PNG
gtsave(tabela_gt,
       filename = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/tabela_resumo_consumo_nova_en_com_min_max.png")


# Classificar voluntários com base no maior consumo
classificacao <- nova_percentages %>%
  group_by(IdVoluntario) %>%
  slice_max(order_by = Percentual, n = 1) %>%
  select(IdVoluntario, NOVA_group, Percentual)

print(classificacao)


library(ggplot2)

# Boxplot
ggplot(nova_percentages, aes(x = as.factor(NOVA_group), y = Percentual)) +
  geom_boxplot() +
  labs(x = "NOVA Group", y = "Consumption Percentage (%)", title = "Distribution of Consumption by NOVA Group") +
  theme_minimal()


library(ggplot2)

boxplot_nova <- ggplot(nova_percentages, aes(x = as.factor(NOVA_group), y = Percentual, fill = as.factor(NOVA_group))) +
  geom_boxplot() +
  labs(
    x = "NOVA Group",
    y = "Consumption Percentage (%)",
    title = "Distribution of Consumption by NOVA Group",
    fill = "NOVA Group"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("1" = "#F8766D", "2" = "#00BA38", "3" = "#619CFF"))

# Salvar o gráfico
ggsave(
  filename = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/boxplot_nova.png",
  plot = boxplot_nova,
  width = 8,
  height = 6,
  dpi = 300
)

# Gráfico de barras das médias
# Carregar biblioteca
library(ggplot2)

# Criar o gráfico
grafico_barras <- ggplot(resumo_consumo, aes(x = as.factor(NOVA_group), y = Mean, fill = as.factor(NOVA_group))) +
  geom_bar(stat = "identity") +
  labs(x = "NOVA Group", y = "Mean Consumption (%)", title = "Mean Consumption by NOVA Group") +
  theme_minimal()

# Salvar o gráfico
ggsave(
  filename = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/grafico_consumo_medio_nova.png",
  plot = grafico_barras,
  width = 8,
  height = 6,
  dpi = 300
)




# Identificar voluntários com alto consumo de ultraprocessados
altos_consumidores <- nova_percentages %>%
  filter(NOVA_group == 3, Percentual > 50)

print(altos_consumidores)


# Identificar voluntários com alto consumo de ultraprocessados
baixos_consumidores <- nova_percentages %>%
  filter(NOVA_group == 3, Percentual < 5)

print(baixos_consumidores)

#pivotar nova_percentages



# Transformar os dados para que cada voluntário tenha apenas uma linha
nova_wide_percentagens <- nova_percentages %>%
  pivot_wider(
    id_cols = IdVoluntario, 
    names_from = NOVA_group, 
    values_from = c(ConsumoGrupo, Percentual),
    names_prefix = "NOVA_group_"
  )

# Ver a tabela final
print(nova_wide_percentagens)


# Substituir NA por 0
nova_wide_percentagens <- nova_wide_percentagens %>%
  mutate(across(starts_with("ConsumoGrupo"), ~ replace_na(., 0))) %>%
  mutate(across(starts_with("Percentual"), ~ replace_na(., 0)))

# Salvar a nova planilha
#write_xlsx(nova_wide_percentagens, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas/nova_wide_sem_NA.xlsx")

# Salvar o objeto em um arquivo CSV
#write.csv(nova_wide_percentagens, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas/nova_wide_sem_NA.csv", row.names = FALSE)




# Alterar o formato de IdVoluntario para "Sxxxxx.F00"
#nova_wide_percentagens <- nova_wide_percentagens %>%
#  mutate(IdVoluntario = paste0("S", IdVoluntario, ".F00"))

# Verificar o resultado
print(nova_wide_percentagens)

# Salvar o novo arquivo
library(writexl)
write_xlsx(nova_wide_percentagens, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas/nova_wide_percentagens_n105.xlsx")




# Realizar a junção
metadata__ASVs_Saude_Dieta_NOVA <- metadata_ASVs_Saude_dieta %>%
  left_join(nova_wide_percentagens, by = c("Sample.id" = "IdVoluntario"))

# Salvar o novo dataframe
#library(writexl)
write_xlsx(metadata__ASVs_Saude_Dieta_NOVA, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas/metadata__ASVs_Saude_Dieta_NOVA.xlsx")

# Visualizar o resultado
head(metadata_com_percentagens)



# Realizar a junção
metadados.all.filtrado <- metadados_completo %>%
  left_join(nova_wide_percentagens, by = c("Sample.id" = "IdVoluntario"))

# Salvar o novo dataframe
library(writexl)
write_xlsx(metadata_com_percentagens, "C:/Users/polia/Desktop/metadata_com_percentagens.xlsx")

colnames(metadados.all.filtrado)


#======> Comparação de Consumo de Ultraprocessados por Sexo



ggplot(metadados.all.filtrado, aes(x = Sex, y = Percentual_NOVA_group_3, fill = Sex)) +
  geom_boxplot() +
  labs(title = "Consumo de Ultraprocessados por Sexo",
       x = "Sexo",
       y = "Percentual de Ultraprocessados (%)") +
  theme_minimal()




# Teste t para comparar os grupos
t.test(Percentual_NOVA_group_3 ~ Sex, data = metadados.all.filtrado, var.equal = FALSE)


# Comparação de Consumo de Ultraprocessados por Faixas Etárias

# Criar faixas etárias
metadados.all.filtrado$Faixa_Etaria <- cut(metadados.all.filtrado$Age,
                                      breaks = c(0, 30, 50, 70, 100),
                                      labels = c("0-30", "31-50", "51-70", "71+"))

# Boxplot por Faixa Etária
ggplot(metadados.all.filtrado, aes(x = Faixa_Etaria, y = Percentual_NOVA_group_3, fill = Faixa_Etaria)) +
  geom_boxplot() +
  labs(title = "Consumo de Ultraprocessados por Faixa Etária",
       x = "Faixa Etária",
       y = "Percentual de Ultraprocessados (%)") +
  theme_minimal()



# ANOVA para testar diferenças entre grupos etários
anova_result <- aov(Percentual_NOVA_group_3 ~ Faixa_Etaria, data = metadados.all.filtrado)
summary(anova_result)



# Teste post-hoc (Tukey)
TukeyHSD(anova_result)

#Correlação de Consumo de Ultraprocessados com Idade e IMC

# Calcular IMC
metadados.all.filtrado$IMC <- metadados.all.filtrado$Weight / (metadados.all.filtrado$Height/100)^2

# Correlação entre Percentual de Ultraprocessados e Idade
cor.test(metadados.all.filtrado$Percentual_NOVA_group_3, metadados.all.filtrado$Age, method = "spearman")



# Correlação entre Percentual de Ultraprocessados e IMC
cor.test(metadados.all.filtrado$Percentual_NOVA_group_3, metadados.all.filtrado$IMC, method = "spearman")


# Scatterplot com linha de tendência
ggplot(metadados.all.filtrado, aes(x = Age, y = Percentual_NOVA_group_3)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Correlação entre Idade e Consumo de Ultraprocessados",
       x = "Idade",
       y = "Percentual de Ultraprocessados (%)") +
  theme_minimal()



# Ajustar um modelo de regressão
modelo <- lm(Percentual_NOVA_group_3 ~ Age + Sex + IMC + Region_type, data = metadados.all.filtrado)
summary(modelo)

# Visualizar os resultados
library(broom)
tidy(modelo)


library(ggplot2)

# Dados dos ajustes e resíduos
modelo_df <- data.frame(
  Ajustados = modelo$fitted.values,
  Residuos = modelo$residuals
)

# Gráfico
ggplot(modelo_df, aes(x = Ajustados, y = Residuos)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Gráfico de Resíduos vs Ajustados",
       x = "Valores Ajustados",
       y = "Resíduos") +
  theme_minimal()

library(ggplot2)
library(ggfortify)

# Gráfico de efeitos parciais
autoplot(modelo, which = 1:6, ncol = 2)

# Mostrando a linha 146 do dataframe
metadados.all.filtrado[146, ]



#scatterplot consumo de UPF e idade

ggplot(metadados.all.filtrado, aes(x = Age, y = Percentual_NOVA_group_3)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(title = "Consumo de Ultraprocessados vs Idade",
       x = "Idade",
       y = "Percentual de Ultraprocessados (%)") +
  theme_minimal()


