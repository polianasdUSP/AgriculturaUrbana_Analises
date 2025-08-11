
#==================================#
# "Dieta e Agricultura Urbana"
#==================================#


# 2024-05-23
# dietary data analysis code V01

# using the tbca and recall data. 

# tbca is standardized within the lab and by the TBCA 
# the recall data must have 4 columns:

# IdVoluntario, R24H, Quantidade, TBCA_code
# must be sequential for all volunteers (one on top of the other):
# must have these colnames

# as such:
# IdVoluntario	R24H	Quantidade	TBCA_code
# 10011	Day 1	200	C0409A
# 10011	Day 2	100	C0044H
# 10011	Day 3	225	C0064T
# 10012	Day 1	75	C0209A
# 10012	Day 2	20	C0145B
# 10012	Day 3	11	C0170F

# these 4 cols have the person id, the day the recall was collected, the food item code and the quantity ingested



# Instala e carrega os pacotes necessários
library(dplyr)
library(readxl)
library(writexl)
library(fmsb)
library(ggplot2)
library(dplyr)
library(reshape2)



#========== analysis ==============#


#import the TBCA

tbca.raw <- read_excel("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/Dieta/Dados_composição_TBCA_Hoffmann11022022v6.xlsx")

# convert negative values that indicate traces of not available data etc to zero
tbca <- tbca.raw %>%
  mutate(across(where(is.numeric), ~ replace(., . < 0, 0)))



# import the recall data

recall.raw <- read_excel("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/Dieta/Planilha_UNICA_R24H_Projeto_Agricultura (3).xlsx")
# convert ID to char
recall.raw$IdVoluntario <- as.character(recall.raw$IdVoluntario)

recall <- merge(x = recall.raw[c("IdVoluntario", "R24H", "Quantidade", "TBCA_code")], y = tbca, by.x = "TBCA_code", by.y = "cod_alimento", all.x =T)
recall$convertionFactor <- recall$Quantidade/100

recall <- cbind(recall[1:5], recall[12:(ncol(recall)-1)]*recall$convertionFactor)

# Salvar o computador
write.table(recall, file = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/Dieta/recall.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write_xlsx(recall, path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/Dieta/recall.xlsx")

# Salvar o computador
write.table(saude_dieta, file = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/Dieta/saude_dieta.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write_xlsx(saude_dieta, path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/Dieta/saude_dieta.xlsx")
str(recall)



recall.calc <- recall[-c(1,4,5)] %>% 
  group_by(IdVoluntario, R24H) %>%
  summarise_all(sum, na.rm=T) %>%
  summarise(across(-R24H, mean))

#colocar S no inicio do ID e .F00 no final
recall.calc$IdVoluntario <- paste0("S", recall.calc$IdVoluntario, ".F00")



write.table(x = recall.calc, file = "recall.calc.txt", col.names = T)

write.table(recall.calc, file = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/recall.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Salva o objeto recall.calc em um arquivo Excel
write_xlsx(recall.calc, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/Dieta/recall_calc.xlsx")


hist(recall.calc$energia2_kcal, main="Distribuição de Energia (kcal)", xlab="Energia (kcal)", col="blue", border="black")


boxplot(recall.calc$proteina_g ~ recall.calc$IdVoluntario, main="Distribuição de Proteína por Voluntário", xlab="ID do Voluntário", ylab="Proteína (g)", col="orange")

barplot(colMeans(recall.calc[,c("proteina_g", "carboidrato_total_g", "lipidios_g")]), 
        main="Média de Macronutrientes Consumidos", 
        ylab="Quantidade (g)", 
        col=c("red", "green", "blue"), 
        names.arg=c("Proteína", "Carboidrato", "Lipídios"))

plot(recall.calc$carboidrato_disponivel_g, recall.calc$energia2_kcal, 
     main="Relação entre Carboidratos Disponíveis e Energia (kcal)", 
     xlab="Carboidratos Disponíveis (g)", ylab="Energia (kcal)", 
     col="purple", pch=19)

plot(recall.calc$IdVoluntario, recall.calc$vitamina_C_mg, type="o", col="blue", 
     main="Tendência de Consumo de Vitamina C", xlab="ID do Voluntário", ylab="Vitamina C (mg)")


df_radar <- recall.calc[,c("proteina_g", "carboidrato_total_g", "lipidios_g")]
df_radar <- rbind(rep(max(df_radar, na.rm = TRUE),3), rep(min(df_radar, na.rm = TRUE),3), df_radar)
radarchart(df_radar, axistype=1, 
           pcol=rainbow(3), pfcol=rainbow(3, alpha=0.5), plwd=2)

# Supondo que recall.calc e saude_raw sejam data frames


# Verificar as mudanças
head(recall.calc$IdVoluntario)


# Junta os data frames utilizando os IDs dos voluntários
saude_dieta <- merge(recall.calc, metadados, by.x = "IdVoluntario", by.y = "Sample.id")

# Exibe as primeiras linhas do data frame resultante
head(saude_dieta)

write_xlsx(saude_dieta, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/saude_dieta130.xlsx")

# Adiciona uma coluna indicando se o consumo de fibra é suficiente
saude_dieta$Consumo_Fibra <- ifelse(saude_dieta$fibra_alimentar_g.x > 25, "Suficiente", "Insuficiente")

saude_dieta$Consumo_Fibra




#=======================================#
#       Consumo de Fibras
#=======================================#


#Contar quantos são suficientes vs. insuficientes
table(saude_dieta$Consumo_Fibra)

library(ggplot2)

#Gráfico de barras
ggplot(saude_dieta, aes(x = Consumo_Fibra)) +
  geom_bar(fill = "steelblue") +
  labs(title = "Fiber Intake", x = "Category", y = "Number of Participants") +
  theme_minimal()

# Ver proporção (em %)
prop.table(table(saude_dieta$Consumo_Fibra)) * 100

#Ver quem são os voluntários com consumo suficiente/insuficiente
saude_dieta %>% 
  select(IdVoluntario, fibra_alimentar_g.x, Consumo_Fibra) %>%
  arrange(desc(Consumo_Fibra))

#Ver distribuição média dos grupos
library(dplyr)

saude_dieta %>%
  group_by(Consumo_Fibra) %>%
  summarise(
    media_fibra = mean(fibra_alimentar_g.x),
    sd_fibra = sd(fibra_alimentar_g),
    n = n()
  )


#Boxplot para comparar os grupos
library(ggplot2)

ggplot(saude_dieta, aes(x = Consumo_Fibra, y = fibra_alimentar_g.x, fill = Consumo_Fibra)) +
  geom_boxplot() +
  labs(title = "Fiber intake vs Category",
       x = "Category",
       y = "Fiber (g)") +
  theme_minimal()




#1. Tabela com média, desvio padrão, mínimo e máximo do consumo de fibras
# Step 1: Load required package
library(dplyr)


#Deixar apenas 105 voluntarios



# Step 2: Create summary table grouped by fiber intake adequacy
fiber_summary <- saude_dieta %>%
  group_by(Consumo_Fibra) %>%  # Groups by "Sufficient" and "Insufficient"
  summarise(
    n = n(),  # Number of participants in each group
    mean = round(mean(fibra_alimentar_g.x), 2),  # Mean fiber intake
    sd = round(sd(fibra_alimentar_g.x), 2),  # Standard deviation
    min = round(min(fibra_alimentar_g.x), 2),  # Minimum intake
    max = round(max(fibra_alimentar_g.x), 2)  # Maximum intake
  )

# Step 3: View the result
fiber_summary


#Gerar a tabela estilizada com gt

# Instalar se ainda não tiver
# install.packages("gt")
library(dplyr)
library(gt)



# Tabela com total
resumo_por_categoria <- fiber_summary %>%
  rename(
    `Fiber Intake Category` = Consumo_Fibra,
    `n` = n,
    `Mean (g)` = mean,
    `SD (g)` = sd,
    `Min (g)` = min,
    `Max (g)` = max
  )

resumo_total <- saude_dieta %>%
  summarise(
    `Fiber Intake Category` = "Total",
    n = n(),
    `Mean (g)` = round(mean(fibra_alimentar_g.x), 2),
    `SD (g)` = round(sd(fibra_alimentar_g.x), 2),
    `Min (g)` = round(min(fibra_alimentar_g.x), 2),
    `Max (g)` = round(max(fibra_alimentar_g.x), 2)
  )

tabela_com_total <- bind_rows(resumo_por_categoria, resumo_total)

# Gerar tabela com o primeiro alinhamento à esquerda
tabela_gt <- tabela_com_total %>%
  gt() %>%
  tab_header(
    title = md("**Descriptive Statistics for Fiber Intake**")
  ) %>%
  cols_align(align = "left", columns = `Fiber Intake Category`) %>%  # Alinha a primeira coluna à esquerda
  cols_align(align = "center", columns = c(n, `Mean (g)`, `SD (g)`, `Min (g)`, `Max (g)`)) %>%  # Alinha o resto ao centro
  tab_options(
    table.font.size = 12,
    heading.title.font.size = 14,
    heading.title.font.weight = "bold"
  )

# Salvar a tabela como PNG
gtsave(tabela_gt,
       filename = "fiber_summary_with_total_aligned105.png",
       path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/")


#2. Pie chart showing proportion of participants by fiber adequacy

# Step 1: Load required package
library(ggplot2)

# Step 1: Translate categories to English
saude_dieta$Fiber_Status <- ifelse(saude_dieta$Consumo_Fibra == "Suficiente", 
                                            "Adequate", "Inadequate")

# Step 2: Create proportion table
fiber_pie <- saude_dieta %>%
  count(Fiber_Status) %>%
  mutate(
    prop = n / sum(n) * 100,
    label = paste0(Fiber_Status, "\n", round(prop, 1), "%")
  )

# Step 3: Create the pie chart
fiber_pie_chart <- ggplot(fiber_pie, aes(x = "", y = prop, fill = Fiber_Status)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 5) +
  labs(
    title = "Distribution of Fiber Intake Adequacy",
    fill = "Fiber Intake"
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

# Step 4: Display the chart
print(fiber_pie_chart)

# Step 5: Save in high resolution
# Salvar o gráfico em inglês em alta resolução na pasta da dissertação
ggsave("fiber_pie_chart_english105.png",
       plot = fiber_pie_chart,
       dpi = 300,
       width = 6,
       height = 6,
       path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/")



# Carrega os data frames, se necessário
# recall.calc <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/Dieta/recall_calc.csv")

# saude_raw <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/Dieta/saude_raw.csv")

# Junta os data frames mantendo todas as linhas de recall.calc
saude_dieta <- merge(recall.calc, saude_raw, by.x = "IdVoluntario", by.y = "Id", all.x = TRUE)

# Exibe as primeiras linhas do data frame resultante
head(saude_dieta)

# Cria o gráfico de barras empilhadas

# Recode values in Sex column
saude_dieta$Sex <- recode(saude_dieta$Sex,
                                   "Masculino" = "Male",
                                   "Feminino" = "Female")


# Criar o gráfico e salvar em um objeto
grafico_fibra_sexo <- ggplot(saude_dieta, aes(x = Sex, fill = Consumo_Fibra)) +
  geom_bar(position = "stack") +
  labs(title = "Fiber Intake",
       x = "Sex",
       y = "Number of Volunteers",
       fill = "Fiber Intake") +
  theme_minimal()


# Salvar o gráfico como PNG (300 dpi)
ggsave(filename = "fiber_intake_by_sex105.png",
       plot = grafico_fibra_sexo,
       dpi = 300,
       width = 6,
       height = 5,
       path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/")




library(dplyr)

# Tabela com número de voluntários por sexo e consumo de fibra
tabela_fibra_sexo <- saude_dieta %>%
  group_by(Sex, Consumo_Fibra) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Sex) %>%
  mutate(
    percent = round(n / sum(n) * 100, 1)
  )

# Visualizar
tabela_fibra_sexo


# Calculando a média das colunas desejadas
medias <- colMeans(saude_dieta[, c("fibra_alimentar_g", "lipidios_g", "proteina_g", "carboidrato_total_g")], na.rm = TRUE)

# Convertendo o resultado em um data frame para facilitar a visualização como tabela
tabela_medias <- data.frame(Nutriente = names(medias), Media = medias)

# Exibindo a tabela
print(tabela_medias)


# Calcula as proporções por sexo

fibra_por_sexo <- saude_dieta %>%
  group_by(Sex, Consumo_Fibra) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Percent = Count / sum(Count) * 100)

# Visualizando o resultado
print(fibra_por_sexo)

# Calculando a porcentagem correta para cada grupo de sexo
fibra_por_sex <- saude_dieta %>%
  filter(!is.na(Sex) & !is.na(Consumo_Fibra)) %>%
  group_by(Sex, Consumo_Fibra) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  group_by(Sex) %>%
  mutate(Percent = Percent / sum(Percent) * 100)

# Criando o gráfico de pizza com as porcentagens dentro das fatias
ggplot(fibra_por_sexo, aes(x = "", y = Percent, fill = Consumo_Fibra)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  facet_wrap(~Sex) +
  geom_text(aes(label = paste0(round(Percent, 1), "%")), 
            position = position_stack(vjust = 0.5), 
            color = "white") +  # Cor do texto para ficar visível sobre as fatias
  labs(title = "Consumo Suficiente de Fibra por Sexo",
       x = "",
       y = "Percentagem",
       fill = "Consumo de Fibra") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        panel.grid  = element_blank())



# Verifica a estrutura dos dados
print(fibra_por_sexo)

# Verificar a estrutura dos dados e amostras
str(saude_dieta_filtrado)
head(saude_dieta_filtrado)

# Remove valores NA das colunas Sex e fibra_alimentar_g
saude_dieta_filtrado <- saude_dieta_filtrado %>%
  filter(!is.na(Sex) & !is.na(fibra_alimentar_g))

summary(saude_dieta_filtrado)







