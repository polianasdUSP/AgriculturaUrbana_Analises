
#==============================================# 
#title: "Summary, Plots and Descriptive Stats"
#==============================================#


#===============================# 
# 1. Install necessary packages
#===============================#
# Pacotes necessários

library(dplyr)
library(gt)
install.packages("pagedown")
install.packages("xfun")
library(pagedown)
install.packages("FSA")
library(FSA)
library(ggplot2)

# Carregando os dados

metadados <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados.raw.csv") # Atualize com o caminho correto do arquivo

#criar uma coluna com IMC

# Adicionando uma coluna de IMC ao metadados
metadados <- metadados %>%
  mutate(
    IMC = Weight / (Height^2) # Cálculo do IMC
  )

# Visualizando as primeiras linhas para verificar o resultado
head(metadados$IMC)

#salvar planilha

write.csv(x = metadados, 
          file = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados.csv", row.names = FALSE)


# Salvando como CSV
write.csv(caracteristicas, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/population_characteristics_corrected.csv", row.names = FALSE)

# Calculando características da população (Mean, SD, Min e Max)
caracteristicas <- data.frame(
  Variable = c(
    "Age (years)", "BMI (Kg/m²)", "Waist-to-hip ratio", "Fasting glucose (mmol/L)",
    "Fasting insulin (μU/mL)", "HbA1c (%)", "Total cholesterol (mmol/L)",
    "LDL-C (mmol/L)", "HDL-C (mmol/L)", "Triglycerides (mmol/L)",
    "Systolic blood pressure (mmHg)", "Diastolic blood pressure (mmHg)",
    "CRP (mg/L)", "TyG", "HOMA-IR"
  ),
  Mean = c(
    round(mean(metadados$Age, na.rm = TRUE), 2),
    round(mean(metadados$IMC, na.rm = TRUE), 2),
    round(mean(metadados$Waist / metadados$Hip, na.rm = TRUE), 2),
    round(mean(metadados$GLICOSE, na.rm = TRUE), 2),
    round(mean(metadados$INSULINA, na.rm = TRUE), 2),
    round(mean(metadados$HbA1c, na.rm = TRUE), 2),
    round(mean(metadados$COLESTEROL, na.rm = TRUE), 2),
    round(mean(metadados$LDL, na.rm = TRUE), 2),
    round(mean(metadados$HDL, na.rm = TRUE), 2),
    round(mean(metadados$TRIGLICERIDES, na.rm = TRUE), 2),
    round(mean(metadados$Systolic, na.rm = TRUE), 2),
    round(mean(metadados$Diastolic, na.rm = TRUE), 2),
    round(mean(metadados$PCR, na.rm = TRUE), 2),
    round(mean(metadados$TyG, na.rm = TRUE), 2),
    round(mean(metadados$HOMA, na.rm = TRUE), 2)
  ),
  SD = c(
    round(sd(metadados$Age, na.rm = TRUE), 2),
    round(sd(metadados$IMC, na.rm = TRUE), 2),
    round(sd(metadados$Waist / metadados$Hip, na.rm = TRUE), 2),
    round(sd(metadados$GLICOSE, na.rm = TRUE), 2),
    round(sd(metadados$INSULINA, na.rm = TRUE), 2),
    round(sd(metadados$HbA1c, na.rm = TRUE), 2),
    round(sd(metadados$COLESTEROL, na.rm = TRUE), 2),
    round(sd(metadados$LDL, na.rm = TRUE), 2),
    round(sd(metadados$HDL, na.rm = TRUE), 2),
    round(sd(metadados$TRIGLICERIDES, na.rm = TRUE), 2),
    round(sd(metadados$Systolic, na.rm = TRUE), 2),
    round(sd(metadados$Diastolic, na.rm = TRUE), 2),
    round(sd(metadados$PCR, na.rm = TRUE), 2),
    round(sd(metadados$TyG, na.rm = TRUE), 2),
    round(sd(metadados$HOMA, na.rm = TRUE), 2)
  ),
  Min = c(
    round(min(metadados$Age, na.rm = TRUE), 2),
    round(min(metadados$IMC, na.rm = TRUE), 2),
    round(min(metadados$Waist / metadados$Hip, na.rm = TRUE), 2),
    round(min(metadados$GLICOSE, na.rm = TRUE), 2),
    round(min(metadados$INSULINA, na.rm = TRUE), 2),
    round(min(metadados$HbA1c, na.rm = TRUE), 2),
    round(min(metadados$COLESTEROL, na.rm = TRUE), 2),
    round(min(metadados$LDL, na.rm = TRUE), 2),
    round(min(metadados$HDL, na.rm = TRUE), 2),
    round(min(metadados$TRIGLICERIDES, na.rm = TRUE), 2),
    round(min(metadados$Systolic, na.rm = TRUE), 2),
    round(min(metadados$Diastolic, na.rm = TRUE), 2),
    round(min(metadados$PCR, na.rm = TRUE), 2),
    round(min(metadados$TyG, na.rm = TRUE), 2),
    round(min(metadados$HOMA, na.rm = TRUE), 2)
  ),
  Max = c(
    round(max(metadados$Age, na.rm = TRUE), 2),
    round(max(metadados$IMC, na.rm = TRUE), 2),
    round(max(metadados$Waist / metadados$Hip, na.rm = TRUE), 2),
    round(max(metadados$GLICOSE, na.rm = TRUE), 2),
    round(max(metadados$INSULINA, na.rm = TRUE), 2),
    round(max(metadados$HbA1c, na.rm = TRUE), 2),
    round(max(metadados$COLESTEROL, na.rm = TRUE), 2),
    round(max(metadados$LDL, na.rm = TRUE), 2),
    round(max(metadados$HDL, na.rm = TRUE), 2),
    round(max(metadados$TRIGLICERIDES, na.rm = TRUE), 2),
    round(max(metadados$Systolic, na.rm = TRUE), 2),
    round(max(metadados$Diastolic, na.rm = TRUE), 2),
    round(max(metadados$PCR, na.rm = TRUE), 2),
    round(max(metadados$TyG, na.rm = TRUE), 2),
    round(max(metadados$HOMA, na.rm = TRUE), 2)
  )
)



# Visualizando o resultado final
print(caracteristicas)

# Salvando como CSV
write.csv(caracteristicas, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/population_characteristics_final.csv", row.names = FALSE)



#===========================================#
#Fazer uma tabela com formatação bonita
#===========================================#

library(gt)

# Criar a tabela bonitinha com gt
TableClinicalTraits <- caracteristicas %>%
  gt() %>%
  tab_header(
    title = "Clinical Traits",
    subtitle = "(n=130)"
  ) %>%
  fmt_number(
    columns = c(Mean, SD),
    decimals = 2
  ) %>%
  cols_label(
    Variable = "Variable",
    Mean = "Mean",
    SD = "Standard Deviation"
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.align = "center"
  )

# Exibir a tabela
TableClinicalTraits

install.packages("pagedown")
install.packages("xfun")
library(pagedown)

gtsave(TableClinicalTraits, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/TableClinicalTraits.html")

#========================================#
#Another option for a pretty table
#========================================#

install.packages("kableExtra")
library(kableExtra)

# Criando a tabela com kableExtra
tabela_características_população2 <- caracteristicas %>%
  knitr::kable("html", col.names = c("Variable", "Mean", "Standard Deviation")) %>%
  kableExtra::kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "responsive"),
    full_width = FALSE,
    font_size = 14
  )

# Exibir a tabela no Viewer
tabela_características_população2

# Salvar a tabela como HTML
tabela_características_população2 %>%
  save_kable("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/population_characteristics_table2.html")

#=========================================#
     #Table with AGE, BMI and sex
#=========================================#


#Kruskal-Wallis para mais de dois grupos
#Para comparar IMC ou Idade entre Sul, Norte, Leste, etc.:


# Kruskal-Wallis para IMC
kruskal.test(IMC ~ Region, data = metadados)


# Kruskal-Wallis para Idade


kruskal.test(Age ~ Region, data = metadados)



####Se o resultado do Kruskal-Wallis for significativo, você pode fazer comparações par-a-par usando Dunn's Test:###

install.packages("FSA")
library(FSA)

# Comparações par-a-par
dunnTest(IMC ~ Region, data = metadados, method = "bonferroni")


######----------------------------------------------####
#Mann-Whitney U para dois grupos (Rural vs Urbano)
######----------------------------------------------####


#Para comparar IMC ou Idade entre Rural e Urbano:

# Mann-Whitney para IMC
wilcox.test(IMC ~ Region_type, data = metadados)

# Mann-Whitney para Idade
wilcox.test(Age ~ Region_type, data = metadados)

######----------------------------------------------------####
#Qui-quadrado ou Fisher para variáveis categóricas (Sexo)
######----------------------------------------------------####

#Se você quiser comparar Sexo entre regiões ou tipos de região:

######Para mais de dois grupos (Sul, Norte, etc.):######

# Tabela cruzada para Sexo e Região
tabela_sexo <- table(metadados$Sex, metadados$Region)

# Qui-quadrado
chisq.test(tabela_sexo)


#Para dois grupos (Rural vs Urbano):

# Tabela cruzada para Sexo e Tipo de Região
tabela_sexo <- table(metadados$Sex, metadados$Region_type)

# Teste exato de Fisher (para pequenas amostras)
fisher.test(tabela_sexo)

###==========================================###
### IMPORTANTE PARA ESCOLHER OS TESTES ###
###===========================================###

#Resumo de Como Escolher o Teste

#Variável: IMC, Idade	Grupos: 2 grupos (Rural/Urbano)	=> Mann-Whitney U

#Variável: IMC, Idade	Grupos: Mais de 2 grupos =>	Kruskal-Wallis

#Variável: Sexo (categórica)	Grupos: 2 grupos (Rural/Urbano)	=> Fisher ou Qui-quadrado

#Variável: Sexo (categórica)	Grupos: Mais de 2 grupos	=> Qui-quadrado


#### !!!!! Não deu diferença estatística !!!!! ####
#O que fazer a seguir?
#Verificar os dados:

#Certifique-se de que os dados de idade estão corretos e completos.
#Visualize os dados por região usando boxplots para identificar diferenças visuais ou outliers.

#================================#
#Exemplo de código para boxplots:
#================================#


library(ggplot2)

#Alterar a paleta de cores. Se quiser mudar as cores, substitua scale_fill_brewer(palette = "Set3") por outra paleta disponível, como: "Paired", "Pastel1", "Dark2" Ou use cores personalizadas: 

boxplot_idade <- ggplot(metadados, aes(x = Region, y = Age, fill = Region)) +
  geom_boxplot() +
  labs(title = "Distribuição da idade por região",
       x = "Região",
       y = "Idade") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")  # Paleta de cores

# Exibir o gráfico
boxplot_idade

# Salvar como PNG
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/boxplot_idade.png", plot = boxplot_idade, width = 8, height = 6, dpi = 300)

# Salvar como PDF
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/boxplot_idade.pdf", plot = boxplot_idade, width = 8, height = 6)

#======================================================================#
#Tabela com : Idade, Sexo, Média de Idade entre os sexos e Número Total. 
#Criar a tabela resumida
#Usaremos o dplyr para calcular os resumos e o gt para formatar a tabela.
#======================================================================#

# Carregar os pacotes necessários
library(dplyr)
library(gt)

# Ajustar os rótulos para inglês (Male e Female)
tabela_resumo <- metadados %>%
  mutate(Sex = ifelse(Sex == "Feminino", "Female", 
                      ifelse(Sex == "Masculino", "Male", Sex))) %>%
  group_by(Sex) %>%
  summarise(
    `N` = n(),  # Contagem total por sexo
    `Mean Age` = round(mean(Age, na.rm = TRUE), 2),  # Média de idade
    `SD` = round(sd(Age, na.rm = TRUE), 2),  # Desvio padrão da idade
    `Min Age` = min(Age, na.rm = TRUE),  # Idade mínima
    `Max Age` = max(Age, na.rm = TRUE)   # Idade máxima
  ) %>%
  ungroup()

# Criar uma tabela formatada com gt
tabela_final <- tabela_resumo %>%
  gt() %>%
  tab_header(
    title = "Age and Sex (n=130)"
  ) %>%
  fmt_number(
    columns = c(`Mean Age`, `SD`, `Min Age`, `Max Age`),
    decimals = 2
  ) %>%
  cols_label(
    Sex = "Sex",
    `N` = "Count",
    `Mean Age` = "Mean Age",
    `SD` = "Standard Deviation",
    `Min Age` = "Minimum Age",
    `Max Age` = "Maximum Age"
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.align = "center"
  )

# Exibir a tabela
tabela_final

# Salvar como HTML
gtsave(tabela_final, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/age_sex_summary.html")

# Salvar como PDF
gtsave(tabela_final, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/age_sex_summary.jpg")

gtsave(
  tabela_final, 
  "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/age_sex_summary.png"
)



### Parasitológico. ###
#Tabela que mostre quantos voluntários têm resultado positivo no parasitológico e listar os tipos de parasitas encontrados


# Carregar os pacotes necessários
library(dplyr)
library(gt)

# Substituir "e" por "and" na coluna Parasitas
#metadados <- metadados %>%
#  mutate(Parasitas = gsub(" e ", " and ", Parasitas))

# Substituir "Nana" por "nana" na coluna Parasitas
#metadados <- metadados %>%
#  mutate(Parasitas = gsub("Nana", "nana", Parasitas))

# Fazer as substituições na coluna Parasitas
#metadados <- metadados %>%
#  mutate(
#    Parasitas = gsub("Prersença", "Presença", Parasitas)  # Corrigir "Prersença" para "Presença"
#  )

# Substituir "/" por vírgula e ajustar o formato
#metadados <- metadados %>%
#  mutate(
#   Parasitas = gsub(" . ", ", ", Parasitas)  # Substituir " / " por ", "
#  )

# Verificar se as alterações foram realizadas corretamente
unique(metadados$Parasitas)

# Contar os resultados do parasitológico
resumo_parasitologico <- metadados %>%
  group_by(ParasitologicoPositivo) %>%
  summarise(N = n()) %>%
  arrange(desc(N))

# Exibir o resumo
print(resumo_parasitologico)


# Substituir "e" por "and" na coluna de tipos de parasitas
#tabela_parasitologico <- metadados %>%
#  filter(ParasitologicoPositivo == "Sim") %>%  # Selecionar apenas positivos no parasitológico
#  group_by(Parasitas) %>%  # Agrupar por tipos de parasitas
#  summarise(`Número de casos` = n()) %>%  # Contar ocorrências de cada parasita
#  arrange(desc(`Número de casos`)) %>%  # Ordenar por número de casos
#  ungroup()

# Criar uma tabela formatada com gt
tabela_final_parasitas <- tabela_parasitologico %>%
  gt() %>%
  tab_header(
    title = "Parasitological Test (n=130) ("
  ) %>%
  cols_label(
    Parasitas = "Parasite Type",
    `Número de casos` = "Number of Cases"
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.align = "center"
  )

# Exibir a tabela
tabela_final_parasitas



# Carregando os dados
metadados <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados.raw.csv") # Atualize com o caminho correto do arquivo

#criar uma coluna com IMC

# Adicionando uma coluna de IMC ao metadados
#metadados <- metadados %>%
#  mutate(
#    IMC = Weight / (Height^2) # Cálculo do IMC
#  )

# Visualizando as primeiras linhas para verificar o resultado
head(metadados$IMC)

write.csv(x = metadados, 
          file = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados.csv", row.names = FALSE)


# Salvando como CSV
write.csv(caracteristicas, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/population_characteristics_corrected.csv", row.names = FALSE)

# Calculando características da população (Mean, SD, Min e Max)
caracteristicas <- data.frame(
  Variable = c(
    "Age (years)", "BMI (Kg/m²)", "Waist-to-hip ratio", "Fasting glucose (mmol/L)",
    "Fasting insulin (μU/mL)", "HbA1c (%)", "Total cholesterol (mmol/L)",
    "LDL-C (mmol/L)", "HDL-C (mmol/L)", "Triglycerides (mmol/L)",
    "Systolic blood pressure (mmHg)", "Diastolic blood pressure (mmHg)",
    "CRP (mg/L)", "TyG", "HOMA-IR"
  ),
  Mean = c(
    round(mean(metadados$Age, na.rm = TRUE), 2),
    round(mean(metadados$IMC, na.rm = TRUE), 2),
    round(mean(metadados$Waist / metadados$Hip, na.rm = TRUE), 2),
    round(mean(metadados$GLICOSE, na.rm = TRUE), 2),
    round(mean(metadados$INSULINA, na.rm = TRUE), 2),
    round(mean(metadados$HbA1c, na.rm = TRUE), 2),
    round(mean(metadados$COLESTEROL, na.rm = TRUE), 2),
    round(mean(metadados$LDL, na.rm = TRUE), 2),
    round(mean(metadados$HDL, na.rm = TRUE), 2),
    round(mean(metadados$TRIGLICERIDES, na.rm = TRUE), 2),
    round(mean(metadados$Systolic, na.rm = TRUE), 2),
    round(mean(metadados$Diastolic, na.rm = TRUE), 2),
    round(mean(metadados$PCR, na.rm = TRUE), 2),
    round(mean(metadados$TyG, na.rm = TRUE), 2),
    round(mean(metadados$HOMA, na.rm = TRUE), 2)
  ),
  SD = c(
    round(sd(metadados$Age, na.rm = TRUE), 2),
    round(sd(metadados$IMC, na.rm = TRUE), 2),
    round(sd(metadados$Waist / metadados$Hip, na.rm = TRUE), 2),
    round(sd(metadados$GLICOSE, na.rm = TRUE), 2),
    round(sd(metadados$INSULINA, na.rm = TRUE), 2),
    round(sd(metadados$HbA1c, na.rm = TRUE), 2),
    round(sd(metadados$COLESTEROL, na.rm = TRUE), 2),
    round(sd(metadados$LDL, na.rm = TRUE), 2),
    round(sd(metadados$HDL, na.rm = TRUE), 2),
    round(sd(metadados$TRIGLICERIDES, na.rm = TRUE), 2),
    round(sd(metadados$Systolic, na.rm = TRUE), 2),
    round(sd(metadados$Diastolic, na.rm = TRUE), 2),
    round(sd(metadados$PCR, na.rm = TRUE), 2),
    round(sd(metadados$TyG, na.rm = TRUE), 2),
    round(sd(metadados$HOMA, na.rm = TRUE), 2)
  ),
  Min = c(
    round(min(metadados$Age, na.rm = TRUE), 2),
    round(min(metadados$IMC, na.rm = TRUE), 2),
    round(min(metadados$Waist / metadados$Hip, na.rm = TRUE), 2),
    round(min(metadados$GLICOSE, na.rm = TRUE), 2),
    round(min(metadados$INSULINA, na.rm = TRUE), 2),
    round(min(metadados$HbA1c, na.rm = TRUE), 2),
    round(min(metadados$COLESTEROL, na.rm = TRUE), 2),
    round(min(metadados$LDL, na.rm = TRUE), 2),
    round(min(metadados$HDL, na.rm = TRUE), 2),
    round(min(metadados$TRIGLICERIDES, na.rm = TRUE), 2),
    round(min(metadados$Systolic, na.rm = TRUE), 2),
    round(min(metadados$Diastolic, na.rm = TRUE), 2),
    round(min(metadados$PCR, na.rm = TRUE), 2),
    round(min(metadados$TyG, na.rm = TRUE), 2),
    round(min(metadados$HOMA, na.rm = TRUE), 2)
  ),
  Max = c(
    round(max(metadados$Age, na.rm = TRUE), 2),
    round(max(metadados$IMC, na.rm = TRUE), 2),
    round(max(metadados$Waist / metadados$Hip, na.rm = TRUE), 2),
    round(max(metadados$GLICOSE, na.rm = TRUE), 2),
    round(max(metadados$INSULINA, na.rm = TRUE), 2),
    round(max(metadados$HbA1c, na.rm = TRUE), 2),
    round(max(metadados$COLESTEROL, na.rm = TRUE), 2),
    round(max(metadados$LDL, na.rm = TRUE), 2),
    round(max(metadados$HDL, na.rm = TRUE), 2),
    round(max(metadados$TRIGLICERIDES, na.rm = TRUE), 2),
    round(max(metadados$Systolic, na.rm = TRUE), 2),
    round(max(metadados$Diastolic, na.rm = TRUE), 2),
    round(max(metadados$PCR, na.rm = TRUE), 2),
    round(max(metadados$TyG, na.rm = TRUE), 2),
    round(max(metadados$HOMA, na.rm = TRUE), 2)
  )
)



# Visualizando o resultado final
print(caracteristicas)

# Salvando como CSV
write.csv(caracteristicas, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/population_characteristics_final.csv", row.names = FALSE)



#####Fazer uma tabela com formatação bonita######

install.packages("gt")

library(gt)

# Criar a tabela bonitinha com gt
TableClinicalTraits <- caracteristicas %>%
  gt() %>%
  tab_header(
    title = "Clinical Traits",
    subtitle = "(n=130)"
  ) %>%
  fmt_number(
    columns = c(Mean, SD),
    decimals = 2
  ) %>%
  cols_label(
    Variable = "Variable",
    Mean = "Mean",
    SD = "Standard Deviation"
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.align = "center"
  )

# Exibir a tabela
TableClinicalTraits

install.packages("pagedown")
install.packages("xfun")
library(pagedown)

gtsave(TableClinicalTraits, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/TableClinicalTraits.html")




install.packages("kableExtra")
library(kableExtra)

# Criando a tabela com kableExtra
tabela_características_população2 <- caracteristicas %>%
  knitr::kable("html", col.names = c("Variable", "Mean", "Standard Deviation")) %>%
  kableExtra::kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "responsive"),
    full_width = FALSE,
    font_size = 14
  )

# Exibir a tabela no Viewer
tabela_características_população2

# Salvar a tabela como HTML
tabela_características_população2 %>%
  save_kable("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/population_characteristics_table2.html")

#Tabela com idade, IMC e sexo

```{r}
#Kruskal-Wallis para mais de dois grupos
#Para comparar IMC ou Idade entre Sul, Norte, Leste, etc.:


# Kruskal-Wallis para IMC
kruskal.test(IMC ~ Region, data = metadados)

# Kruskal-Wallis para Idade


kruskal.test(Age ~ Region, data = metadados)



####Se o resultado do Kruskal-Wallis for significativo, você pode fazer comparações par-a-par usando Dunn's Test:#####


install.packages("FSA")
library(FSA)

# Comparações par-a-par
dunnTest(IMC ~ Region, data = metadados, method = "bonferroni")


#Mann-Whitney U para dois grupos (Rural vs Urbano)
#Para comparar IMC ou Idade entre Rural e Urbano:

# Mann-Whitney para IMC
wilcox.test(IMC ~ Region_type, data = metadados)

# Mann-Whitney para Idade
wilcox.test(Age ~ Region_type, data = metadados)

#Qui-quadrado ou Fisher para variáveis categóricas (Sexo)
#Se você quiser comparar Sexo entre regiões ou tipos de região:

#Para mais de dois grupos (Sul, Norte, etc.):


# Tabela cruzada para Sexo e Região
tabela_sexo <- table(metadados$Sex, metadados$Region)

# Qui-quadrado
chisq.test(tabela_sexo)


#Para dois grupos (Rural vs Urbano):
```{r}
# Tabela cruzada para Sexo e Tipo de Região
tabela_sexo <- table(metadados$Sex, metadados$Region_type)

# Teste exato de Fisher (para pequenas amostras)
fisher.test(tabela_sexo)

```
### IMPORTANTE PARA ESCOLHER OS TESTES ###
#Resumo de Como Escolher o Teste

#Variável: IMC, Idade	Grupos: 2 grupos (Rural/Urbano)	=> Mann-Whitney U

#Variável: IMC, Idade	Grupos: Mais de 2 grupos =>	Kruskal-Wallis

#Variável: Sexo (categórica)	Grupos: 2 grupos (Rural/Urbano)	=> Fisher ou Qui-quadrado

#Variável: Sexo (categórica)	Grupos: Mais de 2 grupos	=> Qui-quadrado


#Não deu diferença estatística!
#O que fazer a seguir?
####Verificar os dados###

#Certifique-se de que os dados de idade estão corretos e completos.
#Visualize os dados por região usando boxplots para identificar diferenças visuais ou outliers.
#Exemplo de código para boxplots:


library(ggplot2)

#Alterar a paleta de cores. Se quiser mudar as cores, substitua scale_fill_brewer(palette = "Set3") por outra paleta disponível, como: "Paired", "Pastel1", "Dark2" Ou use cores personalizadas: 

boxplot_idade <- ggplot(metadados, aes(x = Region, y = Age, fill = Region)) +
  geom_boxplot() +
  labs(title = "Distribuição da idade por região",
       x = "Região",
       y = "Idade") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")  # Paleta de cores

# Exibir o gráfico
boxplot_idade

# Salvar como PNG
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/boxplot_idade.png", plot = boxplot_idade, width = 8, height = 6, dpi = 300)

# Salvar como PDF
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/boxplot_idade.pdf", plot = boxplot_idade, width = 8, height = 6)



#Tabela com : Idade, Sexo, Média de Idade entre os sexos e Número Total. 
#Criar a tabela resumida
#Usaremos o dplyr para calcular os resumos e o gt para formatar a tabela.


# Carregar os pacotes necessários
library(dplyr)
library(gt)

# Ajustar os rótulos para inglês (Male e Female)
tabela_resumo <- metadados %>%
  mutate(Sex = ifelse(Sex == "Feminino", "Female", 
                      ifelse(Sex == "Masculino", "Male", Sex))) %>%
  group_by(Sex) %>%
  summarise(
    `N` = n(),  # Contagem total por sexo
    `Mean Age` = round(mean(Age, na.rm = TRUE), 2),  # Média de idade
    `SD` = round(sd(Age, na.rm = TRUE), 2),  # Desvio padrão da idade
    `Min Age` = min(Age, na.rm = TRUE),  # Idade mínima
    `Max Age` = max(Age, na.rm = TRUE)   # Idade máxima
  ) %>%
  ungroup()

# Criar uma tabela formatada com gt
tabela_final <- tabela_resumo %>%
  gt() %>%
  tab_header(
    title = "Age and Sex (n=130)"
  ) %>%
  fmt_number(
    columns = c(`Mean Age`, `SD`, `Min Age`, `Max Age`),
    decimals = 2
  ) %>%
  cols_label(
    Sex = "Sex",
    `N` = "Count",
    `Mean Age` = "Mean Age",
    `SD` = "Standard Deviation",
    `Min Age` = "Minimum Age",
    `Max Age` = "Maximum Age"
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.align = "center"
  )

# Exibir a tabela
tabela_final

# Salvar como HTML
gtsave(tabela_final, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/age_sex_summary.html")

# Salvar como PDF
gtsave(tabela_final, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/age_sex_summary.jpg")

gtsave(
  tabela_final, 
  "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/age_sex_summary.png"
)




### Parasitológico. ###
#Tabela que mostre quantos voluntários têm resultado positivo no parasitológico e listar os tipos de parasitas encontrados


# Carregar os pacotes necessários
library(dplyr)
library(gt)

# Substituir "e" por "and" na coluna Parasitas
metadados <- metadados %>%
  mutate(Parasitas = gsub(" e ", " and ", Parasitas))

# Substituir "Nana" por "nana" na coluna Parasitas
metadados <- metadados %>%
  mutate(Parasitas = gsub("Nana", "nana", Parasitas))

# Fazer as substituições na coluna Parasitas
metadados <- metadados %>%
  mutate(
    Parasitas = gsub("Prersença", "Presença", Parasitas)  # Corrigir "Prersença" para "Presença"
  )

# Substituir "/" por vírgula e ajustar o formato
metadados <- metadados %>%
  mutate(
    Parasitas = gsub(" . ", ", ", Parasitas)  # Substituir " / " por ", "
  )

# Verificar se as alterações foram realizadas corretamente
unique(metadados$Parasitas)

# Contar os resultados do parasitológico
resumo_parasitologico <- metadados %>%
  group_by(ParasitologicoPositivo) %>%
  summarise(N = n()) %>%
  arrange(desc(N))

# Exibir o resumo
print(resumo_parasitologico)


# Substituir "e" por "and" na coluna de tipos de parasitas
tabela_parasitologico <- metadados %>%
  filter(ParasitologicoPositivo == "Sim") %>%  # Selecionar apenas positivos no parasitológico
  group_by(Parasitas) %>%  # Agrupar por tipos de parasitas
  summarise(`Número de casos` = n()) %>%  # Contar ocorrências de cada parasita
  arrange(desc(`Número de casos`)) %>%  # Ordenar por número de casos
  ungroup()

# Criar uma tabela formatada com gt
tabela_final_parasitas <- tabela_parasitologico %>%
  gt() %>%
  tab_header(
    title = "Parasitological Test (n=130) ("
  ) %>%
  cols_label(
    Parasitas = "Parasite Type",
    `Número de casos` = "Number of Cases"
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.align = "center"
  )

# Exibir a tabela
tabela_final_parasitas


