---
  title: "R Notebook"
output: html_notebook
---
  



#####Carregar pacotes necessários
library(ggplot2)
library(tidyr)
library(dplyr)
install.packages("waffle")
library(waffle)
# Instalar o pacote plotly, se necessário
if (!require(plotly)) install.packages("plotly")
install.packages("plotly", dependencies = TRUE)

library(plotly)


# Transformar o objeto `venn_binary` para formato long
venn_binary_long <- venn_binary %>%
  rownames_to_column(var = "Volunteer_ID") %>%  # Adicionar IDs como coluna, se necessário
  pivot_longer(-Volunteer_ID, names_to = "Parameter", values_to = "Value") %>%
  mutate(Value = as.factor(Value)) # Converter para fator

# Carregar pacotes necessários
library(ggplot2)
library(tidyr)
library(dplyr)

# Transformar o objeto `venn_binary` para formato long
venn_binary_long <- venn_binary %>%
  rownames_to_column(var = "Volunteer_ID") %>%  # Adicionar IDs como coluna
  pivot_longer(-Volunteer_ID, names_to = "Parameter", values_to = "Value") %>%
  mutate(Value = factor(Value, levels = c("1", "0", NA), labels = c("Altered", "Normal", "Missing")))

## Carregar pacotes necessários
library(ggplot2)
library(tidyr)
library(dplyr)

# Transformar o objeto `venn_binary` para formato long
venn_binary_long <- venn_binary %>%
  rownames_to_column(var = "Volunteer_ID") %>%  # Adicionar IDs como coluna, se necessário
  pivot_longer(-Volunteer_ID, names_to = "Parameter", values_to = "Value") %>%
  mutate(Value = as.factor(Value)) # Converter para fator

# Criar o heatmap com ggplot2
ggplot(venn_binary_long, aes(x = Parameter, y = Volunteer_ID, fill = Value)) +
  geom_tile(color = "white" ) + # Criar os quadrados
  scale_fill_manual(
    values = c("1" = "blue", "0" = "red", "NA" = "gray"),
    labels = c("Altered", "Normal", "Missing"),
    na.value = "gray" # Definir a cor para NA
  ) +
  labs(
    title = "Visualização de Parâmetros por Voluntário",
    x = "Parâmetro",
    y = "Voluntário",
    fill = "Status"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6) # Ajusta o tamanho do texto para muitos voluntários
  )


```
```{r}
# Carregar pacotes necessários
library(ggplot2)
library(tidyr)
library(dplyr)

# Transformar o objeto `venn_binary` para formato long
venn_binary_long <- venn_binary %>%
  rownames_to_column(var = "Volunteer_ID") %>%  # Adicionar IDs como coluna
  mutate(Sum = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>% # Calcular soma dos valores numéricos
  arrange(Sum, Volunteer_ID) %>%                # Organizar voluntários pela soma e ID
  select(-Sum) %>%                              # Remover coluna de soma após organizar
  pivot_longer(-Volunteer_ID, names_to = "Parameter", values_to = "Value") %>%
  mutate(Value = as.factor(Value)) # Converter para fator

# Criar o heatmap com ggplot2
ggplot(venn_binary_long, aes(x = Parameter, y = Volunteer_ID, fill = Value)) +
  geom_tile(color = "white", height = 0.8, width = 1.2) + # Ajustar altura e largura dos retângulos
  scale_fill_manual(
    values = c("1" = "red", "0" = "blue", "NA" = "gray"),
    labels = c("Altered", "Normal", "Missing"),
    na.value = "gray" # Definir a cor para NA
  ) +
  labs(
    title = "Visualização de Parâmetros por Voluntário",
    x = "Parâmetro",
    y = "Voluntário",
    fill = "Status"
  ) +
  coord_fixed(ratio = 0.5) + # Ajustar proporção (retângulos mais altos e menos largos)
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(), # Remover rótulos do eixo Y para muitos voluntários
    axis.ticks.y = element_blank() # Remover ticks do eixo Y
  )

```
```{r}
# Carregar pacotes necessários
library(ggplot2)
library(tidyr)
library(dplyr)

# Transformar o objeto `venn_binary` para formato long
venn_binary_long <- venn_binary %>%
  rownames_to_column(var = "Volunteer_ID") %>%  # Adicionar IDs como coluna
  mutate(
    Num_Altered = rowSums(across(where(is.numeric)) == 1, na.rm = TRUE) # Contar número de valores alterados (1)
  ) %>%
  arrange(Num_Altered, Volunteer_ID) %>% # Ordenar por número de alterados e ID do voluntário
  pivot_longer(-c(Volunteer_ID, Num_Altered), names_to = "Parameter", values_to = "Value") %>%
  mutate(Value = as.factor(Value)) # Converter para fator

# Criar o heatmap com ggplot2
ggplot(venn_binary_long, aes(x = Parameter, y = Volunteer_ID, fill = Value)) +
  geom_tile(color = "white", height = 0.9, width = 4) + # Aumentar largura dos retângulos
  scale_fill_manual(
    values = c("1" = "red", "0" = "blue", "NA" = "gray"),
    labels = c("Altered", "Normal", "Missing"),
    na.value = "gray" # Definir a cor para NA
  ) +
  labs(
    title = "Visualização de Parâmetros por Voluntário",
    x = "Parâmetro",
    y = "Voluntário",
    fill = "Status"
  ) +
  coord_fixed(ratio = 0.3) + # Reduzir proporção para "alargar" o gráfico
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(), # Remover rótulos do eixo Y para muitos voluntários
    axis.ticks.y = element_blank(), # Remover ticks do eixo Y
    plot.title = element_text(hjust = 0.5) # Centralizar o título
  )

```

```{r}
# Transformar o objeto `venn_binary` para formato long
venn_binary_long <- venn_binary %>%
  rownames_to_column(var = "Volunteer_ID") %>%  # Adicionar IDs como coluna
  mutate(
    Proportion_Missing = rowMeans(is.na(across(-Volunteer_ID))) # Calcular proporção de missing
  ) %>%
  arrange(Proportion_Missing, Volunteer_ID) %>% # Ordenar: menos missing primeiro, mais missing depois
  pivot_longer(-c(Volunteer_ID, Proportion_Missing), names_to = "Parameter", values_to = "Value") %>%
  mutate(Value = as.factor(Value)) # Converter para fator

# Criar o heatmap com ggplot2
ggplot(venn_binary_long, aes(x = Parameter, y = Volunteer_ID, fill = Value)) +
  geom_tile(color = "white") + # Ajustar proporção das células
  scale_fill_manual(
    values = c("1" = "blue", "0" = "red", "NA" = "gray"),
    labels = c("Altered", "Normal", "Missing"),
    na.value = "gray" # Definir a cor para NA
  ) +
  labs(
    title = "Visualização de Parâmetros por Voluntário",
    x = "Parâmetro",
    y = "Voluntário",
    fill = "Status"
  ) +
  coord_fixed(ratio = 0.5) + # Proporção ajustada para retângulos mais largos
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), # Melhorar visibilidade do texto
    axis.text.y = element_blank(), # Remover texto do eixo Y
    axis.ticks.y = element_blank(), # Remover ticks do eixo Y
    plot.title = element_text(hjust = 0.5) # Centralizar o título
  )

```

```{r}
# Transformar o objeto `venn_binary` para formato long
venn_binary_long <- venn_binary %>%
  rownames_to_column(var = "Volunteer_ID") %>%  # Adicionar IDs como coluna
  mutate(
    Proportion_Missing = rowMeans(is.na(across(-Volunteer_ID))) # Calcular proporção de missing
  ) %>%
  arrange(Proportion_Missing, Volunteer_ID) %>% # Ordenar: menos missing primeiro, mais missing depois
  pivot_longer(-c(Volunteer_ID, Proportion_Missing), names_to = "Parameter", values_to = "Value") %>%
  mutate(Value = as.factor(Value)) # Converter para fator

# Criar o heatmap com ggplot2
ggplot(venn_binary_long, aes(x = Parameter, y = Volunteer_ID, fill = Value)) +
  geom_tile(color = "white", width = 0.9, height = 1.2) + # Ajustar proporção das células
  scale_fill_manual(
    values = c("1" = "blue", "0" = "red", "NA" = "gray"),
    labels = c("Altered", "Normal", "Missing"),
    na.value = "gray" # Definir a cor para NA
  ) +
  labs(
    title = "Visualização de Parâmetros por Voluntário",
    x = "Parâmetro",
    y = "Voluntário",
    fill = "Status"
  ) +
  coord_fixed(ratio = 0.5) + # Proporção ajustada para retângulos mais largos
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), # Melhorar visibilidade do texto
    axis.text.y = element_blank(), # Remover texto do eixo Y
    axis.ticks.y = element_blank(), # Remover ticks do eixo Y
    plot.title = element_text(hjust = 0.5) # Centralizar o título
  )


###=====================================###
      #Waffle Grafic"
###=====================================###

library(waffle)

# Verificando o número de valores não NA em cada coluna
na_counts <- sapply(metadados[, c("GLICOSE", "TRIGLICERIDES", "HDL", "Systolic", "Waist")], function(x) sum(!is.na(x)))

# Contar as linhas completas (sem NAs nas colunas específicas)
complete_rows <- nrow(metadados[complete.cases(metadados[, c("GLICOSE", "TRIGLICERIDES", "HDL", "Systolic", "Waist")]), ])

# Exibir o número de linhas completas
complete_rows

# Exibindo os resultados
na_counts


# Criar grupos com base na completude dos dados
metadados$Group <- ifelse(
  rowSums(is.na(metadados[, c("GLICOSE", "TRIGLICERIDES", "HDL", "Systolic", "Waist")])) > 0,
  "Incomplete Data",
  "Complete Data"
)

# Contar os grupos
group_counts <- table(metadados$Group)

# Criar o gráfico de waffle corrigido
waffle::waffle(
  group_counts,
  rows = 10,
  colors = c("red", "blue"),
  title = "Metabolic Syndrome Overview"
) +
  ggplot2::labs(fill = "Data Completeness") +  # Define o título da legenda
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))  # Centraliza o título


####==========================####
      # Grafico Sunburst
####==========================####


#############====================================================##############
#Defina os critérios da síndrome metabólica:
  
#  Glicose em jejum > 126
# Triglicérides > 150
#  HDL < 40 (homens) ou < 50 (mulheres)
#  Pressão sistólica ≥ 130 ou pressão diastólica ≥ 85
#  Circunferência da cintura > valores de corte (varia com sexo/população)

## Classifique os indivíduos:
  
# Determine quem tem dados completos.
# Classifique quem possui pelo menos 3 critérios elevados.
# Estruture os dados para o gráfico Sunburst:
  
# Crie uma hierarquia que mostra a distribuição de quem cumpre 0, 1, 2, 3 ou mais critérios.
#  e o gráfico usando plotly.
###############=======================================================##########



#=======================================================================================================================#
 # FOI MUITO DIFÍCIL FAZER GRAFICO AQUI. pASSEI A PLANILHA PRO EXCEL E FIZ 
# NO GOOGLE PLANILHAS E NO EXCEL! gRAFICOS EM :
# https://docs.google.com/spreadsheets/d/1IMo5ZZY_mgdf68pnI28XZMJdAAzpceaOTe1pJ2swaCg/edit?gid=370464566#gid=370464566

#salvei details_count em:

write.csv(x = details_count, 
          file = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/details_count_metabolic_parameters.csv", row.names = FALSE)
#=======================================================================================================================#

if (!require(plotly)) install.packages("plotly")


library(plotly)


# Definir os critérios de síndrome metabólica
metadados$High_Glucose <- metadados$GLICOSE > 126
metadados$High_Triglycerides <- metadados$TRIGLICERIDES > 150
metadados$Low_HDL <- ifelse(
  metadados$Sex == "Masculino", metadados$HDL < 40, metadados$HDL < 50
)
metadados$High_Blood_Pressure <- metadados$Systolic >= 130 | metadados$Diastolic >= 85
metadados$High_Waist <- metadados$Waist > ifelse(
  metadados$Sex == "Masculino", 102, 88
)

# Contar o número de critérios elevados por indivíduo
metadados$Criteria_Count <- rowSums(metadados[, c(
  "High_Glucose", "High_Triglycerides", "Low_HDL", "High_Blood_Pressure", "High_Waist"
)], na.rm = TRUE)

# Classificar os indivíduos
metadados$Group <- ifelse(
  metadados$Criteria_Count >= 3, "Metabolic Syndrome",
  ifelse(metadados$Criteria_Count > 0, "Partial Criteria", "Health Metabolic Profile")
)

# Preparar os dados para o Sunburst Chart
sunburst_data <- data.frame(
  labels = c("Total", "Metabolic Syndrome", "Partial Criteria", "Health Metabolic Profile"),
  parents = c("", "Total", "Total", "Total"),
  values = c(
    nrow(metadados),
    sum(metadados$Group == "Metabolic Syndrome"),
    sum(metadados$Group == "Partial Criteria"),
    sum(metadados$Group == "Health Metabolic Profile")
  )
)

#Filtrar "Partial Criteria": Crie um subconjunto do dataset que inclua apenas as pessoas com dados incompletos:
partial_criteria <- metadados[rowSums(is.na(metadados[, c("GLICOSE", "TRIGLICERIDES", "HDL", "Systolic", "Waist")])) > 0, ]

#Contar critérios elevados: Verifique quantas pessoas dentro de "Partial Criteria" têm cada critério elevado, considerando os valores de referência:

elevated_glucose <- sum(partial_criteria$GLICOSE > 126, na.rm = TRUE)
elevated_triglycerides <- sum(partial_criteria$TRIGLICERIDES > 150, na.rm = TRUE)
low_hdl <- sum(partial_criteria$HDL < 40, na.rm = TRUE)
high_systolic <- sum(partial_criteria$Systolic > 130, na.rm = TRUE)
high_waist_women <- sum(partial_criteria$Waist > 88 & partial_criteria$Sex == "Female", na.rm = TRUE)
high_waist_men <- sum(partial_criteria$Waist > 102 & partial_criteria$Sex == "Male", na.rm = TRUE)

high_waist <- high_waist_women + high_waist_men

#Incluir informações no gráfico: Combine os resultados em uma string para exibição no gráfico:
partial_text <- paste(
  "Partial Criteria:\n",
  "High Glucose: ", elevated_glucose, "\n",
  "High Triglycerides: ", elevated_triglycerides, "\n",
  "Low HDL: ", low_hdl, "\n",
  "High Systolic: ", high_systolic, "\n",
  "High Waist (Women > 88 cm): ", high_waist_women, "\n",
  "High Waist (Men > 102 cm): ", high_waist_men
)

#Atualizar o gráfico Sunburst: Inclua o texto no Sunburst ou como legenda adicional:

plot_ly(
  labels = c("Total", "Metabolic Syndrome", "No Criteria", "Partial Criteria", 
             "High Glucose", "High Triglycerides", "Low HDL", "High Systolic", "High Waist"),
  parents = c("", "Total", "Total", "Total", "Partial Criteria", "Partial Criteria", 
              "Partial Criteria", "Partial Criteria", "Partial Criteria"),
  values = c(130, metabolic_syndrome_count, no_criteria_count, nrow(partial_criteria),
             elevated_glucose, elevated_triglycerides, low_hdl, high_systolic, high_waist),
  type = "sunburst"
) %>%
  layout(title = "Health Metabolic Profile")


# Criar o gráfico Sunburst com números visíveis
fig <- plot_ly(
  data = sunburst_data,
  labels = ~labels,
  parents = ~parents,
  values = ~values,
  type = 'sunburst',
  branchvalues = 'total',
  textinfo = "label+value" # Adiciona os números no gráfico
)

# Mostrar o gráfico
fig






# Criar um novo objeto metadados2 como cópia de metadados
metadados2 <- metadados

# Criar uma coluna para contar critérios alterados
metadados2$Criteria_Altered <- rowSums(data.frame(
  Glucose = metadados2$GLICOSE > 126,
  Triglycerides = metadados2$TRIGLICERIDES > 150,
  HDL = ifelse(metadados2$Sex == "Feminino", metadados2$HDL < 50, metadados2$HDL < 40),
  BloodPressure = metadados2$Systolic > 130 | metadados2$Diastolic > 85,
  Waist = ifelse(metadados2$Sex == "Feminino", metadados2$Waist > 88, metadados2$Waist > 102)
), na.rm = TRUE)

# Identificar linhas com NA
metadados2$Complete_Data <- rowSums(is.na(metadados2[, c("GLICOSE", "TRIGLICERIDES", "HDL", "Systolic", "Diastolic", "Waist")])) == 0

# Categorizar indivíduos
metadados2$Group <- ifelse(
  metadados2$Complete_Data & metadados2$Criteria_Altered == 0, "Health Metabolic Profile",
  ifelse(
    metadados2$Complete_Data & metadados2$Criteria_Altered >= 3, "Metabolic Syndrome",
    ifelse(
      metadados2$Complete_Data & metadados2$Criteria_Altered > 0, "Partial Criteria (No NA)",
      ifelse(
        !metadados2$Complete_Data & metadados2$Criteria_Altered > 0, "Partial Criteria (With NA)",
        "Incomplete Data (No Alterations)"
      )
    )
  )
)

# Contar o número de indivíduos em cada grupo
group_counts <- table(metadados2$Group)
print(group_counts)

# Preparar os dados para o gráfico Sunburst
sunburst_data <- data.frame(
  labels = c("Total", 
             "Health Metabolic Profile", 
             "Metabolic Syndrome", 
             "Partial Criteria (No NA)", 
             "Partial Criteria (With NA)", 
             "Incomplete Data (No Alterations)"),
  parents = c("", "Total", "Total", "Total", "Total", "Total"),
  values = c(
    nrow(metadados2),  # Total
    sum(metadados2$Group == "Health Metabolic Profile"),
    sum(metadados2$Group == "Metabolic Syndrome"),
    sum(metadados2$Group == "Partial Criteria (No NA)"),
    sum(metadados2$Group == "Partial Criteria (With NA)"),
    sum(metadados2$Group == "Incomplete Data (No Alterations)")
  )
)

# Criar o gráfico Sunburst
library(plotly)

fig <- plot_ly(
  data = sunburst_data,
  labels = ~labels,
  parents = ~parents,
  values = ~values,
  type = 'sunburst',
  branchvalues = 'total'
)

# Adicionar título
fig <- fig %>%
  layout(
    title = list(
      text = "Health Metabolic Profile",
      font = list(size = 18),
      x = 0.5
    )
  )

# Mostrar o gráfico
fig



########################


# Criar um novo objeto metadados2 como cópia de metadados
metadados2 <- metadados

# Criar uma coluna para contar critérios alterados
metadados2$Criteria_Altered <- rowSums(data.frame(
  Glucose = metadados2$GLICOSE > 126,
  Triglycerides = metadados2$TRIGLICERIDES > 150,
  HDL = ifelse(metadados2$Sex == "Feminino", metadados2$HDL < 50, metadados2$HDL < 40),
  BloodPressure = metadados2$Systolic > 130 | metadados2$Diastolic > 85,
  Waist = ifelse(metadados2$Sex == "Feminino", metadados2$Waist > 88, metadados2$Waist > 102)
), na.rm = TRUE)

# Identificar linhas com NA
metadados2$Complete_Data <- rowSums(is.na(metadados2[, c("GLICOSE", "TRIGLICERIDES", "HDL", "Systolic", "Diastolic", "Waist")])) == 0

# Categorizar indivíduos
metadados2$Group <- ifelse(
  metadados2$Complete_Data & metadados2$Criteria_Altered == 0, "Health Metabolic Profile",
  ifelse(
    metadados2$Complete_Data & metadados2$Criteria_Altered >= 3, "Metabolic Syndrome",
    ifelse(
      metadados2$Complete_Data & metadados2$Criteria_Altered > 0, "<3 Altered Parameters",
      ifelse(
        !metadados2$Complete_Data & metadados2$Criteria_Altered > 0, "Incomplete Data with Alterations",
        "Incomplete Data: Health"
      )
    )
  )
)

# Contar alterações em cada grupo
partial_criteria_with_na <- metadados2[metadados2$Group == "Incomplete Data with Alterations", ]
less_than_3_altered <- metadados2[metadados2$Group == "<3 Altered Parameters", ]

# Criar texto detalhado para os grupos
partial_text <- paste0(
  "Incomplete Data with Alterations: ", nrow(partial_criteria_with_na), "\n",
  "  Glucose Altered: ", sum(partial_criteria_with_na$GLICOSE > 126, na.rm = TRUE), "\n",
  "  Triglycerides Altered: ", sum(partial_criteria_with_na$TRIGLICERIDES > 150, na.rm = TRUE), "\n",
  "  Low HDL: ", sum(ifelse(partial_criteria_with_na$Sex == "Feminino", partial_criteria_with_na$HDL < 50, partial_criteria_with_na$HDL < 40), na.rm = TRUE), "\n",
  "  High Blood Pressure: ", sum(partial_criteria_with_na$Systolic > 130 | partial_criteria_with_na$Diastolic > 85, na.rm = TRUE), "\n",
  "  High Waist: ", sum(ifelse(partial_criteria_with_na$Sex == "Feminino", partial_criteria_with_na$Waist > 88, partial_criteria_with_na$Waist > 102), na.rm = TRUE), "\n",
  "<3 Altered Parameters: ", nrow(less_than_3_altered), "\n",
  "  Glucose Altered: ", sum(less_than_3_altered$GLICOSE > 126, na.rm = TRUE), "\n",
  "  Triglycerides Altered: ", sum(less_than_3_altered$TRIGLICERIDES > 150, na.rm = TRUE), "\n",
  "  Low HDL: ", sum(ifelse(less_than_3_altered$Sex == "Feminino", less_than_3_altered$HDL < 50, less_than_3_altered$HDL < 40), na.rm = TRUE), "\n",
  "  High Blood Pressure: ", sum(less_than_3_altered$Systolic > 130 | less_than_3_altered$Diastolic > 85, na.rm = TRUE), "\n",
  "  High Waist: ", sum(ifelse(less_than_3_altered$Sex == "Feminino", less_than_3_altered$Waist > 88, less_than_3_altered$Waist > 102), na.rm = TRUE)
)

# Contar o número de indivíduos em cada grupo
group_counts <- table(metadados2$Group)
print(group_counts)

# Preparar os dados para o gráfico Sunburst
sunburst_data <- data.frame(
  labels = c("Total", 
             "Health Metabolic Profile", 
             "Metabolic Syndrome", 
             "<3 Altered Parameters", 
             "Incomplete Data with Alterations", 
             "Incomplete Data: Health"),
  parents = c("", "Total", "Total", "Total", "Total", "Total"),
  values = c(
    nrow(metadados2),  # Total
    sum(metadados2$Group == "Health Metabolic Profile"),
    sum(metadados2$Group == "Metabolic Syndrome"),
    sum(metadados2$Group == "<3 Altered Parameters"),
    sum(metadados2$Group == "Incomplete Data with Alterations"),
    sum(metadados2$Group == "Incomplete Data: Health")
  )
)

# Criar o gráfico Sunburst
library(plotly)

fig <- plot_ly(
  data = sunburst_data,
  labels = ~labels,
  parents = ~parents,
  values = ~values,
  type = 'sunburst',
  branchvalues = 'total'
)

# Adicionar título e texto explicativo
fig <- fig %>%
  layout(
    title = list(
      text = "Health Metabolic Profile",
      font = list(size = 18),
      x = 0.5
    ),
    annotations = list(
      list(
        x = 1.3,
        y = 0.5,
        text = partial_text,
        showarrow = FALSE,
        align = "left",
        font = list(size = 12)
      )
    )
  )

# Mostrar o gráfico
fig

# Fatores a serem analisados
criteria <- c("High_Glucose", "Low_HDL", "High_Triglycerides", "High_Blood_Pressure", "High_Waist")


# Criar uma nova coluna com a soma dos critérios alterados por linha
metadados2$AlteredCount <- rowSums(metadados2[, criteria], na.rm = TRUE)

# Visualizar os primeiros resultados
head(metadados2$AlteredCount)


# Filtrar dados completos (sem NAs nos critérios)
filtered_data <- metadados2[complete.cases(metadados2[, criteria]), ]

# Contar o número de pessoas por quantidade de fatores alterados
group_counts <- table(filtered_data$AlteredCount)

# Criar uma tabela com os resultados detalhados
results_table <- data.frame(
  FactorsAltered = as.integer(names(group_counts)),
  Count = as.numeric(group_counts)
)

# Adicionar os nomes das pessoas e quais critérios estão alterados em cada grupo
results_table$Details <- sapply(
  results_table$FactorsAltered,
  function(n) {
    # Identificar pessoas com exatamente 'n' fatores alterados
    subset <- filtered_data[filtered_data$AlteredCount == n, ]
    if (nrow(subset) == 0) return("None")
    apply(subset, 1, function(row) {
      altered_criteria <- criteria[row[criteria] == 1]
      paste0("ID: ", row["Sample.id"], " (", paste(altered_criteria, collapse = ", "), ")")
    }) %>%
      paste(collapse = "; ")
  }
)

# Visualizar a tabela
print(results_table)

# Opcional: Salvar como CSV
write.csv(results_table, "factors_altered_summary.csv", row.names = FALSE)


colnames(metadados2)



#####======================================#####

# Carregar a biblioteca necessária
library(plotly)

# Preparar os dados para o gráfico Sunburst
sunburst_data <- data.frame(
  labels = c(
    "Total", 
    "High Glucose", 
    "High Triglycerides", 
    "Low HDL", 
    "High Blood Pressure", 
    "High Waist", 
    "Combinations"
  ),
  parents = c(
    "",              # "Total" não tem pai
    "Total",         # Critérios principais têm "Total" como pai
    "Total",
    "Total",
    "Total",
    "Total",
    "Total"          # Combinações têm "Total" como pai
  ),
  values = c(
    130,                                        # Total de linhas
    sum(grepl("High Glucose", metadados$Partial_Details, na.rm = TRUE)),  # Contagem de "High Glucose"
    sum(grepl("High Triglycerides", metadados$Partial_Details, na.rm = TRUE)),  # Contagem de "High Triglycerides"
    sum(grepl("Low HDL", metadados$Partial_Details, na.rm = TRUE)),             # Contagem de "Low HDL"
    sum(grepl("High Blood Pressure", metadados$Partial_Details, na.rm = TRUE)), # Contagem de "High Blood Pressure"
    sum(grepl("High Waist", metadados$Partial_Details, na.rm = TRUE)),          # Contagem de "High Waist"
    sum(!is.na(metadados$Partial_Details) & grepl("&", metadados$Partial_Details)) # Combinações
  )
)

# Criar o gráfico Sunburst
fig <- plot_ly(
  data = sunburst_data,
  labels = ~labels,
  parents = ~parents,
  values = ~values,
  type = 'sunburst',
  branchvalues = 'total'
)

# Adicionar título e ajustar layout
fig <- fig %>%
  layout(
    title = list(
      text = "Health Metabolic Sunburst Overview",
      font = list(size = 18),
      x = 0.5  # Centraliza o título horizontalmente
    )
  )

# Mostrar o gráfico
fig


##========================##

# Carregar a biblioteca necessária
library(plotly)

# Contar os detalhes únicos e suas ocorrências
details_count <- as.data.frame(table(metadados$Partial_Details, useNA = "ifany"))
colnames(details_count) <- c("Details", "Count")

# Renomear NA como "No Alterations"
details_count$Details[is.na(details_count$Details)] <- "No Alterations"

# Criar o gráfico de pizza
fig <- plot_ly(
  data = details_count,
  labels = ~Details,
  values = ~Count,
  type = 'pie',
  textinfo = 'label+percent',
  hoverinfo = 'label+value+percent',
  marker = list(colors = colorRampPalette(c("blue", "orange", "green", "red", "purple"))(nrow(details_count)))
)

# Adicionar título ao gráfico
fig <- fig %>%
  layout(
    title = list(
      text = "Health Metabolic Profile - Details",
      font = list(size = 18),
      x = 0.5  # Centraliza o título horizontalmente
    ),
    showlegend = TRUE
  )

# Mostrar o gráfico
fig

#======#

# Substituir valores NA em "Details" por "Not measured"
# Isso garante que categorias sem medição apareçam no gráfico como "Not measured"
details_count$Details[is.na(details_count$Details)] <- "Not measured"

# Ordenar a tabela details_count em ordem decrescente de "Count"
# Isso ajuda a organizar as categorias no gráfico, colocando as mais relevantes no topo
details_count <- details_count[order(-details_count$Count), ]

# Criar o gráfico de pizza usando o pacote plotly
# Aqui utilizamos a tabela "details_count", que contém as categorias e suas contagens
fig <- plot_ly(
  data = details_count,                  # Dados do gráfico
  labels = ~Details,                    # Rótulos das categorias (coluna "Details")
  values = ~Count,                      # Valores das fatias (coluna "Count")
  type = 'pie',                         # Tipo do gráfico: pizza
  textinfo = 'label+percent',           # Mostrar rótulo e porcentagem no gráfico
  hoverinfo = 'label+value+percent',    # Mostrar detalhes ao passar o mouse: rótulo, valor e porcentagem
  marker = list(
    # Paleta de cores dinâmica com 7 cores básicas ajustadas ao número de categorias
    colors = colorRampPalette(c("blue", "orange", "green", "red", "purple", "brown", "pink"))(nrow(details_count))
  )
)

# Personalizar o layout do gráfico
fig <- fig %>%
  layout(
    # Adicionar título ao gráfico
    title = list(
      text = "Health Metabolic Profile - Detailed Breakdown",  # Texto do título
      font = list(size = 18),                                  # Tamanho da fonte do título
      x = 0.5                                                  # Centralizar o título horizontalmente
    ),
    # Personalizar a legenda
    legend = list(
      title = list(text = "Criteria Details"),  # Título da legenda
      font = list(size = 12),                   # Tamanho da fonte da legenda
      x = 1.1,                                  # Posição horizontal da legenda (à direita)
      y = 0.5                                   # Posição vertical da legenda (no centro)
    )
  )

# Mostrar o gráfico
fig


