####Metabolic Syndrome References#####

  #Referências para valores glicemia, insulina e HOMA-IR: American Diabetes Association para glicemia. Diretriz da Sociedade Brasileira de Diabetes – Edição 2024 Aprovado pelo Comitê Central – DOI: 10.29327/5412848 / ISBN: 978-65-272-0704-7


#Cálculo de TyG: O índice TyG (Triglyceride-Glucose Index) é uma fórmula amplamente utilizada como marcador para avaliar a resistência à insulina, especialmente em estudos populacionais e na prática clínica. É melhor do que o HOMA-IR, pq mede glicemia de jejum e também triglicerides



#=========================#
       #Library
#=========================#
library(VennDiagram)
library(ggvenn)
library(ComplexUpset)
library(ggplot2)
install.packages("UpSetR")
library(UpSetR)
library(ComplexUpset)
library(ggplot2)
library(dplyr)


colnames(metadados)

####=================####
        #TyG
####=================####

# Cálculo do índice TyG
metadados$TyG <- log(metadados$TRIGLICERIDES * metadados$GLICOSE) / 2

# Visualizar os primeiros valores calculados
head(metadados$TyG)

# Resumo estatístico da nova coluna TyG
summary(metadados$TyG)

# Exportar o dataframe atualizado, se necessário
write.csv(metadados, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados.csv", row.names = FALSE)


#Criar um objeto só com dados relevantes para Sindrome Metabolica!


# Criar um novo dataframe com as colunas selecionadas
metadados_MetS <- metadados[, c(
  "Sample.id", "Region_type", "Region", "Age", "Sex", "Raca", "Menopausa", "Gravida",
  "Fuma", "Alcool", "Metformina", "Estatinas", "Systolic", "Diastolic", "Weight",
  "Height", "Waist", "Hip", "HbA1c", "Risco_diabetes", "COLESTEROL", "LDL", "HDL",
  "VLDL", "TRIGLICERIDES", "GLICOSE", "INSULINA", "HOMA", "PCR", "IMC", "TyG"
)
]

# Verificar o novo dataframe
head(metadados_MetS)

# Salvar como CSV
write.csv(metadados_MetS, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados_MetS.csv", row.names = FALSE)



# Calcular média e desvio padrão
mean_TyG <- mean(metadados$TyG, na.rm = TRUE)  # Média (ignorar NA)
sd_TyG <- sd(metadados$TyG, na.rm = TRUE)      # Desvio padrão (ignorar NA)

# Exibir os resultados

cat("Média do TyG:", mean_TyG, "\nDesvio Padrão do TyG:", sd_TyG)


#acrescentar coluna Waist/Hip Ratio

metadados$WHR <- metadados$Waist / metadados$Hip


# Exemplo de listas com metabólitos ou associações para cada variável
dados_metabolicos <- list(
  GLICOSE = metadados$GLICOSE[!is.na(metadados$GLICOSE)], # Filtra valores não NA
  INSULINA = metadados$INSULINA[!is.na(metadados$INSULINA)],
  HOMA = metadados$HOMA[!is.na(metadados$HOMA)],
  TyG = metadados$HOMA[!is.na(metadados$HOMA)],
  IMC = metadados$IMC[!is.na(metadados$IMC)],
  WHR = metadados$WHR[!is.na(metadados$WHR)],
  LDL = metadados$LDL[!is.na(metadados$LDL)],
  HDL = metadados$HDL[!is.na(metadados$HDL)],
  VLDL = metadados$VLDL[!is.na(metadados$VLDL)],
  TRIGLICERIDES = metadados$TRIGLICERIDES[!is.na(metadados$TRIGLICERIDES)],
  TyG = metadados$TyG [!is.na(metadados$TyG)],
  COLESTEROL = metadados$COLESTEROL[!is.na(metadados$COLESTEROL)]
)


# Filtrar os dados baseados nos critérios
dados_venn <- list(
  GLICOSE = metadados$Sample.id[metadados$GLICOSE > 126 & !is.na(metadados$GLICOSE)],
  HbA1c = metadados$Sample.id[metadados$HbA1c > 6.5 & !is.na(metadados$HbA1c)],
  INSULINA = metadados$Sample.id[metadados$INSULINA > 29.0 & !is.na(metadados$INSULINA)]
)

# Criar o diagrama de Venn
library(VennDiagram)



# Criar o diagrama de Venn
venn.plot <- draw.triple.venn(
  area1 = sum(metadados$GLICOSE > 126, na.rm = TRUE),
  area2 = sum(metadados$HbA1c > 6.5, na.rm = TRUE),
  area3 = sum(metadados$INSULINA > 29, na.rm = TRUE),
  n12 = sum(metadados$GLICOSE > 126 & metadados$HbA1c > 6.5, na.rm = TRUE),
  n13 = sum(metadados$GLICOSE > 126 & metadados$INSULINA > 29, na.rm = TRUE),
  n23 = sum(metadados$HbA1c > 6.5 & metadados$INSULINA > 29, na.rm = TRUE),
  n123 = sum(metadados$GLICOSE > 126 & metadados$HbA1c > 6.5 & metadados$INSULINA > 29, na.rm = TRUE),
  category = c("GLUCOSE > 126", "HbA1c > 6.5", "INSULIN > 29"),
  fill = c("blue", "yellow", "green"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
)


# Adicionar título e número total de indivíduos
grid.draw(venn.plot) # Desenha o Venn
grid.text("GLYCEMIA, HbA1c, and INSULIN LEVELS", x = 0.5, y =0.99, gp = gpar(fontsize = 12, fontface = "bold"))
grid.text(paste0("N= ", nrow(metadados[!is.na(metadados$GLICOSE) & 
                                         !is.na(metadados$HbA1c) & 
                                         !is.na(metadados$INSULINA), ])), 
          x = 0.5, y = 0.96, gp = gpar(fontsize = 12))

library(VennDiagram)

# Definir os critérios
HbA1c_high <- metadados$HbA1c > 6.5
Insulin_high <- metadados$INSULINA > 29
Glucose_high <- metadados$GLICOSE > 126
Health <- !(HbA1c_high | Insulin_high | Glucose_high) # Sem nenhum critério elevado

# Criar o diagrama de Venn
venn.plot2 <- draw.quad.venn(
  area1 = sum(HbA1c_high, na.rm = TRUE),           # HbA1c > 6.5
  area2 = sum(Insulin_high, na.rm = TRUE),         # Insulina > 29
  area3 = sum(Glucose_high, na.rm = TRUE),         # Glicose > 126
  area4 = sum(Health, na.rm = TRUE),               # Normais
  n12 = sum(HbA1c_high & Insulin_high, na.rm = TRUE),
  n13 = sum(HbA1c_high & Glucose_high, na.rm = TRUE),
  n14 = sum(HbA1c_high & Health, na.rm = TRUE),
  n23 = sum(Insulin_high & Glucose_high, na.rm = TRUE),
  n24 = sum(Insulin_high & Health, na.rm = TRUE),
  n34 = sum(Glucose_high & Health, na.rm = TRUE),
  n123 = sum(HbA1c_high & Insulin_high & Glucose_high, na.rm = TRUE),
  n124 = sum(HbA1c_high & Insulin_high & Health, na.rm = TRUE),
  n134 = sum(HbA1c_high & Glucose_high & Health, na.rm = TRUE),
  n234 = sum(Insulin_high & Glucose_high & Health, na.rm = TRUE),
  n1234 = sum(HbA1c_high & Insulin_high & Glucose_high & Health, na.rm = TRUE),
  category = c("HbA1c > 6.5", "Insulin > 29", "Glucose > 126", "Health"),
  fill = c("red", "yellow", "blue", "gray"),       # Adiciona a cor cinza para "Health"
  alpha = 0.5,                                     # Transparência
  cex = 1.5,                                       # Tamanho do texto
  cat.cex = 1.5                                    # Tamanho das categorias
)

# Adicionar título
grid.text("HbA1c, INSULIN and  GLUCOLE LEVELS", x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))



# Venn diagram interessante, mas pensar em outros dados para mostrar, ou outros parâmetros, talvez.

# Criar uma coluna para identificar "Health" (indivíduos sem nenhuma condição elevada)
metadados$Health <- !(metadados$HbA1c > 6.5 | metadados$INSULINA > 29 | metadados$GLICOSE > 126)

# Criar uma lista para os dados do Venn
venn_data <- list(
  "HbA1c > 6.5" = which(metadados$HbA1c > 6.5),
  "Glucose > 126" = which(metadados$GLICOSE > 126),
  "Insulin > 29" = which(metadados$INSULINA > 29),
  "Health" = which(metadados$Health) # Aqui usamos a coluna lógica corrigida
)

# Criar o gráfico com ggvenn
ggvenn(
  venn_data,
  fill_color = c("red", "blue", "yellow", "gray"),
  show_percentage = TRUE,
  stroke_size = 0.5,
  set_name_size = 6,
  text_size = 4
) +
  ggtitle("HbA1c, Glucose, Insulin Levels") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))


###??? Será que faz sentido? Colocar HOMA com glicemia e insulina ???###


# Instalar e carregar o pacote
install.packages("VennDiagram")
library(VennDiagram)

# Criar o diagrama de Venn
venn.plot3 <- draw.triple.venn(
  area1 = 30,   # Tamanho do conjunto 1
  area2 = 40,   # Tamanho do conjunto 2
  area3 = 50,   # Tamanho do conjunto 3
  n12 = 15,     # Sobreposição entre 1 e 2
  n13 = 5,      # Sobreposição entre 1 e 3
  n23 = 10,     # Sobreposição entre 2 e 3
  n123 = 2,     # Sobreposição entre todos os 3
  category = c("Fasting Glucose", "Insulin", "HOMA-IR"),
  fill = c("blue", "yellow", "green"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5
)

# Exibir o diagrama na visualização do R
grid.draw(venn.plot3)



##Só teste!!! ###

# Instalar e carregar o pacote
install.packages("ggvenn")
library(ggvenn)

# Criar os dados de exemplo
dados <- list(
  GLICOSE = c("A", "B", "C", "D"),
  INSULINA = c("B", "C", "E", "F"),
  HOMA = c("C", "D", "G", "H")
)

# Gerar o diagrama de Venn e exibir
ggvenn(dados, fill_color = c("blue", "yellow", "green"), stroke_size = 0.5, text_size = 4)


# Criar colunas para cada critério baseado nas condições
metadados.alpha.all <- metadados.alpha.all %>%
  mutate(
    # Critério 1: Obesidade abdominal (agora usando cm diretamente)
    Obesidade_Abdominal = ifelse(
      (Sex == "Masculino" & Waist > 101.6) | 
        (Sex == "Feminino" & Waist > 88.9),
      TRUE, FALSE
    ),
    
    # Critério 2: Triglicerídeos elevados
    Triglicerideos_Altos = ifelse(TRIGLICERIDES >= 150, TRUE, FALSE),
    
    # Critério 3: HDL baixo
    HDL_Baixo = ifelse(
      (Sex == "Masculino" & HDL < 40) | 
        (Sex == "Feminino" & HDL < 50),
      TRUE, FALSE
    ),
    
    # Critério 4: Pressão arterial elevada
    Pressao_Alta = ifelse(
      Systolic >= 130 | Diastolic >= 85,
      TRUE, FALSE
    ),
    
    # Critério 5: Glicemia de jejum elevada
    Glicemia_Alta = ifelse(GLICOSE >= 100, TRUE, FALSE)
  )



#Incluindo saudáveis

# Adicionar a coluna 'Saudavel' ao dataframe
metadados.alpha.all <- metadados.alpha.all %>%
  mutate(
    Healthy = ifelse(
      !Obesidade_Abdominal & 
        !Triglicerideos_Altos & 
        !HDL_Baixo & 
        !Pressao_Alta & 
        !Glicemia_Alta,
      TRUE, FALSE
    )
  )

colnames(metadados.alpha.all)

# Criar listas baseadas nos critérios
venn_data2 <- list(
  "Abdominal Obesity" = metadados.alpha.all$Sample.id[metadados.alpha.all$Obesidade_Abdominal],
  "High Triglycerides" = metadados.alpha.all$Sample.id[metadados.alpha.all$Triglicerideos_Altos],
  "Low HDL" = metadados.alpha.all$Sample.id[metadados.alpha.all$HDL_Baixo],
  "High Blood Pressure" = metadados.alpha.all$Sample.id[metadados.alpha.all$Pressao_Alta],
  "High Fasting Glucose" = metadados.alpha.all$Sample.id[metadados.alpha.all$Glicemia_Alta]
)

# Instalar e carregar ggvenn
install.packages("ggvenn")
library(ggvenn)

# Criar o diagrama de Venn
ggvenn(venn_data2, 
       fill_color = c("blue", "yellow", "green", "red", "purple"), # Cores para cada critério
       stroke_size = 0.5, 
       text_size = 4,
       show_percentage = TRUE) # Mostrar porcentagem




# 1. Remover linhas com valores NA
venn_binary_clean <- na.omit(venn_binary)  # Remover voluntários com NA

# 2. Criar o diagrama de Venn
venn_plot <- draw.quintuple.venn(
  area1 = sum(venn_binary_clean$Abdominal_Obesity),
  area2 = sum(venn_binary_clean$High_Triglycerides),
  area3 = sum(venn_binary_clean$Low_HDL),
  area4 = sum(venn_binary_clean$High_Blood_Pressure),
  area5 = sum(venn_binary_clean$High_Fasting_Glucose),
  
  n12 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Triglycerides),
  n13 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$Low_HDL),
  n14 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Blood_Pressure),
  n15 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Fasting_Glucose),
  n23 = sum(venn_binary_clean$High_Triglycerides & venn_binary_clean$Low_HDL),
  n24 = sum(venn_binary_clean$High_Triglycerides & venn_binary_clean$High_Blood_Pressure),
  n25 = sum(venn_binary_clean$High_Triglycerides & venn_binary_clean$High_Fasting_Glucose),
  n34 = sum(venn_binary_clean$Low_HDL & venn_binary_clean$High_Blood_Pressure),
  n35 = sum(venn_binary_clean$Low_HDL & venn_binary_clean$High_Fasting_Glucose),
  n45 = sum(venn_binary_clean$High_Blood_Pressure & venn_binary_clean$High_Fasting_Glucose),
  n123 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Triglycerides & venn_binary_clean$Low_HDL),
  n124 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Triglycerides & venn_binary_clean$High_Blood_Pressure),
  n125 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Triglycerides & venn_binary_clean$High_Fasting_Glucose),
  n134 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$Low_HDL & venn_binary_clean$High_Blood_Pressure),
  n135 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$Low_HDL & venn_binary_clean$High_Fasting_Glucose),
  n145 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Blood_Pressure & venn_binary_clean$High_Fasting_Glucose),
  n234 = sum(venn_binary_clean$High_Triglycerides & venn_binary_clean$Low_HDL & venn_binary_clean$High_Blood_Pressure),
  n235 = sum(venn_binary_clean$High_Triglycerides & venn_binary_clean$Low_HDL & venn_binary_clean$High_Fasting_Glucose),
  n245 = sum(venn_binary_clean$High_Triglycerides & venn_binary_clean$High_Blood_Pressure & venn_binary_clean$High_Fasting_Glucose),
  n345 = sum(venn_binary_clean$Low_HDL & venn_binary_clean$High_Blood_Pressure & venn_binary_clean$High_Fasting_Glucose),
  n1234 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Triglycerides & venn_binary_clean$Low_HDL & venn_binary_clean$High_Blood_Pressure),
  n1235 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Triglycerides & venn_binary_clean$Low_HDL & venn_binary_clean$High_Fasting_Glucose),
  n1245 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Triglycerides & venn_binary_clean$High_Blood_Pressure & venn_binary_clean$High_Fasting_Glucose),
  n1345 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$Low_HDL & venn_binary_clean$High_Blood_Pressure & venn_binary_clean$High_Fasting_Glucose),
  n2345 = sum(venn_binary_clean$High_Triglycerides & venn_binary_clean$Low_HDL & venn_binary_clean$High_Blood_Pressure & venn_binary_clean$High_Fasting_Glucose),
  n12345 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Triglycerides & venn_binary_clean$Low_HDL & venn_binary_clean$High_Blood_Pressure & venn_binary_clean$High_Fasting_Glucose),
  
  category = c("Abdominal Obesity", "Triglycerides.", "HDL", 
               "Blood Pressure", "Glucose"),
  fill = c("blue", "yellow", "green", "red", "purple"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5)
)

# 3. Visualizar o diagrama de Venn no RStudio
grid.newpage()
grid.draw(venn_plot)

# 4. Salvar em PDF
pdf("venn_diagram_metabolic_syndrome.pdf", width = 8, height = 6)
grid.draw(venn_plot)
dev.off()

# 5. Salvar em JPG
jpeg("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/venn_diagram_metabolic_syndrome.jpg", width = 800, height = 600, units = "px", res = 300)
grid.draw(venn_plot)
dev.off()


colnames(venn_binary)

# Instalar e carregar UpSetR
install.packages("UpSetR")
library(UpSetR)

# Converter os dados de Venn para uma estrutura binária
venn_binary <- data.frame(
  Central_Obesity = metadados.alpha.all$Obesidade_Abdominal,
  High_Triglycerides = metadados.alpha.all$Triglicerideos_Altos,
  Low_HDL = metadados.alpha.all$HDL_Baixo,
  High_Blood_Pressure = metadados.alpha.all$Pressao_Alta,
  High_Fasting_Glucose = metadados.alpha.all$Glicemia_Alta,
  Healthy= metadados.alpha.all$Healthy
)

# Adicionar a coluna 'Healthy' e lidar com os valores NA
venn_binary <- venn_binary %>%
  mutate(
    Central_Obesity = ifelse(is.na(Central_Obesity), "Not Measured", Central_Obesity),
    High_Triglycerides = ifelse(is.na(High_Triglycerides), "Not Measured", High_Triglycerides),
    Low_HDL = ifelse(is.na(Low_HDL), "Not Measured", Low_HDL),
    High_Blood_Pressure = ifelse(is.na(High_Blood_Pressure), "Not Measured", High_Blood_Pressure),
    High_Fasting_Glucose = ifelse(is.na(High_Fasting_Glucose), "Not Measured", High_Fasting_Glucose),
    
    # Adicionando a categoria 'Healthy'
    Healthy = ifelse(
      Central_Obesity == FALSE & 
        High_Triglycerides == FALSE & 
        Low_HDL == FALSE & 
        High_Blood_Pressure == FALSE & 
        High_Fasting_Glucose == FALSE, 
      TRUE, 
      FALSE
    )
  )

# Instalar e carregar os pacotes necessários
if (!requireNamespace("ComplexUpset", quietly = TRUE)) install.packages("ComplexUpset")
library(ComplexUpset)
library(ggplot2)
library(dplyr)

# Remover linhas com "Not Measured" em qualquer coluna
venn_binary_clean <- venn_binary %>%
  filter(across(everything(), ~ . != "Not Measured"))

# Verificar o dataset após a remoção
print(venn_binary_clean)

# Criar o gráfico UpSet com o dataset limpo
ComplexUpset::upset(
  venn_binary_clean,
  intersect = c("Central_Obesity", "High_Triglycerides", "Low_HDL", "High_Blood_Pressure", "High_Fasting_Glucose"),
  mode = "exclusive_intersection", # Exclusivo para interseções
  base_annotations = list(
    'Intersection Size' = ComplexUpset::intersection_size(
      counts = TRUE,
      text = list(size = 3)
    )
  ),
  set_sizes = ComplexUpset::upset_set_size()
)

# Criar o gráfico UpSet
upset(venn_binary_clean, 
      sets = c("Central_Obesity", "High_Triglycerides", "Low_HDL", "High_Blood_Pressure", "High_Fasting_Glucose", ),
      order.by = "freq", # Ordenar pelas frequências
      main.bar.color = "darkblue", 
      sets.bar.color = "darkred",
      text.scale = 1.5)


library(UpSetR)

# Verificar as primeiras linhas do dataset
head(venn_binary_clean)

# Transformar as colunas em formato binário (0/1) ou lógico (TRUE/FALSE)
venn_binary_clean <- venn_binary %>%
  mutate(
    Central_Obesity = ifelse(Central_Obesity == "Present", 1, 0),
    High_Triglycerides = ifelse(High_Triglycerides == "Present", 1, 0),
    Low_HDL = ifelse(Low_HDL == "Present", 1, 0),
    High_Blood_Pressure = ifelse(High_Blood_Pressure == "Present", 1, 0),
    High_Fasting_Glucose = ifelse(High_Fasting_Glucose == "Present", 1, 0)
  )


# Transformar os valores de 'Present', 'Absent', e 'Not Measured'
venn_binary <- venn_binary %>%
  mutate(
    Central_Obesity = case_when(
      Central_Obesity == "Present" ~ 1,
      Central_Obesity == "Absent" ~ 0,
      Central_Obesity == "Not Measured" ~ NA_real_
    ),
    High_Triglycerides = case_when(
      High_Triglycerides == "Present" ~ 1,
      High_Triglycerides == "Absent" ~ 0,
      High_Triglycerides == "Not Measured" ~ NA_real_
    ),
    Low_HDL = case_when(
      Low_HDL == "Present" ~ 1,
      Low_HDL == "Absent" ~ 0,
      Low_HDL == "Not Measured" ~ NA_real_
    ),
    High_Blood_Pressure = case_when(
      High_Blood_Pressure == "Present" ~ 1,
      High_Blood_Pressure == "Absent" ~ 0,
      High_Blood_Pressure == "Not Measured" ~ NA_real_
    ),
    High_Fasting_Glucose = case_when(
      High_Fasting_Glucose == "Present" ~ 1,
      High_Fasting_Glucose == "Absent" ~ 0,
      High_Fasting_Glucose == "Not Measured" ~ NA_real_
    )
  )

# Exibir a tabela transformada
print(venn_binary)

# Criar o objeto venn_binary_clean removendo linhas com NA
venn_binary_clean <- venn_binary %>%
  drop_na()

# Verificar o objeto resultante
print(venn_binary_clean)

UpSetR::upset(
  venn_binary_clean, 
  sets = c("Central_Obesity", "High_Triglycerides", "Low_HDL", "High_Blood_Pressure", "High_Fasting_Glucose"),
  order.by = "freq", # Ordenar pelas frequências
  main.bar.color = "darkblue", 
  sets.bar.color = "darkred",
  text.scale = 1.5
)


hist(metadados$TRIGLICERIDES)

hist(metadados$GLICOSE)

# Verificar se a substituição foi realizada
table(metadados$GLICOSE)  # Frequência dos valores para conferir

hist(dados$HDL)

hist(metadados$Waist)

str(metadados)

hist(metadados$Systolic)

hist(metadados$Diastolic)

#corrigir diastolic! Está com valores errados!

# Substituir valores 9 por 90 na coluna Diastolic
metadados$Diastolic <- ifelse(metadados$Diastolic == 9, 90, metadados$Diastolic)

# Verificar se a substituição foi realizada
table(metadados$Diastolic)  # Frequência dos valores para conferir

hist(metadados$TyG)

hist(metadados$HOMA)

hist(metadados$HbA1c)

hist(metadados$PCR)

hist(metadados$LDL)


venn_binary_clean <- venn_binary_clean %>%
  mutate(All_Absent = ifelse(
    Central_Obesity == 0 & 
      High_Triglycerides == 0 & 
      Low_HDL == 0 & 
      High_Blood_Pressure == 0 & 
      High_Fasting_Glucose == 0, 
    1, 
    0
  ))


# Gerar o gráfico UpSet
UpSetR::upset(
  venn_binary_clean, 
  sets = c("Central_Obesity", "High_Triglycerides", "Low_HDL", "High_Blood_Pressure", "High_Fasting_Glucose"),
  queries = list(
    list(query = intersects, 
         params = list("Central_Obesity", "High_Triglycerides", "Low_HDL", "High_Blood_Pressure", "High_Fasting_Glucose"),
         color = "gray", 
         active = TRUE, 
         query.name = "All Parameters Absent")
  ),
  order.by = "freq", # Ordenar pelas frequências
  keep.order = TRUE, # Manter a ordem dos conjuntos
  main.bar.color = "darkblue", 
  sets.bar.color = "darkred",
  text.scale = 1.5
)


# Gerar o gráfico UpSet com as barras horizontais ordenadas
#grafico com 74 voluntarios que tem todos os fatores medidos


UpSetR::upset(
  venn_binary_clean, 
  sets = c("Central_Obesity", "High_Triglycerides", "Low_HDL", "High_Blood_Pressure", "High_Fasting_Glucose"),
  queries = list(
    list(query = intersects, 
         params = list("Central_Obesity", "High_Triglycerides", "Low_HDL", "High_Blood_Pressure", "High_Fasting_Glucose"),
         color = "lightgray", 
         active = TRUE, 
         query.name = "All Parameters Absent")
  ),
  order.by = "freq", # Ordenar pelas frequências
  main.bar.color = "darkblue", 
  sets.bar.color = "darkred",
  text.scale = 1.5
)



# Carregar biblioteca
library(ComplexUpset)
library(ggplot2)

# Criar uma coluna para identificar o grupo com todos os parâmetros ausentes
venn_binary_clean$All_Absent <- rowSums(venn_binary_clean[, c("Central_Obesity", "High_Triglycerides", 
                                                              "Low_HDL", "High_Blood_Pressure", 
                                                              "High_Fasting_Glucose")]) == 0

# Preparar o gráfico usando ComplexUpset
ComplexUpset::upset(
  venn_binary_clean,
  intersect = c("Central_Obesity", "High_Triglycerides", "Low_HDL", "High_Blood_Pressure", "High_Fasting_Glucose"),
  annotations = list(
    'Intersection Size' = ComplexUpset::intersection_size(counts = TRUE, text = list(size = 3)),
    'All Parameters Absent' = (
      ggplot(mapping = aes(x = intersection, fill = All_Absent)) +
        geom_bar(stat = 'count') +
        scale_fill_manual(values = c("TRUE" = "gray", "FALSE" = "darkblue")) +
        labs(y = "Count", fill = "All Absent")
    )
  ),
  mode = "exclusive_intersection", # Apenas interseções exclusivas
  base_annotations = list(
    'Intersection Size' = ComplexUpset::intersection_size(counts = TRUE, text = list(size = 3))
  ),
  queries = list(
    upset_query(set = "All_Absent", color = "gray", fill = "gray", only_components = TRUE)
  ),
  set_sizes = ComplexUpset::upset_set_size(),
  keep_empty_groups = FALSE # Remover grupos vazios
)

# Criar o diagrama de Venn com voluntários sem NAs


venn_plot <- draw.quintuple.venn(
  area1 = sum(venn_binary_clean$Abdominal_Obesity),
  area2 = sum(venn_binary_clean$High_Triglycerides),
  area3 = sum(venn_binary_clean$Low_HDL),
  area4 = sum(venn_binary_clean$High_Blood_Pressure),
  area5 = sum(venn_binary_clean$High_Fasting_Glucose),
  
  n12 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Triglycerides),
  n13 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$Low_HDL),
  n14 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Blood_Pressure),
  n15 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Fasting_Glucose),
  n23 = sum(venn_binary_clean$High_Triglycerides & venn_binary_clean$Low_HDL),
  n24 = sum(venn_binary_clean$High_Triglycerides & venn_binary_clean$High_Blood_Pressure),
  n25 = sum(venn_binary_clean$High_Triglycerides & venn_binary_clean$High_Fasting_Glucose),
  n34 = sum(venn_binary_clean$Low_HDL & venn_binary_clean$High_Blood_Pressure),
  n35 = sum(venn_binary_clean$Low_HDL & venn_binary_clean$High_Fasting_Glucose),
  n45 = sum(venn_binary_clean$High_Blood_Pressure & venn_binary_clean$High_Fasting_Glucose),
  n123 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Triglycerides & venn_binary_clean$Low_HDL),
  n124 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Triglycerides & venn_binary_clean$High_Blood_Pressure),
  n125 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Triglycerides & venn_binary_clean$High_Fasting_Glucose),
  n134 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$Low_HDL & venn_binary_clean$High_Blood_Pressure),
  n135 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$Low_HDL & venn_binary_clean$High_Fasting_Glucose),
  n145 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Blood_Pressure & venn_binary_clean$High_Fasting_Glucose),
  n234 = sum(venn_binary_clean$High_Triglycerides & venn_binary_clean$Low_HDL & venn_binary_clean$High_Blood_Pressure),
  n235 = sum(venn_binary_clean$High_Triglycerides & venn_binary_clean$Low_HDL & venn_binary_clean$High_Fasting_Glucose),
  n245 = sum(venn_binary_clean$High_Triglycerides & venn_binary_clean$High_Blood_Pressure & venn_binary_clean$High_Fasting_Glucose),
  n345 = sum(venn_binary_clean$Low_HDL & venn_binary_clean$High_Blood_Pressure & venn_binary_clean$High_Fasting_Glucose),
  n1234 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Triglycerides & venn_binary_clean$Low_HDL & venn_binary_clean$High_Blood_Pressure),
  n1235 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Triglycerides & venn_binary_clean$Low_HDL & venn_binary_clean$High_Fasting_Glucose),
  n1245 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Triglycerides & venn_binary_clean$High_Blood_Pressure & venn_binary_clean$High_Fasting_Glucose),
  n1345 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$Low_HDL & venn_binary_clean$High_Blood_Pressure & venn_binary_clean$High_Fasting_Glucose),
  n2345 = sum(venn_binary_clean$High_Triglycerides & venn_binary_clean$Low_HDL & venn_binary_clean$High_Blood_Pressure & venn_binary_clean$High_Fasting_Glucose),
  n12345 = sum(venn_binary_clean$Abdominal_Obesity & venn_binary_clean$High_Triglycerides & venn_binary_clean$Low_HDL & venn_binary_clean$High_Blood_Pressure & venn_binary_clean$High_Fasting_Glucose),
  
  category = c("Abdominal Obesity", "High Triglycerides", "Low HDL", 
               "High Blood Pressure", "High Fasting Glucose"),
  fill = c("blue", "yellow", "green", "red", "purple"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5)
)

# Visualizar o gráfico
grid.draw(venn_plot)



# Calcular a soma dos critérios atendidos
metadados.alpha.all <- metadados.alpha.all %>%
  mutate(
    Sindrome_Metabolica = rowSums(
      cbind(Obesidade_Abdominal, Triglicerideos_Altos, HDL_Baixo, Pressao_Alta, Glicemia_Alta)
    ) >= 3
  )


# Visualizar as primeiras linhas
head(metadados.alpha.all[, c("Obesidade_Abdominal", "Triglicerideos_Altos", "HDL_Baixo", 
                             "Pressao_Alta", "Glicemia_Alta", "Sindrome_Metabolica")])



table(metadados.alpha.all$Sindrome_Metabolica) / nrow(metadados.alpha.all) * 100


# Instale o pacote VennDiagram, caso ainda não tenha
if (!require("VennDiagram")) {
  install.packages("VennDiagram")
}
library(VennDiagram)

# Dados de exemplo: Número de indivíduos em cada conjunto
venn_counts <- list(
  "Abdominal Obesity" = c(1:50),
  "Glucose" = c(10:30),
  "Blood Pressure" = c(20:40),
  "HDL" = c(30:60),
  "Triglycerides" = c(15:35)
)

# Criando o diagrama
venn <- venn.diagram(
  x = venn_counts,
  filename = NULL,  # Para salvar como objeto e customizar depois
  fill = c("grey70", "grey70", "grey70", "grey70", "grey70"), # Cinza por padrão
  alpha = 0.5,
  lty = "dotted",  # Linha tracejada
  cex = 1.5,       # Tamanho do texto
  fontface = "bold",
  cat.cex = 1.2,   # Tamanho do texto das categorias
  cat.fontface = "bold",
  cat.pos = 0,     # Ajusta a posição dos rótulos
  cat.dist = 0.05,
  main = "Metabolic Syndrome",
  main.cex = 2
)

# Customização manual: Colorir áreas com 3 ou mais parâmetros
# Identifique os subconjuntos manualmente ou usando funções como `get.venn.partitions`

# Adicionar as cores onde há interseção de 3 ou mais
highlight_indices <- c(10, 15, 20) # Exemplo de índices onde deve colorir
for (i in highlight_indices) {
  venn[[i]]$gp$fill <- "red"
}

# Plotando o diagrama final
grid.newpage()
grid.draw(venn)


