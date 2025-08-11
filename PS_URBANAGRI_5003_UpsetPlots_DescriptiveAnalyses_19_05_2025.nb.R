library(ComplexUpset)
library(ggplot2)




###============= Grafico Upset ====================###


metadados_upset <- metadados.all.filtrado %>%
  mutate(
    High_Triglycerides = ifelse(!is.na(Triglycerides) & Triglycerides > 150, 1, 0),
    High_Fasting_Glucose = ifelse(!is.na(Glucose) & Glucose > 100, 1, 0),
    Obesity = ifelse(BMI >= 30, 1, 0),  # IMC tem 0 NA, sem ifelse extra
    Low_HDL = ifelse(!is.na(HDL) & HDL < 40, 1, 0),
    High_Blood_Pressure = ifelse(!is.na(Systolic) & !is.na(Diastolic) &
                                   (Systolic >= 130 | Diastolic >= 85), 1, 0)
  )





# Criar o gráfico com cores ajustadas
grafico_upset <- upset(
  metadados_upset,
  intersect = c('High_Triglycerides', 'High_Fasting_Glucose', 'Obesity', 'Low_HDL', 'High_Blood_Pressure'),
  name = 'Metabolic Risk Criteria',
  base_annotations = list(
    'Intersection size' = intersection_size(
      mapping = aes(fill = "azul"),
      text = list(size = 4)
    ) + scale_fill_manual(values = c("azul" = "#003399"))  # Azul forte
  ),
  set_sizes = upset_set_size(
    mapping = aes(fill = "vermelho")
  ) + scale_fill_manual(values = c("vermelho" = "#8B0000")),  # Vermelho escuro
  width_ratio = 0.2
)+
  labs(title = "Metabolic Health Profile of Farmers") +  # TÍTULO em inglês
  theme(legend.position = "none
        ")  # REMOVE a legenda 'fill'

ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/upset_metabolic_farmers.png", plot = grafico_upset, width = 10, height = 6, dpi = 300)


#========== Separados por idade ==============#

media_idade <- mean(metadados_upset$Age, na.rm = TRUE)

metadados_upset_jovens <- metadados_upset %>%
  filter(Age < media_idade)

metadados_upset_idosos <- metadados_upset %>%
  filter(Age > media_idade)


library(ComplexUpset)
library(ggplot2)

# Gráfico final = Jovens
u <-upset(
  metadados_upset_jovens,
  intersect = c('High_Triglycerides', 'High_Fasting_Glucose', 'Obesity', 'Low_HDL', 'High_Blood_Pressure'),
  name = 'Metabolic Risk Criteria',
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
  labs(title = "Metabolic Health Profile of Farmers Under the Mean Age (n= 60)")

ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/upset_metabolic_farmers_underMean_age.png", plot = u, width = 10, height = 6, dpi = 300)


# Gráfico final = idosos
u2 <-upset(
  metadados_upset_idosos,
  intersect = c('High_Triglycerides', 'High_Fasting_Glucose', 'Obesity', 'Low_HDL', 'High_Blood_Pressure'),
  name = 'Metabolic Risk Criteria',
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
  labs(title = "Metabolic Health Profile of Farmers Over the Mean Age (n= 68)")

ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/upset_metabolic_farmers_OverMean_age.png", plot = u2, width = 10, height = 6, dpi = 300)




#install.packages("writexl")  # só se ainda não tiver instalado
library(writexl)
write_xlsx(metadados_upset, path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados_upset.xlsx")






colunas_binarias <- c("High_Blood_Pressure", "High_Fasting_Glucose", "HbA1c_bin", "Obese",
                      "High_Trig", "High_Col", "Low_HDL", "High_LDL", "High_GGT")

# Gráfico final = idosos
u2 <-upset(
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
  labs(title = "Metabolic Health Profile")

print(u2)

ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/upset_metabolic_farmers_OverMean_age.png", plot = u2, width = 10, height = 6, dpi = 300)

