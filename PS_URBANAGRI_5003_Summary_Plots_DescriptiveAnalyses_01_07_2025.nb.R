#------------------------------------------------------------------------------#
#              title: "Summary, Plots and Descriptive Stats"
#------------------------------------------------------------------------------#
  

# Pacotes necessários
library(dplyr)
#install.packages("gt")
library(gt)
library(tidyr)
library(readr)
library(webshot2)
#webshot::install_phantomjs()
library(ggplot2)
#install.packages("pagedown")
#install.packages("xfun")
library(pagedown)
#install.packages("kableExtra")
library(kableExtra)
#install.packages("FSA")
library(FSA)









# Carregando os dados
metadados.all <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/metadados_130.csv") # Atualize com o caminho correto do arquivo

metadados.all.filtrado <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/metadadosN128_english.csv")

#tirar voluntaris S30092.F00 que tinha diagnostico de cancer e 40142 nao tem nenhum exame ou recordatorios
metadados.all.filtrado <- metadados.all[!metadados.all$Sample.id %in% c("S40142.F00", "S30092.F00"), ]

diet_data <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/diet_data.csv") # Atualize com o caminho correto do arquivo


library(dplyr)


#=======================================================#
#General Information and Anthropometric Markers
#=======================================================#



#write.csv(x = metadados.all.filtrado, 
 #         file = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/metadadosN128_english.csv", row.names = FALSE)

metadados.saude.alpha <- metadados.all.filtrado[, c(
  "Sample.id", "IL17A", "IFNGamma", "TNF", "IL10", "IL6", "IL4", "IL2",
  "Age", "Sex", "Systolic", "Diastolic", "Weight", "Height", "BMI",
  "Waist_Circumference", "Hip_Circumference", "HbA1c", "Cholesterol",
  "LDL", "HDL", "VLDL", "Triglycerides", "TGO", "TGP", "GGT",
  "Glucose", "Insulin", "CRP"
)]




summary(metadados.all.filtrado)






resumo_completo <- data.frame(
  Variavel = names(metadados.all.filtrado),
  Media = sapply(metadados.all.filtrado, function(x) if(is.numeric(x)) mean(x, na.rm = TRUE) else NA),
  DesvioPadrao = sapply(metadados.all.filtrado, function(x) if(is.numeric(x)) sd(x, na.rm = TRUE) else NA),
  Minimo = sapply(metadados.all.filtrado, function(x) if(is.numeric(x)) min(x, na.rm = TRUE) else NA),
  Maximo = sapply(metadados.all.filtrado, function(x) if(is.numeric(x)) max(x, na.rm = TRUE) else NA)
)

# Visualizar
head(resumo_completo, 10)

# (opcional) salvar como planilha
#write.csv(resumo_completo, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/resumo_numerico.csv", row.names = FALSE)

#idade homens e mulheres


metadados.all.filtrado %>%
  group_by(Sex) %>%
  summarise(
    Mean = mean(Age, na.rm = TRUE),
    SD = sd(Age, na.rm = TRUE),
    Min = min(Age, na.rm = TRUE),
    Max = max(Age, na.rm = TRUE)
  )



#tabela com dados basicos



# Suponha que df_sex_age seja o seu resumo:
df_sex_age <- metadados.all.filtrado %>%
  group_by(Sex) %>%
  summarise(
    Count = n(),
    `Mean Age` = round(mean(Age, na.rm = TRUE), 2),
    `Standard Deviation` = round(sd(Age, na.rm = TRUE), 2),
    `Minimum Age` = min(Age, na.rm = TRUE),
    `Maximum Age` = max(Age, na.rm = TRUE)
  )



df_sex_age %>%
  gt() %>%
  tab_header(
    title = md("**General Information (n=128)**")
  ) %>%
  cols_label(
    Sex = "Sex",
    Count = "Count",
    `Mean Age` = "Mean Age",
    `Standard Deviation` = "St Dev",
    `Minimum Age` = "Min Age",
    `Maximum Age` = "Max Age"
  )




summary_stats <- metadados.all.filtrado %>%
  select(where(is.numeric)) %>%
  pivot_longer(everything(), names_to = "Variável", values_to = "Valor") %>%
  group_by(Variável) %>%
  summarise(
    Mean = round(mean(Valor, na.rm = TRUE), 2),
    SD = round(sd(Valor, na.rm = TRUE), 2),
    Min = min(Valor, na.rm = TRUE),
    Max = max(Valor, na.rm = TRUE),
    N = sum(!is.na(Valor)),
    .groups = "drop"
  )
View(summary_stats)

write.csv(summary_stats, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/summary_stats.csv", row.names = FALSE)





# Lê a tabela corrigida
tabela_info <- read_csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/tabela_info_corrigida.csv", 
                        col_types = cols(.default = col_character()))

tabela_gt <- tabela_info %>%
  gt() %>%
  tab_header(title = md("**General Information and Anthropometric Markers**")) %>%
  cols_label(
    Characteristic = "Characteristic",
    Mean = "Mean",
    SD = "SD",
    Min = "Min",
    Max = "Max"
  ) %>%
  opt_table_font(font = list(
    google_font("Arial"),
    default_fonts()
  )) %>%
  fmt_number(
    columns = c(Mean, SD, Min, Max),
    decimals = 2,
    rows = !grepl("/", tabela_info$Mean)
  ) %>%
  fmt_missing(columns = everything(), missing_text = "") %>%
  cols_width(
    Characteristic ~ px(220),
    everything() ~ px(70)
  ) %>%
  tab_options(
    table.font.size = px(14),
    column_labels.font.weight = "bold",
    data_row.padding = px(6),
    table.width = pct(100)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = "Characteristic",
      rows = Characteristic %in% c(
        "Characteristic", 
        "GENERAL INFORMATION (n=128)",
        "ANTHROPOMETRIC MARKERS (n=127)"
      )
    )
  )


# Para salvar como imagem:
gtsave(tabela_gt, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/tabela_general_info.png", zoom = 2, expand = 10)

#salvar em alta. zoom = 2 aumenta a resolução expand = 5 adiciona margem ao redor da imagem
gtsave(tabela_gt, "tabela_info.png", zoom = 2, expand = 5)








#=======================================================#
#              Inflamatory Markers
#=======================================================#



# Lista dos marcadores inflamatórios
inflamatorios <- c("IL17A", "IFNGamma", "TNF", "IL10", "IL6", "IL4", "IL2", "PCR")

# Calcular resumo estatístico
summary_inflamatorios <- metadados.all.filtrado %>%
  select(any_of(inflamatorios)) %>%
  summarise(across(
    everything(),
    list(
      Mean = ~mean(., na.rm = TRUE),
      SD = ~sd(., na.rm = TRUE),
      Min = ~min(., na.rm = TRUE),
      Max = ~max(., na.rm = TRUE),
      N = ~sum(!is.na(.))
    ),
    .names = "{.col}_{.fn}"
  )) %>%
  pivot_longer(
    everything(),
    names_to = c("Marker", ".value"),
    names_sep = "_"
  ) %>%
  select(Marker, Mean, SD, Min, Max, N)

# Ver resultado
print(summary_inflamatorios)

write.csv(summary_inflamatorios, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/tabela_inflammatory_markers.csv", row.names = FALSE)


# Lê a planilha
tabela_infla <- read_csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/tabela_inflammatory_markers.csv", col_types = cols(.default = col_character()))

# Converte colunas numéricas para número real
tabela_infla <- tabela_infla %>%
  mutate(across(c(Mean, SD, Min, Max), as.numeric))


tabela_gt_infla <- tabela_infla %>%
  gt() %>%
  tab_header(title = md("**Inflammatory Markers (n=92)**")) %>%
  cols_label(
    Marker = "Marker",
    Mean = "Mean",
    SD = "SD",
    Min = "Min",
    Max = "Max",
    N = "N"
  ) %>%
  opt_table_font(font = list(
    google_font("Arial"),
    default_fonts()
  )) %>%
  fmt_number(
    columns = c(Mean, SD, Min, Max),
    decimals = 2
  ) %>%
  fmt_missing(columns = everything(), missing_text = "") %>%
  cols_width(
    Marker ~ px(160),
    everything() ~ px(70)
  ) %>%
  tab_options(
    table.font.size = px(14),
    column_labels.font.weight = "bold",
    data_row.padding = px(6),
    table.width = pct(100)
  )



# Salvar como PNG em alta resolução
gtsave(
  data = tabela_gt_infla,  # substitua pelo nome do seu objeto gt
  filename = "tabela_inflammatory_markers.png",
  path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/",
  expand = 10          # aumenta a resolução (quanto maior, mais nítido)
)



#=======================================================#
#                        Diet
#=======================================================#
#"BHEI_R_Score_Total", "Percentual_NOVA_group_1", "Percentual_NOVA_group_2", "Percentual_NOVA_group_3" 



metadados.dieta.all <- diet_data[, c("Sample.id","carboidrato_total_g", "proteina_g", "lipidios_g", "fibra_alimentar_g",   "colesterol_mg", "acidos_graxos_saturados_g", "acidos_graxos_monoinsaturados_g", "acidos_graxos_poliinsaturados_g", "calcio_mg", "ferro_mg", "sodio_mg", "magnesio_mg", "fosforo_mg",  "potassio_mg", "manganes_mg", "zinco_mg", "cobre_mg", "selenio_mcg", "vitamina_A_RAE_mcg", "vitamina_D_mcg",  "vitamina_E_mg", "tiamina_mg", "riboflavina_mg",  "niacina_mg", "vitamina_C_mg",  "equivalente_de_folato_mcg", "sal_de_adicao_g", "acucar_de_adicao_g"          
)]



summary_dieta <- as.data.frame (summary (metadados.dieta.all))

as.data.frame(summary_dieta)

#===========================================#
#                Region
#===========================================#



region_counts <- metadados.all.filtrado %>%
  count(Region, name = "N")

#mudar nomes pro ingles

region_counts <- region_counts %>%
  mutate(Region_EN = case_when(
    Region == "Leste" ~ "East",
    Region == "Leste Adjacente" ~ "Adjacent East",
    Region == "Norte" ~ "North",
    Region == "Sul" ~ "South",
    Region == "Sul Adjacente" ~ "Adjacent South",
    TRUE ~ Region  # mantém original se não houver correspondência
  ))

#S40121 esta na regiao errada
metadados.all.filtrado <- metadados.all.filtrado %>%
  mutate(Region = ifelse(Sample.id == "S40121.F00", "Sul Adjacente", Region))

metadados.all <- metadados.all %>%
  mutate (Region = ifelse(Sample.id == "S40121.F00", "Sul Adjacente", Region))


grafico_regiao <-ggplot(region_counts, aes(x = reorder(Region_EN, -N), y = N, fill = Region_EN)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of Participants by Region",
       x = "Region",
       y = "Number of Participants") +
  theme_minimal() +
  theme(legend.position = "none")


# Salvar como PNG em alta resolução
ggsave(
  filename = "grafico_regiao.png",
  plot = grafico_regiao,  # seu objeto ggplot
  path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/",
  width = 8,
  height = 6,
  dpi = 300
)



# Calculando características da população (Mean, SD, Min e Max)
caracteristicas <- data.frame(
  Variable = c(
    "Fasting glucose (mmol/L)",
    "Fasting insulin (μU/mL) (n=73)", "HbA1c (%)", "Total cholesterol (mmol/L)",
    "LDL-C (mmol/L)", "HDL-C (mmol/L)", "Triglycerides (mmol/L)",
    "Systolic blood pressure (mmHg)", "Diastolic blood pressure (mmHg)",
    "CRP (mg/L) (n=73)"
  ),
  Mean = c(
    round(mean(metadados.all.filtrado$Glucose, na.rm = TRUE), 2),
    round(mean(metadados.all.filtrado$Insulin, na.rm = TRUE), 2),
    round(mean(metadados.all.filtrado$HbA1c, na.rm = TRUE), 2),
    round(mean(metadados.all.filtrado$Cholesterol, na.rm = TRUE), 2),
    round(mean(metadados.all.filtrado$LDL, na.rm = TRUE), 2),
    round(mean(metadados.all.filtrado$HDL, na.rm = TRUE), 2),
    round(mean(metadados.all.filtrado$Triglycerides, na.rm = TRUE), 2),
    round(mean(metadados.all.filtrado$Systolic, na.rm = TRUE), 2),
    round(mean(metadados.all.filtrado$Diastolic, na.rm = TRUE), 2),
    round(mean(metadados.all.filtrado$CRP, na.rm = TRUE), 2)
  ),
  SD = c(
    round(sd(metadados.all.filtrado$Glucose, na.rm = TRUE), 2),
    round(sd(metadados.all.filtrado$Insulin, na.rm = TRUE), 2),
    round(sd(metadados.all.filtrado$HbA1c, na.rm = TRUE), 2),
    round(sd(metadados.all.filtrado$Cholesterol, na.rm = TRUE), 2),
    round(sd(metadados.all.filtrado$LDL, na.rm = TRUE), 2),
    round(sd(metadados.all.filtrado$HDL, na.rm = TRUE), 2),
    round(sd(metadados.all.filtrado$Triglycerides, na.rm = TRUE), 2),
    round(sd(metadados.all.filtrado$Systolic, na.rm = TRUE), 2),
    round(sd(metadados.all.filtrado$Diastolic, na.rm = TRUE), 2),
    round(sd(metadados.all.filtrado$CRP, na.rm = TRUE), 2)
  ),
  Min = c(
    round(min(metadados.all.filtrado$Glucose, na.rm = TRUE), 2),
    round(min(metadados.all.filtrado$Insulin, na.rm = TRUE), 2),
    round(min(metadados.all.filtrado$HbA1c, na.rm = TRUE), 2),
    round(min(metadados.all.filtrado$Cholesterol, na.rm = TRUE), 2),
    round(min(metadados.all.filtrado$LDL, na.rm = TRUE), 2),
    round(min(metadados.all.filtrado$HDL, na.rm = TRUE), 2),
    round(min(metadados.all.filtrado$Triglycerides, na.rm = TRUE), 2),
    round(min(metadados.all.filtrado$Systolic, na.rm = TRUE), 2),
    round(min(metadados.all.filtrado$Diastolic, na.rm = TRUE), 2),
    round(min(metadados.all.filtrado$CRP, na.rm = TRUE), 2)
    
  ),
  Max = c(
    round(max(metadados.all.filtrado$Glucose, na.rm = TRUE), 2),
    round(max(metadados.all.filtrado$Insulin, na.rm = TRUE), 2),
    round(max(metadados.all.filtrado$HbA1c, na.rm = TRUE), 2),
    round(max(metadados.all.filtrado$Cholesterol, na.rm = TRUE), 2),
    round(max(metadados.all.filtrado$LDL, na.rm = TRUE), 2),
    round(max(metadados.all.filtrado$HDL, na.rm = TRUE), 2),
    round(max(metadados.all.filtrado$Triglycerides, na.rm = TRUE), 2),
    round(max(metadados.all.filtrado$Systolic, na.rm = TRUE), 2),
    round(max(metadados.all.filtrado$Diastolic, na.rm = TRUE), 2),
    round(max(metadados.all.filtrado$CRP, na.rm = TRUE), 2)
    
  )
)



# Visualizando o resultado final
print(caracteristicas)

# Salvando como CSV
write.csv(caracteristicas, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/population_characteristics_final.csv", row.names = FALSE)





#Fazer uma tabela com formatação bonita


column_labels.font.weight = "bold"


# Criar a tabela bonitinha com gt
TableClinicalTraits <- caracteristicas %>%
  gt() %>%
  tab_header(
    title = "Clinical Traits",
    subtitle = "(n=117)"
  ) %>%
  fmt_number(
    columns = c(Mean, SD),
    decimals = 2
  ) %>%
  cols_label(
    Variable = "Variable",
    Mean = "Mean",
    SD = "SD",
    Min = "Min",
    Max = "Max"
  ) %>%
  cols_width(
    Variable ~ px(250),
    everything() ~ px(90)
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.align = "center",
    table.width = pct(100),
    column_labels.font.weight = "bold"  # ⬅️ aqui ativa o negrito do cabeçalho
  )



# Exibir a tabela
TableClinicalTraits



gtsave(TableClinicalTraits, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/TableClinicalTraits.html")

gtsave(
  data = TableClinicalTraits,
  filename = "TableClinicalTraits.png",
  path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/",  # ou outro caminho da sua preferência
  expand = 10,   # aumenta a resolução do fundo
  zoom = 2       # dobra a nitidez da tabela
)





#Tabela com idade, IMC e sexo


#Kruskal-Wallis para mais de dois grupos
#Para comparar IMC ou Idade entre Sul, Norte, Leste, etc.:


# Kruskal-Wallis para IMC
kruskal.test(BMI ~ Region, data = metadados.all.filtrado)



# Kruskal-Wallis para Idade


kruskal.test(Age ~ Region, data = metadados.all.filtrado)



#Se o resultado do Kruskal-Wallis for significativo, você pode fazer comparações par-a-par usando Dunn's Test:


# Comparações par-a-par
dunnTest(BMI ~ Region, data = metadados.all.filtrado, method = "bonferroni")


#Mann-Whitney U para dois grupos (Rural vs Urbano)
#Para comparar IMC ou Idade entre Rural e Urbano:

# Mann-Whitney para IMC
wilcox.test(BMI ~ Region_type, data = metadados.all.filtrado)



```{r}
# Mann-Whitney para Idade
wilcox.test(Age ~ Region_type, data = metadados)


#Qui-quadrado ou Fisher para variáveis categóricas (Sexo)
#Se você quiser comparar Sexo entre regiões ou tipos de região:

#Para mais de dois grupos (Sul, Norte, etc.):


# Tabela cruzada para Sexo e Região
tabela_sexo <- table(metadados.all.filtrado$Sex, metadados.all.filtrado$Region)

# Qui-quadrado
chisq.test(tabela_sexo)



#Para dois grupos (Rural vs Urbano):

# Tabela cruzada para Sexo e Tipo de Região
tabela_sexo <- table(metadados.all.filtrado$Sex, metadados.all.filtrado$Region_type)

# Teste exato de Fisher (para pequenas amostras)
fisher.test(tabela_sexo)

### IMPORTANTE PARA ESCOLHER OS TESTES ###
#Resumo de Como Escolher o Teste

#Variável: IMC, Idade	Grupos: 2 grupos (Rural/Urbano)	=> Mann-Whitney U

#Variável: IMC, Idade	Grupos: Mais de 2 grupos =>	Kruskal-Wallis

#Variável: Sexo (categórica)	Grupos: 2 grupos (Rural/Urbano)	=> Fisher ou Qui-quadrado

#Variável: Sexo (categórica)	Grupos: Mais de 2 grupos	=> Qui-quadrado


#Não deu diferença estatística!


#Alterar a paleta de cores. Se quiser mudar as cores, substitua scale_fill_brewer(palette = "Set3") por outra paleta disponível, como: "Paired", "Pastel1", "Dark2" Ou use cores personalizadas: 

boxplot_idade <- ggplot(metadados.all.filtrado, aes(x = Region, y = Age, fill = Region)) +
  geom_boxplot() +
  labs(title = "Age vs Region",
       x = "Region",
       y = "Age") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")  # Paleta de cores



# Salvar como PNG
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/boxplot_idade.png", plot = boxplot_idade, width = 8, height = 6, dpi = 300)

# Salvar como PDF
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/boxplot_idade.pdf", plot = boxplot_idade, width = 8, height = 6)




#Tabela com : Idade, Sexo, Média de Idade entre os sexos e Número Total. 
#Criar a tabela resumida
#Usaremos o dplyr para calcular os resumos e o gt para formatar a tabela.



# Ajustar os rótulos para inglês (Male e Female)
tabela_resumo <- metadados.all.filtrado %>%
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
gtsave(tabela_final, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/age_sex_summary.html")

# Salvar como PDF
gtsave(tabela_final, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/age_sex_summary.jpg")

gtsave(
  tabela_final, 
  "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/age_sex_summary.png"
)




### Parasitológico. ###
#Tabela que mostre quantos voluntários têm resultado positivo no parasitológico e listar os tipos de parasitas encontrados





# Substituir "e" por "and" na coluna Parasitas
metadados <- metadados_comp_Corrigido_130 %>%
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
  group_by("Parasitologico +") %>%
  summarise(N = n()) %>%
  arrange(desc(N))

# Exibir o resumo
print(resumo_parasitologico)


# Substituir "e" por "and" na coluna de tipos de parasitas
tabela_parasitologico <- metadados %>%
  filter(`Parasitologico +` == "Sim") %>%  # Selecionar apenas positivos no parasitológico
  group_by(Parasitas) %>%  # Agrupar por tipos de parasitas
  summarise(`Número de casos` = n()) %>%  # Contar ocorrências de cada parasita
  arrange(desc(`Número de casos`)) %>%  # Ordenar por número de casos
  ungroup()

tabela_final_parasitas <- tabela_parasitologico %>%
  gt() %>%
  tab_header(
    title = md("**Parasitological Test (n=130)**")
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

gtsave(tabela_final_parasitas, 
       filename = "tabela_parasitologico.png", 
       path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/",  # Ajuste o caminho se quiser outro
       expand = 10,  # aumenta a resolução
       vwidth = 1000,  # largura em pixels
       vheight = 800   # altura em pixels
)

