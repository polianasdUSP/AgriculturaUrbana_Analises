#=========================================#
  # title: "Alpha Diversity Analisys"
#=========================================#

  

# Carregar as bibliotecas necessárias
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggcorrplot)
library(pheatmap)
library(grid)
library(patchwork)
library(ggpubr)
library(ggpmisc)
library(microbiome)



# Substituir valores 9 por 90 na coluna Diastolic
metadados.shannon$Diastolic <- ifelse(metadados.shannon$Diastolic == 9, 90, metadados.shannon$Diastolic)

# Substituir valores 9 por 90 na coluna Diastolic
metadados.alpha.all$Diastolic <- ifelse(metadados.alpha.all$Diastolic == 9, 90, metadados.alpha.all$Diastolic)

#criar uma coluna com IMC

# Adicionando uma coluna de IMC ao metadados
metadados.shannon <- metadados.shannon %>%
  mutate(
    IMC = Weight / (Height^2) # Cálculo do IMC
  )

# Visualizando as primeiras linhas para verificar o resultado
head(metadados.shannon$BMI)

# Renomear a coluna "IMC" para "BMI" no objeto metadados.shannon
colnames(metadados.shannon)[colnames(metadados.shannon) == "IMC"] <- "BMI"


# Adicionando uma coluna de Waist-to_Hip Ratio ao metadados.shannon
metadados.shannon <- metadados.shannon %>%
  mutate(
    WHR = Waist / Hip # Cálculo do IMC
  )



# Visualizando as primeiras linhas para verificar o resultado
head(metadados.shannon$WHR)

# Cálculo do índice TyG
metadados.shannon$TyG <- log(metadados.shannon$TRIGLICERIDES * metadados.shannon$GLICOSE) / 2

# Visualizar os primeiros valores calculados
head(metadados.shannon$TyG)

# Resumo estatístico da nova coluna TyG
summary(metadados.shannon$TyG)


# Selecionar as colunas numéricas de interesse e renomeá-las para metadados_shannon_selected
metadados_shannon_blood_work <- metadados.shannon %>%
  select(shannon_entropy, IL17A, IFNGamma, TNF, IL10, IL6, IL4, IL2, 
         Age, Systolic, Diastolic, Weight, Waist, HbA1c, COLESTEROL, 
         LDL, HDL, VLDL, TRIGLICERIDES, TGO, GGT, TGP, INSULINA, 
         HOMA, PCR, TyG, BMI, WHR)

colnames(metadados_shannon_blood_work)

# Função para calcular os p-valores das correlações
cor.mtest_shannon <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(p.mat)
}


# Calcular a correlação de Pearson entre as variáveis selecionadas
cor_matrix_shannon <- cor(metadados_shannon_blood_work, method = "spearman", use = "pairwise.complete.obs")

# Calcular os p-valores
p_values_shannon <- cor.mtest_shannon(metadados_shannon_blood_work, conf.level = 0.95)

# Criar uma matriz de asteriscos para significância
asterisks_shannon <- ifelse(p_values_shannon < 0.001, "***", 
                            ifelse(p_values_shannon < 0.01, "**", 
                                   ifelse(p_values_shannon < 0.05, "*", "")))


#Definir os limites para centralizar o zero na escala de cores
breaks <- seq(-1, 1, length.out = 51) # Define uma escala de -1 a 1 com 50 intervalos


# Gerar o heatmap com o pheatmap
pheatmap(cor_matrix_shannon, 
         display_numbers = asterisks_shannon, # Mostra os asteriscos no mapa
         color = colorRampPalette(c("blue", "white", "red"))(50), # Escala de cores com branco no centro
         breaks = breaks, # Aplica os limites para centralizar o zero
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         main = "Shannon Entropy, Inflamatory and Metabolic Parameters")





#===== ALPHA.ALL ====#

# Adicionando uma coluna de IMC ao metadados
metadados.alpha.all <- metadados.alpha.all %>%
  mutate(
    BMI = Weight / (Height^2) # Cálculo do IMC
  )

# Visualizando as primeiras linhas para verificar o resultado
head(metadados.alpha.all$BMI)


# Adicionando uma coluna de Waist-to_Hip Ratio ao metadados.shannon
metadados.alpha.all <- metadados.alpha.all %>%
  mutate(
    WHR = Waist / Hip # Cálculo do Waist-to-hip_ratio
  )



# Visualizando as primeiras linhas para verificar o resultado
head(metadados.alpha.all$WHR)

# Cálculo do índice TyG
metadados.alpha.all$TyG <- log(metadados.alpha.all$TRIGLICERIDES * metadados.alpha.all$GLICOSE) / 2

# Visualizar os primeiros valores calculados
head(metadados.alpha.all$TyG)

# Resumo estatístico da nova coluna TyG
summary(metadados.alpha.all$TyG)


# Selecionar as colunas numéricas de interesse e renomeá-las para metadados_faith_pd_selected
metadados_alpha_all_blood_work <- metadados.alpha.all %>%
  select(shannon_entropy, observed_features, pielou_evenness,faith_pd, IL17A, IFNGamma, TNF, IL10, IL6, IL4, IL2, 
         Age, Systolic, Diastolic, Weight, Waist, HbA1c, COLESTEROL, 
         LDL, HDL, VLDL, TRIGLICERIDES, TGO, GGT, TGP, INSULINA, 
         HOMA, PCR, TyG, BMI, WHR)


# Função para calcular os p-valores das correlações
cor.mtest_alpha <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(p.mat)
}

# Calcular a correlação de Pearson entre as variáveis selecionadas
cor_matrix_alpha <- cor(metadados_alpha_all_blood_work, method = "spearman", use = "pairwise.complete.obs")

# Calcular os p-valores
p_values_alpha <- cor.mtest_alpha(metadados_alpha_all_blood_work, conf.level = 0.95)

# Criar uma matriz de asteriscos para significância
asterisks_alpha <- ifelse(p_values_alpha < 0.001, "***", 
                          ifelse(p_values_alpha < 0.01, "**", 
                                 ifelse(p_values_alpha < 0.05, "*", "")))


#Definir os limites para centralizar o zero na escala de cores
breaks <- seq(-1, 1, length.out = 51) # Define uma escala de -1 a 1 com 50 intervalos

# Gerar o heatmap com o pheatmap
pheatmap(cor_matrix_alpha, 
         display_numbers = asterisks_alpha, # Mostra os asteriscos no mapa
         color = colorRampPalette(c("blue", "white", "red"))(50),
         breaks = breaks, # Aplica os limites para centralizar o zero
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         main = "Alpha Diversity, Inflamatory and Metabolic Parameters")





#==== Parasitologico X Alpha ==========#


# Selecionar as colunas de interesse
metadados_alpha_Parasito <- metadados.alpha.all %>%
  select(ParasitologicoPositivo, pielou_evenness, faith_pd, observed_features, shannon_entropy)




# Realizar o Teste de Wilcoxon para cada índice de diversidade alfa

# Pielou Evenness
wilcox_pielou <- wilcox.test(pielou_evenness ~ ParasitologicoPositivo, data = metadados_alpha_Parasito)
p_value_pielou <- wilcox_pielou$p.value

# Faith PD
wilcox_faith <- wilcox.test(faith_pd ~ ParasitologicoPositivo, data = metadados_alpha_Parasito)
p_value_faith <- wilcox_faith$p.value

# Observed Features
wilcox_observed <- wilcox.test(observed_features ~ ParasitologicoPositivo, data = metadados_alpha_Parasito)
p_value_observed <- wilcox_observed$p.value

# Shannon Entropy
wilcox_shannon <- wilcox.test(shannon_entropy ~ ParasitologicoPositivo, data = metadados_alpha_Parasito)
p_value_shannon <- wilcox_shannon$p.value

# Exibir os p-valores
p_value_pielou



# Gráfico para pielou_evenness
plot_pielou <- ggplot(metadados_alpha_Parasito, aes(x = ParasitologicoPositivo, y = pielou_evenness)) +
  geom_boxplot(aes(fill = ParasitologicoPositivo), position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Pielou Evenness vs Parasitologico", y = "Pielou Evenness", x = "Parasitologico Positivo") +
  theme_minimal() +
  annotate("text", x = 1.5, y = max(metadados_alpha_Parasito$pielou_evenness, na.rm = TRUE), 
           label = paste("p =", round(p_value_pielou, 4)), size = 5)

print(plot_pielou)


# Gráfico para faith_pd
plot_faith <- ggplot(metadados_alpha_Parasito, aes(x = ParasitologicoPositivo, y = faith_pd)) +
  geom_boxplot(aes(fill = ParasitologicoPositivo), position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Faith PD vs Parasitologico", y = "Faith PD", x = "Parasitologico Positivo") +
  theme_minimal() +
  annotate("text", x = 1.5, y = max(metadados_alpha_Parasito$faith_pd, na.rm = TRUE), 
           label = paste("p =", round(p_value_faith, 4)), size = 5)

print(plot_faith)


# Gráfico para observed_features
plot_observed <- ggplot(metadados_alpha_Parasito, aes(x = ParasitologicoPositivo, y = observed_features)) +
  geom_boxplot(aes(fill = ParasitologicoPositivo), position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Richness vs Parasitologico", y = "Richness", x = "Parasitologico Positivo") +
  theme_minimal() +
  annotate("text", x = 1.5, y = max(metadados_alpha_Parasito$observed_features, na.rm = TRUE), 
           label = paste("p =", round(p_value_observed, 4)), size = 5)

print(plot_observed)

# Gráfico para shannon_entropy
plot_shannon <- ggplot(metadados_alpha_Parasito, aes(x = ParasitologicoPositivo, y = shannon_entropy)) +
  geom_boxplot(aes(fill = ParasitologicoPositivo), position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Shannon Entropy vs Parasitologico", y = "Shannon Entropy", x = "Parasitologico Positivo") +
  theme_minimal() +
  annotate("text", x = 1.5, y = max(metadados_alpha_Parasito$shannon_entropy, na.rm = TRUE), 
           label = paste("p =", round(p_value_shannon, 4)), size = 5)

print(plot_shannon)



# Organizar os gráficos usando ggarrange
combined_plot <- ggarrange(plot_pielou, plot_faith, plot_observed, plot_shannon,
                           labels = c("A", "B", "C", "D"),
                           ncol = 2, nrow = 2, 
                           common.legend = TRUE, legend = "right")

# Exibir o gráfico combinado
print(combined_plot)




#========= TREPONEMA x ALPHA =============#

# Calcular os p-valores usando o teste de Wilcoxon para a coluna Treponema
p_value_pielou <- wilcox.test(pielou_evenness ~ Treponema, data = metadados.alpha.all)$p.value

p_value_faith <- wilcox.test(faith_pd ~ Treponema, data = metadados.alpha.all)$p.value

p_value_observed <- wilcox.test(observed_features ~ Treponema, data = metadados.alpha.all)$p.value

p_value_shannon <- wilcox.test(shannon_entropy ~ Treponema, data = metadados.alpha.all)$p.value


# Gráfico para pielou_evenness

# Gráfico para pielou_evenness com "Presença de Treponema"
plot_pielou_trep <- ggplot(metadados.alpha.all, aes(x = factor(Treponema, levels = c(0, 1), labels = c("Ausente", "Presente")), y = pielou_evenness)) +
  geom_boxplot(aes(fill = factor(Treponema, levels = c(0, 1), labels = c("Ausente", "Presente"))), position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7, color = "black") +
  labs(title = "Pielou Evenness vs Presença de Treponema", y = "Pielou Evenness", x = "Presença de Treponema") +
  theme_minimal() +
  scale_fill_discrete(name = "Presença de Treponema") +
  annotate("text", x = 1.5, y = max(metadados.alpha.all$pielou_evenness, na.rm = TRUE) * 1.05, 
           label = paste("p =", round(p_value_pielou, 4)), size = 5)



print (plot_pielou_trep)


# Gráfico para faith_pd
plot_faith_trep <- ggplot(metadados.alpha.all, aes(x = factor(Treponema, levels = c(0, 1), labels = c("Ausente", "Presente")), y = faith_pd)) +
  geom_boxplot(aes(fill = factor(Treponema, levels = c(0, 1), labels = c("Ausente", "Presente"))), position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7, color = "black") +
  labs(title = "Faith PD vs Presença de Treponema", y = "Faith PD", x = "Presença de Treponema") +
  theme_minimal() +
  scale_fill_discrete(name = "Presença de Treponema") +
  annotate("text", x = 1.5, y = max(metadados.alpha.all$faith_pd, na.rm = TRUE) * 1.05, 
           label = paste("p =", round(p_value_faith, 4)), size = 5)

print (plot_faith_trep)

# Gráfico para observed_features
plot_observed_trep <- ggplot(metadados.alpha.all, aes(x = factor(Treponema, levels = c(0, 1), labels = c("Ausente", "Presente")), y = observed_features)) +
  geom_boxplot(aes(fill = factor(Treponema, levels = c(0, 1), labels = c("Ausente", "Presente"))), position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7, color = "black") +
  labs(title = "Richness vs Presença de Treponema", y = "Richness", x = "Presença de Treponema") +
  theme_minimal() +
  scale_fill_discrete(name = "Presença de Treponema") +
  annotate("text", x = 1.5, y = max(metadados.alpha.all$observed_features, na.rm = TRUE) * 1.05, 
           label = paste("p =", round(p_value_observed, 4)), size = 5)

print (plot_observed_trep)


# Gráfico para shannon_entropy
plot_shannon_trep <- ggplot(metadados.alpha.all, aes(x = factor(Treponema, levels = c(0, 1), labels = c("Ausente", "Presente")), y = shannon_entropy)) +
  geom_boxplot(aes(fill = factor(Treponema, levels = c(0, 1), labels = c("Ausente", "Presente"))), position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7, color = "black") +
  labs(title = "Shannon Entropy vs Presença de Treponema", y = "Shannon Entropy", x = "Presença de Treponema") +
  theme_minimal() +
  scale_fill_discrete(name = "Presença de Treponema") +
  annotate("text", x = 1.5, y = max(metadados.alpha.all$shannon_entropy, na.rm = TRUE) * 1.05, 
           label = paste("p =", round(p_value_shannon, 4)), size = 5)


print(plot_shannon_trep)

# Organizar os gráficos em uma única imagem
combined_plot_trep <- plot_pielou_trep + plot_faith_trep + plot_observed_trep + plot_shannon_trep + 
  plot_layout(ncol = 2)  # Organiza em 2 colunas

print(combined_plot_trep)



#=========== Parasitologico e Treponema ===============#



# Realizar o teste exato de Fisher diretamente
resultado_fisher <- fisher.test(table(metadados.alpha.all$Treponema, metadados.alpha.all$ParasitologicoPositivo))

# Exibir o resultado do teste exato de Fisher
print(resultado_fisher)


library(ggplot2)

# Criar o gráfico de barras diretamente
grafico_fisher <- ggplot(metadados.alpha.all, aes(x = factor(Treponema, levels = c(0, 1), labels = c("Ausente", "Presente")), 
                                          fill = factor(ParasitologicoPositivo, levels = c("Nao", "Sim")))) +
  geom_bar(position = "dodge") +
  labs(title = "Associação entre Treponema e Parasitológico Positivo",
       subtitle = paste("Teste Exato de Fisher\np-valor =", round(resultado_fisher$p.value, 4),
                        "\nOdds Ratio =", round(resultado_fisher$estimate, 4),
                        "\nIC 95%:", round(resultado_fisher$conf.int[1], 4), "-", round(resultado_fisher$conf.int[2], 4)),
       x = "Presença de Treponema",
       y = "Frequência") +
  theme_minimal() +
  scale_fill_discrete(name = "Parasitológico Positivo")

# Exibir o gráfico
print(grafico_fisher)



#========= DIETA x ALPHA ==============#



saude_dieta <- merge(saude_dieta, alpha.all, by.x = "IdVoluntario", by.y = "Row.names", all = TRUE)


colnames(saude_dieta)

saude_dieta <- merge(saude_dieta, metadados, by.x = "IdVoluntario", by.y = "Sample.id", all.x = TRUE)

saude_dieta <- merge(saude_dieta, metadados, by.x = "IdVoluntario", by.y = "Sample.id", all.x = TRUE)

# Selecionar as colunas numéricas de interesse e renomeá-las para metadados_shannon_selected
metadados.dieta.alpha <- saude_dieta %>%
  select(shannon_entropy, pielou_evenness, faith_pd, observed_features, residual_energia2_kcal, residual_carboidrato_total_g, residual_proteina_g,residual_lipidios_g, residual_fibra_alimentar_g, residual_colesterol_mg, residual_acidos_graxos_saturados_g, residual_acidos_graxos_monoinsaturados_g, residual_acidos_graxos_poliinsaturados_g, residual_acidos_graxos_trans_g, residual_calcio_mg, residual_ferro_mg, residual_sodio_mg, residual_magnesio_mg, residual_fosforo_mg, residual_potassio_mg, residual_manganes_mg, residual_zinco_mg, residual_cobre_mg, residual_selenio_mcg, residual_vitamina_A_RE_mcg, residual_vitamina_A_RAE_mcg, residual_vitamina_D_mcg, residual_vitamina_E_mg, residual_tiamina_mg, residual_riboflavina_mg, residual_niacina_mg, residual_vitamina_B6_mg, residual_vitamina_B12_mcg, residual_vitamina_C_mg, residual_equivalente_de_folato_mcg, residual_sal_de_adicao_g, residual_acucar_de_adicao_g)

```




```{r}
# Função para calcular os p-valores das correlações
cor.mtest_shannon <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(p.mat)
}


# Calcular a correlação de Pearson entre as variáveis selecionadas
cor_matrix_alpha.diet <- cor(metadados.dieta.alpha, method = "spearman", use = "pairwise.complete.obs")

# Calcular os p-valores
p_values_alpha.diet <- cor.mtest_shannon(metadados.dieta.alpha, conf.level = 0.95)

# Criar uma matriz de asteriscos para significância
asterisks_alpha.diet <- ifelse(p_values_alpha.diet < 0.001, "***", 
                               ifelse(p_values_alpha.diet < 0.01, "**", 
                                      ifelse(p_values_alpha.diet < 0.05, "*", "")))
```




```{r}
# Aumentar o tamanho da figura
options(repr.plot.width = 60, repr.plot.height = 80)

# Definir a altura e largura das células
cell_height <- 10  # Ajuste este valor conforme a necessidade
cell_width <- 15   # Ajuste este valor conforme a necessidade

# Gerar o heatmap com pheatmap, ajustando o tamanho das células e removendo a legenda do título
pheatmap(cor_matrix_alpha.diet, 
         display_numbers = asterisks_alpha.diet, # Mostra os asteriscos no mapa
         color = colorRampPalette(c("blue", "white", "red"))(100), # Paleta de cores ajustada
         breaks = seq(-1, 1, length.out = 101),  # Ajusta os breaks para a paleta de cores
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         cellheight = cell_height,  # Ajuste da altura das células
         cellwidth = cell_width,    # Ajuste da largura das células
         main = "Correlação entre Alpha Diversity e Dietas")




```

```{r}

colnames(saude_dieta)

# Selecionar as colunas numéricas de interesse e renomeá-las para metadados_shannon_selected
metadados.saude.alpha <- saude_dieta %>%
  select(shannon_entropy, pielou_evenness, faith_pd, observed_features, IL17A, IFNGamma, TNF, IL10, IL6, IL4, IL2, Age, Systolic, Diastolic, Weight, Height, Waist, Hip, HbA1c, COLESTEROL, LDL, HDL, VLDL, TRIGLICERIDES, TGO, TGP, GGT, GLICOSE, INSULINA, HOMA, PCR )

```




```{r}
# Função para calcular os p-valores das correlações
cor.mtest_shannon <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(p.mat)
}

```

```{r}
# Calcular a correlação de Pearson entre as variáveis selecionadas
cor_matrix_alpha.saude <- cor(metadados.saude.alpha, method = "spearman", use = "pairwise.complete.obs")

# Calcular os p-valores
p_values_alpha.saude <- cor.mtest_shannon(metadados.saude.alpha, conf.level = 0.95)

# Criar uma matriz de asteriscos para significância
asterisks_alpha.saude <- ifelse(p_values_alpha.diet < 0.001, "***", 
                                ifelse(p_values_alpha.diet < 0.01, "**", 
                                       ifelse(p_values_alpha.diet < 0.05, "*", "")))
```




```{r}
# Aumentar o tamanho da figura
options(repr.plot.width = 60, repr.plot.height = 80)

# Definir a altura e largura das células
cell_height <- 10  # Ajuste este valor conforme a necessidade
cell_width <- 15   # Ajuste este valor conforme a necessidade

# Gerar o heatmap com pheatmap, ajustando o tamanho das células e removendo a legenda do título
pheatmap(cor_matrix_alpha.saude, 
         display_numbers = asterisks_alpha.saude, # Mostra os asteriscos no mapa
         color = colorRampPalette(c("blue", "white", "red"))(100), # Paleta de cores ajustada
         breaks = seq(-1, 1, length.out = 101),  # Ajusta os breaks para a paleta de cores
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         cellheight = cell_height,  # Ajuste da altura das células
         cellwidth = cell_width,    # Ajuste da largura das células
         main = "Correlação entre Alpha Diversity e Saúde")




```
```{r}

colnames(saude_dieta)
#pielou_evenness x triglicerides

# usei esse mesmo código para o pelou, Faith, observed, shannon

ggplot(saude_dieta, aes(x = TRIGLICERIDES, y = shannon_entropy)) +
  
  geom_point() + 
  
  labs(x = "TRIGLICERIDES", y = "Shannon") + 
  
  theme_minimal() + 
  
  geom_smooth(method = "lm") +
  
  stat_poly_eq(aes(label = paste(..p.value.label.., sep = "~~~")),
               
               formula = y ~ x, 
               
               parse = TRUE) 


```

```{r}

colnames(saude_dieta)
#pielou_evenness x triglicerides

# usei esse mesmo código para o pelou, Faith, observed, shannon

ggplot(saude_dieta, aes(x = TRIGLICERIDES, y = faith_pd)) +
  
  geom_point() + 
  
  labs(x = "TRIGLICERIDES", y = "Faith") + 
  
  theme_minimal() + 
  
  geom_smooth(method = "lm") +
  
  stat_poly_eq(aes(label = paste(..p.value.label.., sep = "~~~")),
               
               formula = y ~ x, 
               
               parse = TRUE) 

```
```{r}

colnames(saude_dieta)
#pielou_evenness x triglicerides

# usei esse mesmo código para o pelou, Faith, observed, shannon

ggplot(saude_dieta, aes(x = TRIGLICERIDES, y = observed_features)) +
  
  geom_point() + 
  
  labs(x = "TRIGLICERIDES", y = "Richness") + 
  
  theme_minimal() + 
  
  geom_smooth(method = "lm") +
  
  stat_poly_eq(aes(label = paste(..p.value.label.., sep = "~~~")),
               
               formula = y ~ x, 
               
               parse = TRUE) 

```


```{r}

colnames(saude_dieta)
#pielou_evenness x triglicerides

# usei esse mesmo código para o pelou, Faith, observed, shannon

ggplot(saude_dieta, aes(x = LDL, y = shannon_entropy)) +
  
  geom_point() + 
  
  labs(x = "LDL", y = "Shannon") + 
  
  theme_minimal() + 
  
  geom_smooth(method = "lm") +
  
  stat_poly_eq(aes(label = paste(..p.value.label.., sep = "~~~")),
               
               formula = y ~ x, 
               
               parse = TRUE) 


```


```{r}

colnames(saude_dieta)
#pielou_evenness x triglicerides

# usei esse mesmo código para o pelou, Faith, observed, shannon

ggplot(saude_dieta, aes(x = COLESTEROL, y = shannon_entropy)) +
  
  geom_point() + 
  
  labs(x = "LDL", y = "Shannon") + 
  
  theme_minimal() + 
  
  geom_smooth(method = "lm") +
  
  stat_poly_eq(aes(label = paste(..p.value.label.., sep = "~~~")),
               
               formula = y ~ x, 
               
               parse = TRUE) 


```
```{r}
saude_dieta$RazaoWH <- saude_dieta$Waist / saude_dieta$Hip
# Renomeando a coluna RazaoWH para W2HRatio
colnames(saude_dieta)[colnames(saude_dieta) == "RazaoWH"] <- "W2HRatio"

```

```{r}

colnames(saude_dieta)
#pielou_evenness x triglicerides

# usei esse mesmo código para o pelou, Faith, observed, shannon

ggplot(saude_dieta, aes(x = Hip, y = shannon_entropy)) +
  
  geom_point() + 
  
  labs(x = "Hip", y = "Shannon") + 
  
  theme_minimal() + 
  
  geom_smooth(method = "lm") +
  
  stat_poly_eq(aes(label = paste(..p.value.label.., sep = "~~~")),
               
               formula = y ~ x, 
               
               parse = TRUE) 


```
```{r}

colnames(saude_dieta)
#pielou_evenness x triglicerides

# usei esse mesmo código para o pelou, Faith, observed, shannon

ggplot(saude_dieta, aes(x = Waist, y = shannon_entropy)) +
  
  geom_point() + 
  
  labs(x = "Waist", y = "Shannon") + 
  
  theme_minimal() + 
  
  geom_smooth(method = "lm") +
  
  stat_poly_eq(aes(label = paste(..p.value.label.., sep = "~~~")),
               
               formula = y ~ x, 
               
               parse = TRUE) 


```

```{r}

colnames(saude_dieta)
#pielou_evenness x triglicerides

# usei esse mesmo código para o pelou, Faith, observed, shannon

ggplot(saude_dieta, aes(x = W2HRatio, y = shannon_entropy)) +
  
  geom_point() + 
  
  labs(x = "Waist to Hip Ratio", y = "Shannon") + 
  
  theme_minimal() + 
  
  geom_smooth(method = "lm") +
  
  stat_poly_eq(aes(label = paste(..p.value.label.., sep = "~~~")),
               
               formula = y ~ x, 
               
               parse = TRUE) 


```

```{r}

colnames(saude_dieta)
#pielou_evenness x triglicerides

# usei esse mesmo código para o pelou, Faith, observed, shannon

ggplot(saude_dieta, aes(x = residual_magnesio_mg, y = shannon_entropy)) +
  
  geom_point() + 
  
  labs(x = "Magnesio", y = "Shannon") + 
  
  theme_minimal() + 
  
  geom_smooth(method = "lm") +
  
  stat_poly_eq(aes(label = paste(..p.value.label.., sep = "~~~")),
               
               formula = y ~ x, 
               
               parse = TRUE) 


```


```{r}

colnames(saude_dieta)
#pielou_evenness x triglicerides

# usei esse mesmo código para o pelou, Faith, observed, shannon

ggplot(saude_dieta, aes(x = residual_acidos_graxos_trans_g, y = shannon_entropy)) +
  
  geom_point() + 
  
  labs(x = "Ácidos Graxos Trans", y = "Shannon") + 
  
  theme_minimal() + 
  
  geom_smooth(method = "lm") +
  
  stat_poly_eq(aes(label = paste(..p.value.label.., sep = "~~~")),
               
               formula = y ~ x, 
               
               parse = TRUE) 



colnames(saude_dieta)
#pielou_evenness x triglicerides

# usei esse mesmo código para o pelou, Faith, observed, shannon

ggplot(saude_dieta, aes(x = residual_colesterol_mg , y = shannon_entropy)) +
  
  geom_point() + 
  
  labs(x = "Colesterol", y = "Shannon") + 
  
  theme_minimal() + 
  
  geom_smooth(method = "lm") +
  
  stat_poly_eq(aes(label = paste(..p.value.label.., sep = "~~~")),
               
               formula = y ~ x, 
               
               parse = TRUE) 


```