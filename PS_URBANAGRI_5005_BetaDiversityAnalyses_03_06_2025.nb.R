#======================================#
#   title: "Beta Diversity Analisys"
#     N=105
#======================================#

#Beta diversidade com Health Status

#=========== Carregar as bibliotecas necessárias=================#

library(ggplot2)
library(reshape2)
library(dplyr)
library(ggcorrplot)
library(pheatmap)
library(grid)
library(patchwork)
library(ggpubr)
library(ggpmisc)
library(phyloseq)
library(microbiomeMarker)
library(tidyverse)
library(vegan)
library(DESeq2)
library(data.table)
library(qiime2R)
library(viridis)
library(ape)




#============ Importar o arquivo CSV ==============================#

#metadados.metabolic.score <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados.metabolic.score.csv", stringsAsFactors = FALSE)

#============ Weighted Unifrac =================#


# Realizar PCoA na matriz de distância weighted_unifrac


# Executa a PCoA
wei_unifrac.pcoa <- cmdscale(weighted.unifrac, eig = TRUE, k = 3)

# Verifica a saída
print(wei_unifrac.pcoa)


#============================= Unweighted Unifrac =================#

# Realizar PCoA na matriz de distância unweighted_unifrac

unwei_unifrac.pcoa <- cmdscale(unweighted.unifrac, eig = TRUE, k = 3)  # k = número de dimensões

ids_unifrac <- colnames(unweighted.unifrac)  # ou rownames()
ids_metadata <- metadados.metabolic.score$Sample.id  # substitua pela coluna correta se o nome for diferente

# Interseção dos dois conjuntos
ids <- intersect(ids_unifrac, ids_metadata)

length(ids)  # só pra conferir quantos ficaram

# Certifique-se de que os dados estejam na mesma ordem para ambos


permanova_unweighted <- adonis2(
  unweighted.unifrac[ids, ids] ~ metadados.metabolic.score[match(ids, metadados.metabolic.score$Sample.id), "Health_Status"],
  permutations = 999
)

permanova_unweighted



# Criar texto do resultado
res_unweighted <- paste0("PERMANOVA: R² = ", round(permanova_unweighted$R2[1], 3),
                         ", p = ", permanova_unweighted$`Pr(>F)`[1])




# Rodar PERMANOVA com os dados alinhados
permanova_weighted <- adonis2(weighted.unifrac[ids, ids] ~ metadados.metabolic.score[ids, "Health_Status"], permutations = 999)

# Criar texto do resultado
res_weighted <- paste0("PERMANOVA: R² = ", round(permanova_weighted$R2[1], 3),
                       ", p = ", permanova_weighted$`Pr(>F)`[1])



#=====> Weighted

# Realizar PCoA na matriz de distância unweighted_unifrac

wei_unifrac.pcoa <- cmdscale(weighted.unifrac, eig = TRUE, k = 3)  # k = número de dimensões

ids_unifrac <- colnames(weighted.unifrac)  # ou rownames()
ids_metadata <- metadados.metabolic.score$Sample.id  # substitua pela coluna correta se o nome for diferente

# Interseção dos dois conjuntos
ids <- intersect(ids_unifrac, ids_metadata)

length(ids)  # só pra conferir quantos ficaram

# Certifique-se de que os dados estejam na mesma ordem para ambos


permanova_weighted <- adonis2(
  weighted.unifrac[ids, ids] ~ metadados.metabolic.score[match(ids, metadados.metabolic.score$Sample.id), "Health_Status"],
  permutations = 999
)

permanova_weighted



# Criar texto do resultado
res_weighted <- paste0("PERMANOVA: R² = ", round(permanova_weighted$R2[1], 3),
                         ", p = ", permanova_weighted$`Pr(>F)`[1])




# Rodar PERMANOVA com os dados alinhados
permanova_weighted <- adonis2(weighted.unifrac[ids, ids] ~ metadados.metabolic.score[ids, "Health_Status"], permutations = 999)

# Criar texto do resultado
res_weighted <- paste0("PERMANOVA: R² = ", round(permanova_weighted$R2[1], 3),
                       ", p = ", permanova_weighted$`Pr(>F)`[1])



#====== Plot Weighted ===========#
#mudar os nomes das regioes para ingles
#metadados.all.filtrado$Region <- recode(metadados.all.filtrado$Region,
#                                           "Norte" = "North",
#                                          "Sul" = "South",
#                                         "Leste" = "East",
#                                        "Sul Adjacente" = "South Adjacent",
#                                       "Leste Adjacente" = "East Adjacent"
#)




### OBS: USAR ESSE PLOT!!!!


ggplot(merge(wei_unifrac.pcoa$points, metadados.metabolic.score, by.x = "row.names", by.y = "Sample.id")) + 
  geom_point(aes(x = V1, y = V2, color = Health_Status), size = 3) + 
  scale_color_viridis_d(option = "C", name = "Health Status") +
  labs(title = paste("Weighted UniFrac — Health Status         ", res_weighted),
       x = "PCoA1", y = "PCoA2") +
  theme_minimal()

ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/weighted_unifrac_pcoa_Health_Status.png", dpi = 600, width = 8, height = 6, units = "in")



#====== Plot UnWeighted ===========#

### OBS: Usar esse Plot

ggplot(merge(unwei_unifrac.pcoa$points, metadados.metabolic.score, by.x = "row.names", by.y = "Sample.id")) + 
  geom_point(aes(x = V1, y = V2, color = Health_Status), size = 3) + 
  scale_color_viridis_d(option = "C", name = "Health Status") +
  labs(title = paste("Unweighted UniFrac — Health Status         ", res_unweighted),
       x = "PCoA1", y = "PCoA2") +
  theme_minimal()

ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/unweighted_unifrac_pco_health_statusa.png", dpi = 600, width = 8, height = 6, units = "in")

#==============


library(vegan)


# Loop por variável
for (v in vars) {
  cat("\n### PERMANOVA for:", v, "###\n")
  
  # Remover NAs só da variável e manter IDs que estão na matriz
  ids <- rownames(metadados.all.filtrado[!is.na(metadados.all.filtrado[[v]]), ])
  ids <- intersect(ids, rownames(weighted.unifrac))
  
  # Rodar PERMANOVA direto
  result_weighted_permanova <- adonis2(weighted.unifrac[ids, ids] ~ metadados.all.filtrado[ids, v], permutations = 999)
  print(result_weighted_permanova)
}



#######
permanova_weighted <- adonis2(weighted.unifrac[ids, ids] ~ metadados.all.filtrado[ids, "Region"], permutations = 999)



#######============ Permanova Weighted =================#
library(vegan)

vars <- c("Region", "IL17A", "IFNGamma", "TNF", "IL10", "IL6", "IL4", "IL2", "Age", "Sex", "Raca", 
          "Systolic", "Diastolic", "IMC", "W.H", 
          "UREIA", "CREATININA", "HbA1c", "COLESTEROL", "LDL", "HDL", "VLDL", "TRIGLICERIDES", 
          "TGO", "TGP", "GGT", "GLICOSE", "INSULINA", "HOMA.IR", "PCR", 
          "carboidrato_total_g", "proteina_g", "lipidios_g", "fibra_alimentar_g", "colesterol_mg", 
          "acidos_graxos_saturados_g", "acidos_graxos_monoinsaturados_g", "acidos_graxos_poliinsaturados_g", 
          "acidos_graxos_trans_g", "calcio_mg", "ferro_mg", "sodio_mg", "magnesio_mg", "fosforo_mg", 
          "potassio_mg", "manganes_mg", "zinco_mg", "cobre_mg", "selenio_mcg", "vitamina_A_RAE_mcg", 
          "vitamina_D_mcg", "vitamina_E_mg", "tiamina_mg", "riboflavina_mg", "niacina_mg", 
          "vitamina_B6_mg", "vitamina_B12_mcg", "vitamina_C_mg", "equivalente_de_folato_mcg", 
          "sal_de_adicao_g", "acucar_de_adicao_g", 
          "BHEI_R_Score_Total", "Percentual_NOVA_group_1", "Percentual_NOVA_group_2", 
          "Percentual_NOVA_group_3", "BMI", "TyG", "VAI")

# Criar lista para armazenar resultados
result_list <- list()

# Fixar aleatoriedade ANTES do loop
set.seed(123)  

# Loop PERMANOVA
for (v in vars) {
  df <- metadados.all.filtrado[!is.na(metadados.all.filtrado[[v]]), ]
  ids <- intersect(rownames(df), rownames(weighted.unifrac))
  df <- df[ids, , drop = FALSE]
  
  if (length(ids) < 3) next
  
  res <- tryCatch({
    adonis2(weighted.unifrac[ids, ids] ~ df[[v]], permutations = 999)
  }, error = function(e) return(NULL))
  
  if (!is.null(res)) {
    result_list[[v]] <- data.frame(
      Variável = v,
      F_value = res$F[1],
      R2 = res$R2[1],
      p_value = res$`Pr(>F)`[1]
    )
  }
}

# Juntar tudo em um data.frame
permanova_results <- do.call(rbind, result_list)

# Visualizar
print(permanova_results)

#Filtrar apenas os significativos (p < 0.05)
permanova_sig <- permanova_results[permanova_results$p_value < 0.05, ]
print(permanova_sig)

#Ordenar por p-valor
permanova_sorted <- permanova_results[order(permanova_results$p_value), ]
head(permanova_sorted, 10)  # top 10


#library(ggplot2)

ggplot(permanova_sig, aes(x = reorder(Variável, -R2), y = R2)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "PERMANOVA WEIGHTED – R² of significant variables",
    x = "Variable",
    y = expression(R^2)
  )

ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/permanova_significant_variables.png",
       width = 20, height = 15, units = "cm", dpi = 300)


#Volcano plot

# Suponha que seu data frame se chame permanova_results
# Ele precisa ter colunas: Variável, R2 e p_value

permanova_results$log10_p <- -log10(permanova_results$p_value)
permanova_results$significativo <- permanova_results$p_value < 0.05

ggplot(permanova_results, aes(x = R2, y = log10_p, label = Variável)) +
  geom_point(aes(color = significativo), size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text(data = subset(permanova_results, significativo == TRUE),
            hjust = 0.5, vjust = -0.5, size = 3) +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "steelblue"),
                     labels = c("Not significant", "Significant")) +
  theme_minimal() +
  labs(
    title = "PERMANOVA WEIGHTED – Effect size vs. Significance",
    x = expression(R^2),
    y = expression(-log[10](p~value)),
    color = "Significance"
  ) +
  theme(legend.position = "bottom")


#Com nome de todas as variaveis
# Instale ggrepel se ainda não tiver
# install.packages("ggrepel")
library(ggrepel)

ggplot(permanova_results, aes(x = R2, y = -log10(p_value), label = Variável)) +
  geom_point(aes(color = significativo), size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(size = 3, max.overlaps = Inf) +  # Adiciona o nome de todas as variáveis
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "steelblue"),
                     labels = c("Not significant", "Significant")) +
  theme_minimal() +
  labs(
    title = "PERMANOVA WEIGHTED – Effect size vs. Significance",
    x = expression(R^2),
    y = expression(-log[10](p~value)),
    color = "Significance"
  ) +
  theme(legend.position = "bottom")



ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/volcanoPlot_significant_variables_withLabel.png",
       width = 20, height = 15, units = "cm", dpi = 300)

####================== Permanova Unweighted ====#
# Criar lista para armazenar resultados da PERMANOVA com unweighted UniFrac
data_list_unw <- list()

# Fixar aleatoriedade para reprodutibilidade
tidyverse::set.seed(123)

# Loop PERMANOVA com matriz unweighted.unifrac
for (v in vars) {
  df_unw <- metadados.all.filtrado[!is.na(metadados.all.filtrado[[v]]), ]
  ids_unw <- intersect(rownames(df_unw), rownames(unweighted.unifrac))
  df_unw <- df_unw[ids_unw, , drop = FALSE]
  
  if (length(ids_unw) < 3) next
  
  res <- tryCatch({
    adonis2(unweighted.unifrac[ids_unw, ids_unw] ~ df_unw[[v]], permutations = 999)
  }, error = function(e) return(NULL))
  
  if (!is.null(res)) {
    data_list_unw[[v]] <- data.frame(
      Variable = v,
      F_value = res$F[1],
      R2 = res$R2[1],
      p_value = res$`Pr(>F)`[1]
    )
  }
}

# Juntar resultados em um data frame
permanova_results_unw <- do.call(rbind, data_list_unw)

# Filtrar variáveis significativas (p < 0.05)
permanova_sig_unw <- permanova_results_unw[permanova_results_unw$p_value < 0.05, ]

# Ordenar por p-valor crescente
permanova_sorted_unw <- permanova_results_unw[order(permanova_results_unw$p_value), ]

# Barplot das variáveis significativas
library(ggplot2)
library(dplyr)



# Dicionário de tradução
variaveis_ingles <- c(
  "acidos_graxos_trans_g" = "Trans fatty acids",
  "Percentual_NOVA_group_1" = "NOVA group 1 %",
  "acucar_de_adicao_g" = "Added sugar",
  "HbA1c" = "HbA1c",
  "IL4" = "IL-4",
  "TGP" = "TGP",
  "IMC" = "BMI",
  "selenio_mcg" = "Selenium",
  "TyG" = "TyG index",
  "BMI" = "BMI"
)

# Aplica a tradução
permanova_sig_unw <- permanova_sig_unw %>%
  mutate(Variable = recode(Variable, !!!variaveis_ingles))


# Gráfico com nomes em inglês
ggplot(permanova_sig_unw, aes(x = reorder(Variable, -R2), y = R2)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "PERMANOVA UNWEIGHTED – R² of Significant Variables",
    x = "Variable",
    y = expression(R^2)
  )

# Salvar em alta resolução
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/permanova_unweighted_significant_variables.png",
       width = 20, height = 15, units = "cm", dpi = 300)





# Volcano plot para todos os resultados
permanova_results_unw$log10_p <- -log10(permanova_results_unw$p_value)
permanova_results_unw$significant <- permanova_results_unw$p_value < 0.05

# Aplica a tradução
permanova_results_unw <- permanova_results_unw %>%
  mutate(Variable = recode(Variable, !!!variaveis_ingles))


ggplot(permanova_results_unw, aes(x = R2, y = log10_p, label = Variable)) +
  geom_point(aes(color = significant), size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text(data = subset(permanova_results_unw, significant == TRUE),
            hjust = 0.5, vjust = -0.5, size = 3) +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "steelblue"),
                     labels = c("Not significant", "Significant")) +
  theme_minimal() +
  labs(
    title = "PERMANOVA UNWEIGHTED – Effect size vs. Significance",
    x = expression(R^2),
    y = expression(-log[10](p~value)),
    color = "Significance"
  ) +
  theme(legend.position = "bottom")

#Com nome de todas as variaveis
# Instale ggrepel se ainda não tiver
# install.packages("ggrepel")
library(ggrepel)

ggplot(permanova_results_unw, aes(x = R2, y = -log10(p_value), label = Variable)) +
  geom_point(aes(color = significant), size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(size = 3, max.overlaps = Inf) +  # ← Nome em todas as variáveis
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "steelblue"),
                     labels = c("Not significant", "Significant")) +
  theme_minimal() +
  labs(
    title = "PERMANOVA UNWEIGHTED – Effect size vs. Significance",
    x = expression(R^2),
    y = expression(-log[10](p~value)),
    color = "Significance"
  ) +
  theme(legend.position = "bottom")


# Salvar volcano plot
# Salvar em alta resolução
ggsave("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/volcanoPlot_unweighted_significant_variables_withlabel.png",
       width = 20, height = 15, units = "cm", dpi = 300)


# Opcional: salvar em CSV
write.csv(permanova_results, "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/permanova_results_weighted.csv", row.names = FALSE)

##########====== Jaccard e Bray Curtis(está com N109) =========#########


library(ape)

pcoa_jaccard <- pcoa(jaccard)

# Coordenadas das amostras nos dois primeiros eixos
pcoa_points_jaccard <- as.data.frame(pcoa_jaccard$vectors[, 1:3])

expl_var_jaccard <- round(100 * pcoa_jaccard$values$Relative_eig[1:3], 1)
# Ex: PCoA1 = expl_var_jaccard[1]%, PCoA2 = expl_var_jaccard[2]%



pcoa_bray <- pcoa(bray.curtis)

# Coordenadas das amostras nos dois primeiros eixos
pcoa_points_bray <- as.data.frame(pcoa_bray$vectors[, 1:3])
# Pronto: sem coluna SampleID, só rownames com IDs

expl_var_bray <- round(100 * pcoa_bray$values$Relative_eig[1:3], 1)
expl_var_bray




ggplot(pcoa_points_jaccard, aes(x = Axis.1, y = Axis.2)) +
  geom_point(size = 2) +
  labs(title = "PCoA — Jaccard",
       x = paste0("PCoA1 (", expl_var_jaccard[1], "%)"),
       y = paste0("PCoA2 (", expl_var_jaccard[2], "%)")) +
  theme_minimal()

library(viridis)

# Obter a cor mais escura da paleta viridis "C"
cor_escura <- viridis(1, option = "C")

# Gráfico com cor fixa
ggplot(pcoa_points_jaccard, aes(x = Axis.1, y = Axis.2)) +
  geom_point(color = cor_escura, size = 3) +
  labs(title = "PCoA — Jaccard",
       x = paste0("PCoA1 (", expl_var_jaccard[1], "%)"),
       y = paste0("PCoA2 (", expl_var_jaccard[2], "%)")) +
  theme_minimal()



library(viridis)

# Obter a cor mais escura da paleta viridis "C"
cor_escura <- viridis(1, option = "C")

# Gráfico Bray-Curtis
ggplot(pcoa_points_bray, aes(x = Axis.1, y = Axis.2)) +
  geom_point(color = cor_escura, size = 3) +
  labs(title = "PCoA — Bray-Curtis",
       x = paste0("PCoA1 (", expl_var_bray[1], "%)"),
       y = paste0("PCoA2 (", expl_var_bray[2], "%)")) +
  theme_minimal()


pcoa_points_jaccard$ID <- rownames(pcoa_points_jaccard)
metadata$ID <- rownames(metadata)

merged_jaccard <- merge(pcoa_points_jaccard, metadata, by = "ID")

ggplot(merged_jaccard, aes(x = Axis.1, y = Axis.2, color = Region)) +
  geom_point(size = 3) +
  scale_color_viridis_d(option = "C", name = "Region") +
  labs(title = "PCoA — Jaccard",
       x = paste0("PCoA1 (", expl_var_jaccard[1], "%)"),
       y = paste0("PCoA2 (", expl_var_jaccard[2], "%)")) +
  theme_minimal()


pcoa_points_bray$ID <- rownames(pcoa_points_bray)

merged_bray <- merge(pcoa_points_bray, metadata, by = "ID")

ggplot(merged_bray, aes(x = Axis.1, y = Axis.2, color = Region)) +
  geom_point(size = 3) +
  scale_color_viridis_d(option = "C", name = "Region") +
  labs(title = "PCoA — Bray-Curtis",
       x = paste0("PCoA1 (", expl_var_bray[1], "%)"),
       y = paste0("PCoA2 (", expl_var_bray[2], "%)")) +
  theme_minimal()


#======== Statistical Power Calculation: TyG ===========#

library(vegan)

#primeiro ver quais rows nao estao nos dois objetos
# Checar os IDs
length(rownames(unweighted.unifrac))  # deve ser 105
length(dados_tyG_bin$Sample.id)       # deve ser 97

# Pegar os IDs presentes nos dois
ids_comuns <- intersect(rownames(unweighted.unifrac), dados_tyG_bin$Sample.id)

# Filtrar a matriz de distância
unifrac_filtrado <- unweighted.unifrac[ids_comuns, ids_comuns]

# Filtrar os metadados
dados_tyG_bin_filtrado <- subset(dados_tyG_bin, Sample.id %in% ids_comuns)

all(rownames(unifrac_filtrado) == dados_tyG_bin_filtrado$Sample.id)  # deve ser TRUE

# Reordenar os metadados para alinhar com a matriz
dados_tyG_bin_filtrado <- dados_tyG_bin_filtrado[match(rownames(unifrac_filtrado), dados_tyG_bin_filtrado$Sample.id), ]


# Rodar PERMANOVA com 999 permutações
permanova_result <- adonis2(unifrac_filtrado ~ TyG_group2, data = dados_tyG_bin_filtrado, permutations = 999)

# Ver R²
r2_value <- permanova_result$R2[1]  # pega o R² do fator TyG_group2
print(r2_value)


# Função para simular PERMANOVA várias vezes e estimar o poder
simula_power_permanova <- function(dist_matrix, metadados.bioquimicos.alpha, group_var, r2_obs = NULL, n_perm = 999, n_sim = 1000, alpha = 0.05) {
  p_vals <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    # Copia os metadados e embaralha a variável de grupo
    metadados.bioquimicos.alpha$grupo_embaralhado <- sample(metadados.bioquimicos.alpha[[group_var]])
    
    # Rodar PERMANOVA com a variável embaralhada
    p <- adonis2(dist_matrix ~ grupo_embaralhado, data = metadados.bioquimicos.alpha, permutations = n_perm)$`Pr(>F)`[1]
    
    # Armazena o p-valor
    p_vals[i] <- p
  }
  
  # Calcula o poder estatístico: proporção de p-valores < alpha
  power <- mean(p_vals < alpha)
  return(power)
}

library(pwr)


power_result <- simula_power_permanova(unifrac_filtrado, dados_tyG_bin_filtrado, "TyG_group2")
print(paste("Estimated power =", round(power_result * 100, 1), "%"))


#===== Apresentação do grafico com numeros =======#

# 1. Calcular a PCoA
pcoa_result <- cmdscale(unifrac_filtrado, k = 2, eig = TRUE)
pcoa_df <- as.data.frame(pcoa_result$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$Sample.id <- rownames(pcoa_df)

# 2. Juntar com os metadados filtrados
pcoa_df <- merge(pcoa_df, dados_tyG_bin_filtrado[, c("Sample.id", "TyG_group2")], by = "Sample.id")

# 3. Calcular % de variância explicada pelos eixos
eig_vals <- pcoa_result$eig
var_exp <- round(100 * eig_vals[1:2] / sum(eig_vals), 1)

# 4. Inserir valores do PERMANOVA
p_val <- permanova_result$`Pr(>F)`[1]
r2_val <- permanova_result$R2[1]
power_val <- power_result

# 5. Plotar com ggplot2
library(ggplot2)

ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = TyG_group2)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCoA - Unweighted UniFrac by TyG Group",
    subtitle = paste0("PERMANOVA: p = ", format(p_val, digits = 3),
                      " | R² = ", round(r2_val, 3),
                      " | Power = ", round(power_val * 100, 1), "%"),
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    color = "TyG Group"
  ) +
  scale_color_manual(values = c("Baixo" = "#1f77b4", "Alto" = "#ff7f0e"))


# Salvar o gráfico como PNG em alta resolução
ggsave(
  filename = "PCoA_TyG_UnweightedUniFrac.png",   # nome do arquivo
  path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação",      # pasta onde salvar
  dpi = 300,                                     # resolução
  width = 8, height = 6,                         # em polegadas
  units = "in"                                   # unidade de medida
)



#======== Statistical Power Calculation: Region x Unweighted ===========#

library(vegan)

#primeiro ver quais rows nao estao nos dois objetos
# Checar os IDs
length(rownames(unweighted.unifrac))  # deve ser 106
length(metadados.all.filtrado$Sample.id)       # deve ser 105

# Pegar os IDs presentes nos dois
ids_comuns <- intersect(rownames(unweighted.unifrac), metadados.all.filtrado$Sample.id)

# Filtrar a matriz de distância
unweighted.unifrac <- unweighted.unifrac[ids_comuns, ids_comuns]

# Filtrar os metadados
metadados.all.filtrado <- subset(metadados.all.filtrado, Sample.id %in% ids_comuns)

all(rownames(unifrac_filtrado) == dados_tyG_bin_filtrado$Sample.id)  # deve ser TRUE

# Reordenar os metadados para alinhar com a matriz
metadados.all.filtrado <- metadados.all.filtrado[match(rownames(unweighted.unifrac), metadados.all.filtrado$Sample.id), ]


# Rodar PERMANOVA com 999 permutações
permanova_result <- adonis2(unweighted.unifrac ~ Region, data = metadados.all.filtrado, permutations = 999)

# Ver R²
r2_value <- permanova_result$R2[1]  # pega o R² do fator TyG_group2
print(r2_value)


# Função segura para simular PERMANOVA e estimar o poder
simula_power_permanova <- function(dist_matrix, metadata, group_var, r2_obs = NULL, n_perm = 999, n_sim = 1000, alpha = 0.05) {
  p_vals <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    # Faz uma cópia dos metadados e embaralha a variável de grupo
    temp_metadata <- metadata
    temp_metadata$grupo_embaralhado <- sample(temp_metadata[[group_var]])
    
    # Rodar PERMANOVA com a variável embaralhada
    p <- adonis2(dist_matrix ~ grupo_embaralhado, data = temp_metadata, permutations = n_perm)$`Pr(>F)`[1]
    
    # Armazena o p-valor
    p_vals[i] <- p
  }
  
  # Calcula o poder estatístico: proporção de p-valores < alpha
  power <- mean(p_vals < alpha)
  return(power)
}



library(pwr)

power_result <- simula_power_permanova(unweighted.unifrac, metadados.all.filtrado, "Region")
print(paste("Estimated power =", round(power_result * 100, 1), "%"))


#===== Apresentação do grafico com numeros =======#

# 1. Calcular a PCoA
pcoa_result <- cmdscale(unweighted.unifrac, k = 2, eig = TRUE)
pcoa_df <- as.data.frame(pcoa_result$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$Sample.id <- rownames(pcoa_df)




ggplot(pcoa_df) + 
  geom_point(aes(x = PCoA1, y = PCoA2), size = 3, color = "steelblue") + 
  labs(
    title = "PCoA — UnWeighted UniFrac",
    x = "PCoA1", 
    y = "PCoA2"
  ) +
  theme_minimal(base_size = 14)

# Salvar o gráfico como PNG em alta resolução
ggsave(
  filename = "PCoA_UnWeightedUniFrac.png",   # nome do arquivo
  path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação",      # pasta onde salvar
  dpi = 300,                                     # resolução
  width = 8, height = 6,                         # em polegadas
  units = "in"                                   # unidade de medida
)




# 2. Juntar com os metadados filtrados
pcoa_df <- merge(pcoa_df, dados_tyG_bin_filtrado[, c("Sample.id", "TyG_group2")], by = "Sample.id")

# 3. Calcular % de variância explicada pelos eixos
eig_vals <- pcoa_result$eig
var_exp <- round(100 * eig_vals[1:2] / sum(eig_vals), 1)

# 4. Inserir valores do PERMANOVA
p_val <- permanova_result$`Pr(>F)`[1]
r2_val <- permanova_result$R2[1]
power_val <- power_result

# 5. Plotar com ggplot2
library(ggplot2)

ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = TyG_group2)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCoA - Unweighted UniFrac by TyG Group",
    subtitle = paste0("PERMANOVA: p = ", format(p_val, digits = 3),
                      " | R² = ", round(r2_val, 3),
                      " | Power = ", round(power_val * 100, 1), "%"),
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    color = "TyG Group"
  ) +
  scale_color_manual(values = c("Baixo" = "#1f77b4", "Alto" = "#ff7f0e"))


# Salvar o gráfico como PNG em alta resolução
ggsave(
  filename = "PCoA_TyG_UnweightedUniFrac.png",   # nome do arquivo
  path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação",      # pasta onde salvar
  dpi = 300,                                     # resolução
  width = 8, height = 6,                         # em polegadas
  units = "in"                                   # unidade de medida
)


# Supondo que você já tenha:
# unwei_unifrac.pcoa → resultado da PCoA
# unweighted.unifrac → matriz de distância
# metadados.all.filtrado → metadados com coluna Region

# 1. Rodar PERMANOVA (caso não tenha feito)
permanova_result <- adonis2(unweighted.unifrac ~ Region, data = metadados.all.filtrado, permutations = 999)

# 2. Extrair valores para o gráfico
p_val <- format(permanova_result$`Pr(>F)`[1], digits = 3)
r2_val <- round(permanova_result$R2[1], 3)
power_val <- 4.6  # resultado da sua simulação anterior

# 3. Preparar dados da PCoA
pcoa_df <- merge(unwei_unifrac.pcoa$points, metadados.all.filtrado, 
                 by.x = "row.names", by.y = "Sample.id")

# 4. Criar gráfico com título personalizado
library(ggplot2)
library(viridis)

ggplot(pcoa_df) + 
  geom_point(aes(x = V1, y = V2, color = Region), size = 3) + 
  scale_color_viridis_d(option = "C", name = "Region") +
  labs(
    title = paste0(
      "PCoA — Unweighted UniFrac (Region)\n",
      "PERMANOVA: p = ", p_val,
      " | R² = ", r2_val,
      " | Power = ", power_val, "%"
    ),
    x = "PCoA1", 
    y = "PCoA2"
  ) +
  theme_minimal(base_size = 14)

# Salvar o gráfico como PNG em alta resolução
ggsave(
  filename = "PCoA_Region_UnweightedUniFrac.png",   # nome do arquivo
  path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação",      # pasta onde salvar
  dpi = 300,                                     # resolução
  width = 8, height = 6,                         # em polegadas
  units = "in"                                   # unidade de medida
)





#======== Statistical Power Calculation: Region x Weighted ===========#

library(vegan)


# Rodar PERMANOVA com 999 permutações
permanova_result <- adonis2(weighted.unifrac ~ Region, data = metadados.all.filtrado, permutations = 999)

# Ver R²
r2_value <- permanova_result$R2[1]  # pega o R² do fator TyG_group2
print(r2_value)


# Função segura para simular PERMANOVA e estimar o poder
simula_power_permanova <- function(dist_matrix, metadata, group_var, r2_obs = NULL, n_perm = 999, n_sim = 1000, alpha = 0.05) {
  p_vals <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    # Faz uma cópia dos metadados e embaralha a variável de grupo
    temp_metadata <- metadata
    temp_metadata$grupo_embaralhado <- sample(temp_metadata[[group_var]])
    
    # Rodar PERMANOVA com a variável embaralhada
    p <- adonis2(dist_matrix ~ grupo_embaralhado, data = temp_metadata, permutations = n_perm)$`Pr(>F)`[1]
    
    # Armazena o p-valor
    p_vals[i] <- p
  }
  
  # Calcula o poder estatístico: proporção de p-valores < alpha
  power <- mean(p_vals < alpha)
  return(power)
}



library(pwr)

power_result_weighted <- simula_power_permanova(weighted.unifrac, metadados.all.filtrado, "Region")
print(paste("Estimated power =", round(power_result_weighted * 100, 1), "%"))


#===== Apresentação do grafico com numeros =======#

# 1. Calcular a PCoA
pcoa_result <- cmdscale(weighted.unifrac, k = 2, eig = TRUE)
pcoa_df <- as.data.frame(pcoa_result$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$Sample.id <- rownames(pcoa_df)

# 2. Juntar com os metadados filtrados
pcoa_df <- merge(pcoa_df, metadados.all.filtrado[, c("Sample.id", "Region")], by = "Sample.id")

# 3. Calcular % de variância explicada pelos eixos
eig_vals <- pcoa_result$eig
var_exp <- round(100 * eig_vals[1:2] / sum(eig_vals), 1)

# 4. Inserir valores do PERMANOVA
p_val <- round(permanova_result$`Pr(>F)`[1], 2)
r2_val <- round(permanova_result$R2[1], 2)
power_val <- paste0(round(power_result_weighted * 100, 1), "%")


# 5. Plotar com ggplot2
library(ggplot2)




ggplot(pcoa_df) + 
  geom_point(aes(x = PCoA1, y = PCoA2, color = Region), size = 3) + 
  scale_color_viridis_d(option = "C", name = "Region") +
  labs(
    title = paste0(
      "PCoA — Weighted UniFrac (Region)\n",
      "PERMANOVA: p = ", p_val,
      " | R² = ", r2_val,
      " | Power = ", power_val
    ),
    x = "PCoA1", 
    y = "PCoA2"
  ) +
  theme_minimal(base_size = 14)

# Salvar o gráfico como PNG em alta resolução
ggsave(
  filename = "PCoA_Region_WeightedUniFrac.png",   # nome do arquivo
  path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação",      # pasta onde salvar
  dpi = 300,                                     # resolução
  width = 8, height = 6,                         # em polegadas
  units = "in"                                   # unidade de medida
)


# Supondo que você já tenha:
# unwei_unifrac.pcoa → resultado da PCoA
# unweighted.unifrac → matriz de distância
# metadados.all.filtrado → metadados com coluna Region

# 1. Rodar PERMANOVA (caso não tenha feito)
permanova_result <- adonis2(weighted.unifrac ~ Region, data = metadados.all.filtrado, permutations = 999)

# 2. Extrair valores para o gráfico
p_val <- format(permanova_result$`Pr(>F)`[1], digits = 3)
r2_val <- round(permanova_result$R2[1], 3)
power_val <- 5.1  # resultado da sua simulação anterior

# 3. Preparar dados da PCoA
pcoa_df <- merge(wei_unifrac.pcoa$points, metadados.all.filtrado, 
                 by.x = "row.names", by.y = "Sample.id")

# 4. Criar gráfico com título personalizado
library(ggplot2)
library(viridis)

ggplot(pcoa_df) + 
  geom_point(aes(x = V1, y = V2, color = Region), size = 3) + 
  scale_color_viridis_d(option = "C", name = "Region") +
  labs(
    title = paste0(
      "PCoA — Weighted UniFrac (Region)\n",
      "PERMANOVA: p = ", p_val,
      " | R² = ", r2_val,
      " | Power = ", power_val, "%"
    ),
    x = "PCoA1", 
    y = "PCoA2"
  ) +
  theme_minimal(base_size = 14)

# Salvar o gráfico como PNG em alta resolução
ggsave(
  filename = "PCoA_Region_WeightedUniFrac.png",   # nome do arquivo
  path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação",      # pasta onde salvar
  dpi = 300,                                     # resolução
  width = 8, height = 6,                         # em polegadas
  units = "in"                                   # unidade de medida
)

#===============

ggplot(pcoa_df) + 
  geom_point(aes(x = PCoA1, y = PCoA2), size = 3, color = "steelblue") + 
  labs(
    title = "PCoA — Weighted UniFrac",
    x = "PCoA1", 
    y = "PCoA2"
  ) +
  theme_minimal(base_size = 14)

# Salvar o gráfico como PNG em alta resolução
ggsave(
  filename = "PCoA_WeightedUniFrac.png",   # nome do arquivo
  path = "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação",      # pasta onde salvar
  dpi = 300,                                     # resolução
  width = 8, height = 6,                         # em polegadas
  units = "in"                                   # unidade de medida
)
