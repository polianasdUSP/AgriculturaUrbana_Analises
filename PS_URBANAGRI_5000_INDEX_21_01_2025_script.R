---
  #INDEX Urban Agriculture Analisys
  # January 2025
  
  ---
  
# getwd() # Mostra o diretório atual
# setwd("caminho/desejado") # Define um novo diretório

#### WORKING WITH DATA ####

# ====================================== #
# 1. Cleaning and organizing metadata
# ====================================== #

# This script imports raw data:
# Cleaning and organization of metadata: 
#(output object: metadados.raw, alpha.shannon, alpha.pd, alpha.evenness, alpha.richness, alpha.all, unweighted_unifrac, weighted_unifrac, metadados.shannon, metadados.alpha.all, metadados, unweighted_unifrac)
# outputs objects are saved in: "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri"


source ("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/AgriculturaUrbana_Analises/PS_5001_LimpezaDados_08_07_2024.nb.Rmd")


# ================================================= #
# 2. POPULATION DATA - Population caracterization
# ================================================= #


# This script is for the characterization of the population. Summarizing clinical traits, age, sex, BMI, Metabolic Syndrome:

#(output object: population_characteristics_table_, population_characteristics_table2, boxplot_idade, age_sex_summary, ParasitologicalTest)
#all output images are saved in: "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Tables_Images"

source("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/AgriculturaUrbana_Analises/PS_URBANAGRI_5003_Summary_Plots_DescriptiveAnalyses_09_07_2024.nb.Rmd")


# ======================================= #
# 3. POPULATION DATA - METABOLIC SYNDROME
# ======================================= #


# This script is for the characterization of the population regarding Metabolic Syndrome:

#(output object: Glycemia_HbA1c_Insulin, HbA1_glucose_Insulin_High Levels,  )

# output em "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Tables_Images"

source("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/AgriculturaUrbana_Analises/PS_URBANAGRI_5014_VENN_14_01_2025.Rmd")


# ======================================= #
# 4. Alpha Diversity N=105
# ======================================= #

# This script is for alpha diversity analyses. With inflamatory and metabolic parameters, with parasitology testes, with presence of treponema, with diet, NOVA classification


#Output objects em "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/..." : 

#  heatmap_cytokines_raw_pvalues.png
# heatmap_cytokines_FDR_adjusted.png
# heatmap_Clinical Biomarkers_FDRadjusted.png
#  heatmap_Clinical BiomarkersXalpha_p_raw.png
#  heatmap_dietXalpha_FDRadjusted.png
#  heatmap_dietXalpha_p_raw.png
#  heatmap_MetabolicMarkers1_FDRadjusted.png
#  heatmap_MetabolicMarkers1Xalpha_p_raw.png
# heatmap_NOVAXalpha_FDRadjusted.png
#  heatmap_NOVAXalpha_p_raw.png

#Scatterplots:
#  Scatterplot_TGP_vs_alpha_diversity.png
#  Scatterplot_Triglycerides_vs_alpha.png
#  Scatterplot_TyG_vs_alpha
# ScatterPlot_VLDL_vs_alpha.png
# Scatterplot_GGT_vs_shannon.png
#  Scatterplot_IMC_vs_alpha.png





source ("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/AgriculturaUrbana_Analises/PS_URBANAGRI_5004_AlphaDiversityAnalyses_14_08_2024.nb.R")


# ======================================= #
# 5. Beta Diversity #N=105
# ======================================= #

#This script is for beta diversity analyses. With inflamatory and metabolic parameters, with, with diet, NOVA classification


source ("C:\Users\polia\OneDrive\Desktop\EstatisticaR\AgrUrbana\16S_AgriUrbana\AgriculturaUrbana_Analises\PS_URBANAGRI_5005_BetaDiversityAnalyses_29_04_2025.nb.R")

#(Output objects: em "C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Para_Dissertação/..." :  permanova_significant_variables.png, permanova_unweighted_significant_variables.png, unweighted_unifrac_pcoa.png, volcanoPlot_significant_variables.png, volcanoPlot_significant_variables_withLabel.png, volcanoPlot_unweighted_significant_variables.png, volcanoPlot_unweighted_significant_variables_withlabel.png, weighted_unifrac_pcoa.png")



