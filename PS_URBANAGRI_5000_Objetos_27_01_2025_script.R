install.packages("remotes")
remotes::install_github("jbisanz/qiime2R")




library(qiime2R)


#======================================= #
#  All Objects for analyses
# ======================================= #

metadados <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados.csv", stringsAsFactors = FALSE)

Cluster <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/Clusters - Sheet1.csv", stringsAsFactors = FALSE)

metadados_MetS <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados_MetS.csv", stringsAsFactors = FALSE)

metadados_parametros_metabolicos <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados_parametros_metabolicos.csv", stringsAsFactors = FALSE)

weighted_unifrac <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/weighted_unifrac.csv", stringsAsFactors = FALSE, row.names = 1, )

alpha.all <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/alpha.all.csv", stringsAsFactors = FALSE)

alpha.evenness <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/alpha.evenness.csv", stringsAsFactors = FALSE)

alpha.pd <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/alpha.pd.csv", stringsAsFactors = FALSE)

alpha.richness <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/alpha.richness.csv", stringsAsFactors = FALSE)

alpha.shannon <- read.csv ("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/alpha.shannon.csv", stringsAsFactors = FALSE)


metadados.alpha.all <- read.csv ("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados.alpha.all.csv", stringsAsFactors = FALSE)


metadados <- read.csv ("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados.csv", stringAsFactors = FALSE)


metadados.shannon <- read.csv("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadados.shannon.csv", stringsAsFactors = FALSE)


unweighted_unifrac <- read.csv ("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/unweighted_unifrac.csv", stringsAsFactors = FALSE)

metadata__ASVs_Saude_Dieta_NOVA <-read.csv ("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/metadata__ASVs_Saude_Dieta_NOVA.csv", stringsAsFactors = FALSE)

SVs_normalized_log_transformed <-read.csv ("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/SVs_normalized_log_transformed.csv", stringsAsFactors = FALSE)

diet_residual <-read.csv ("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Planilhas_UrbanAgri/diet_data.csv", stringsAsFactors = FALSE)




# Importar os arquivos QZA
otu_table <- read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Phyloseq/table.qza")$data

taxonomy <- read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Phyloseq/taxonomy_silva.qza")$data

tree <- read_qza("C:/Users/polia/OneDrive/Desktop/EstatisticaR/AgrUrbana/16S_AgriUrbana/Phyloseq/rooted-tree.qza")$data

