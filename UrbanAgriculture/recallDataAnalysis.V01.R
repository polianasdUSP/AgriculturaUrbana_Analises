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

# libraries needed
library(dplyr)
library(readxl)

### analysis ###

# import the TBCA
tbca.raw <- read_excel("Dados_composicao_TBCA_Hoffmann11022022v6.xlsx")
# convert negative values that indicate traces of not available data etc to zero
tbca <- tbca.raw %>%
  mutate(across(where(is.numeric), ~ replace(., . < 0, 0)))


# import the recall data
recall.raw <- read_excel("Planilha_UNICA_R24H_Projeto_Agricultura.xlsx")
# convert ID to char
recall.raw$IdVoluntario <- as.character(recall.raw$IdVoluntario)


recall <- merge(x = recall.raw[c("IdVoluntario", "R24H", "Quantidade", "TBCA_code")], y = tbca, by.x = "TBCA_code", by.y = "cod_alimento", all.x =T)
recall$convertionFactor <- recall$Quantidade/100

recall <- cbind(recall[1:5], recall[12:(ncol(recall)-1)]*recall$convertionFactor)


recall.calc <- recall[-c(1,4,5)] %>% 
  group_by(IdVoluntario, R24H) %>%
  summarise_all(sum, na.rm=T) %>%
  summarise(across(-R24H, mean))


write.table(x = recall.calc, file = "recall.calc.txt", col.names = T)






