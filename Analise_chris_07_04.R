
holder <- data.frame(metadados.saude.alpha)

#remove repeated samples
holder <- metadados.saude.alpha[-c(84,96),]

#rename rownames
rownames(holder) <- holder$Sample.id
holder$Sample.id <- NULL

#make na = 0 and reast = 1

# rememeber alpha div is from col 1 to col 6

holder.bin <- apply(X = holder, MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1))
my.rows <- rowSums(holder.bin[,-c(1:6)]) > 4


my.rows <- rowSums(holder.bin[,c("INSULINA", "PCR", "IL17A", "IFNGamma", "TNF", "IL2", "IL4", "IL6", "IL10", "VLDL", "LDL", "HDL", "BMI", "GLICOSE", "HbA1c", "COLESTEROL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "Systolic", "Diastolic")]) >21


library(pheatmap)

pheatmap(mat = holder.bin, clustering_method = "ward.D2")


#calculo VLDL

metadados.all$VLDL_calc <- metadados.all$TRIGLICERIDES / 5


#biplot(princomp(scale(holder[my.rows,c("INSULINA", "PCR", "IL17A", "IFNGamma", "TNF", "IL2", "IL4", "IL6", "IL10", "VLDL", "LDL", "HDL", "BMI", "GLICOSE", "HbA1c", "COLESTEROL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "Systolic", "Diastolic")])))

#princomp(scale(holder[my.rows,c("INSULINA", "PCR", "IL17A", "IFNGamma", "TNF", "IL2", "IL4", "IL6", "IL10", "VLDL", "LDL", "HDL", "BMI", "GLICOSE", "HbA1c", "COLESTEROL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "Systolic", "Diastolic")]))

#summary(as.matrix(dist(scale(holder[my.rows,c("INSULINA", "PCR", "IL17A", "IFNGamma", "TNF", "IL2", "IL4", "IL6", "IL10", "VLDL", "LDL", "HDL", "BMI", "GLICOSE", "HbA1c", "COLESTEROL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "Systolic", "Diastolic")]))))

#as.matrix(dist(scale(holder[my.rows,c("INSULINA", "PCR", "IL17A", "IFNGamma", "TNF", "IL2", "IL4", "IL6", "IL10", "VLDL", "LDL", "HDL", "BMI", "GLICOSE", "HbA1c", "COLESTEROL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "Systolic", "Diastolic")])))

#scale(holder[my.rows,c("INSULINA", "PCR", "IL17A", "IFNGamma", "TNF", "IL2", "IL4", "IL6", "IL10", "VLDL", "LDL", "HDL", "BMI", "GLICOSE", "HbA1c", "COLESTEROL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "Systolic", "Diastolic")])

#holder[my.rows,c("INSULINA", "PCR", "IL17A", "IFNGamma", "TNF", "IL2", "IL4", "IL6", "IL10", "VLDL", "LDL", "HDL", "BMI", "GLICOSE", "HbA1c", "COLESTEROL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "Systolic", "Diastolic")]

#holder.bin[my.rows,c("INSULINA", "PCR", "IL17A", "IFNGamma", "TNF", "IL2", "IL4", "IL6", "IL10", "VLDL", "LDL", "HDL", "BMI", "GLICOSE", "HbA1c", "COLESTEROL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "Systolic", "Diastolic")]

#sort(colSums(holder.bin[,c("INSULINA", "PCR", "IL17A", "IFNGamma", "TNF", "IL2", "IL4", "IL6", "IL10", "VLDL", "LDL", "HDL", "BMI", "GLICOSE", "HbA1c", "COLESTEROL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "Systolic", "Diastolic")]))

#rowSums(holder.bin[,-c(1:6)]) > 4

#holder.bin <- apply(X = holder, MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1))

#pheatmap(apply(X = holder[rowSums(apply(X = holder[,-c(1:6, 8, 10)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)))<26,-c(1:6, 8, 10)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)), clustering_method = "ward.D2")

#colSums(apply(X = holder[rowSums(apply(X = holder[,-c(1:6, 8, 10)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)))<26,-c(1:6, 8, 10)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)))

#table(sort(rowSums(apply(X = holder[,-c(1:6, 8, 10)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)))))

#colSums(apply(X = holder[rowSums(apply(X = holder[,-c(1:6)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)))<26,-c(1:6)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)))

#rowSums(apply(X = holder[,-c(1:6)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)))<26

#sort(rowSums(apply(X = holder[,-c(1:6)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1))))

#pheatmap(apply(X = holder[,-c(1:6)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)), clustering_method = "ward.D2")

#heatmap.2(apply(X = holder, MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)))

#heatmap.2(apply(X = holder, MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)))

#
