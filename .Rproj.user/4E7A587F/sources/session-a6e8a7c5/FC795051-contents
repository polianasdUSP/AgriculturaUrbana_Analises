


holder <- data.frame(metadados.saude.alpha)

#remove repeated samples
holder <- metadados.saude.alpha[-c(84,96),]

holder <- metadados.saude.alpha

#rename rownames
rownames(holder) <- holder$Sample.id
holder$Sample.id <- NULL

#make na = 0 and rest = 1

# remember alpha div is from col 38 to col 43

holder.bin <- apply(X = holder, MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1))
my.rows <- rowSums(holder.bin[,-c(9:10, 31:43)]) 


my.rows <- rowSums(holder.bin[,c("INSULINA", "PCR", "IL17A", "IFNGamma", "TNF", "IL2", "IL4", "IL6", "IL10", "VLDL", "LDL", "HDL", "IMC", "GLICOSE", "HbA1c", "COLESTEROL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "Systolic", "Diastolic")]) >21

sort(my.rows) > 7


library(pheatmap)

pheatmap(mat = holder.bin , clustering_method = "ward.D2")

pheatmap(mat = holder(,9:10, 31:43), clustering_method = "ward.D2")

#calculo VLDL

metadados.all$VLDL_calc <- metadados.all$TRIGLICERIDES / 5

metadados.all[, c("VLDL", "VLDL_calc")]


#biplot(princomp(scale(holder[my.rows,c("INSULINA", "PCR", "IL17A", "IFNGamma", "TNF", "IL2", "IL4", "IL6", "IL10", "VLDL", "LDL", "HDL", "BMI", "GLICOSE", "HbA1c", "COLESTEROL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "Systolic", "Diastolic")])))

#princomp(scale(holder[my.rows,c("INSULINA", "PCR", "IL17A", "IFNGamma", "TNF", "IL2", "IL4", "IL6", "IL10", "VLDL", "LDL", "HDL", "BMI", "GLICOSE", "HbA1c", "COLESTEROL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "Systolic", "Diastolic")]))

#summary(as.matrix(dist(scale(holder[my.rows,c("INSULINA", "PCR", "IL17A", "IFNGamma", "TNF", "IL2", "IL4", "IL6", "IL10", "VLDL", "LDL", "HDL", "BMI", "GLICOSE", "HbA1c", "COLESTEROL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "Systolic", "Diastolic")]))))

#as.matrix(dist(scale(holder[my.rows,c("INSULINA", "PCR", "IL17A", "IFNGamma", "TNF", "IL2", "IL4", "IL6", "IL10", "VLDL", "LDL", "HDL", "BMI", "GLICOSE", "HbA1c", "COLESTEROL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "Systolic", "Diastolic")])))

#scale(holder[my.rows,c("INSULINA", "PCR", "IL17A", "IFNGamma", "TNF", "IL2", "IL4", "IL6", "IL10", "VLDL", "LDL", "HDL", "BMI", "GLICOSE", "HbA1c", "COLESTEROL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "Systolic", "Diastolic")])

#holder[my.rows,c("INSULINA", "PCR", "IL17A", "IFNGamma", "TNF", "IL2", "IL4", "IL6", "IL10", "VLDL", "LDL", "HDL", "BMI", "GLICOSE", "HbA1c", "COLESTEROL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "Systolic", "Diastolic")]

#holder.bin[my.rows,c("INSULINA", "PCR", "IL17A", "IFNGamma", "TNF", "IL2", "IL4", "IL6", "IL10", "VLDL", "LDL", "HDL", "BMI", "GLICOSE", "HbA1c", "COLESTEROL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "Systolic", "Diastolic")]

sort(colSums(holder.bin[,c("INSULINA", "PCR", "IL17A", "IFNGamma", "TNF", "IL2", "IL4", "IL6", "IL10", "VLDL", "LDL", "HDL", "IMC", "GLICOSE", "HbA1c", "COLESTEROL", "TRIGLICERIDES", "TGO", "TGP", "GGT", "Systolic", "Diastolic")]))

rowSums(holder.bin[,-c(31:43)]) > 4

holder.bin <- apply(X = holder, MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1))


pheatmap(apply(X = holder[rowSums(apply(X = holder[,-c(8:9, 12:13, 15:16, 30:42)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1))) >7,-c(8:9, 12:13, 15:16, 30:42)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)), clustering_method = "ward.D2") 



pheatmap(
  apply(
    X = holder[
      apply(holder[, -c(8:9, 12:13, 15:16, 30:42)], 1, function(x) sum(!is.na(x))) > 7,
      -c(8:9, 12:13, 15:16, 30:42)
    ],
    2,
    function(x) ifelse(is.na(x), 0, 1)
  ),
  clustering_method = "ward.D2"
)



rowSums(apply(X = holder[rowSums(apply(X = holder[,-c(8:9, 12:13, 15:16, 30:42)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)))<13,-c(8:9, 12:13, 15:16, 30:42)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)))

table(sort(rowSums(apply(X = holder[,-c(8:9, 12:13, 15:16, 30:42)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)))))

rowSums(apply(X = holder[rowSums(apply(X = holder[,-c(8:9, 12:13, 15:16, 30:42)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1))),-c(8:9, 12:13, 15:16, 30:42)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)))

#rowSums(apply(X = holder[,-c(1:6)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)))<26

#sort(rowSums(apply(X = holder[,-c(1:6)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1))))

pheatmap(apply(X = holder[,-c(8:9, 12:13, 15:16, 30:42)], MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)), clustering_method = "ward.D2")

colnames(holder.bin)[colSums(holder.bin) == nrow(holder.bin)]
#heatmap.2(apply(X = holder, MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)))

#heatmap.2(apply(X = holder, MARGIN = 2, FUN = function(x)ifelse(is.na(x),0,1)))

#
