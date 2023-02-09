args=commandArgs(trailingOnly=TRUE)
# print(length(args))
if (length(args) != 2 )
  stop("Invalid number of arguments to Rscript.")

script_path <- args[1]
input_file <- args[2]

# print(script_path)
# print(input_file)

library(reshape2, quietly=TRUE)
library(caret, quietly=TRUE)

A_col <- unlist(read.table(paste0(script_path, "A_col.txt")))
svm_R <- readRDS(paste0(script_path, "svm_R_model.rds"))

In <- read.table(input_file, comment.char = "", sep = "\t", header = TRUE)
Ind <- In[,c("POS", "REF", "ALT", "QUAL", "SAMPLE", "AO", "DP")]
Ind$AO[Ind$AO=="."] <- 0
Ind$DP[Ind$DP=="."] <- 0
Ind <- transform(Ind, DP = as.numeric(DP))
Ind <- transform(Ind, AO = as.numeric(AO))
Ind$MUT <- paste(Ind$REF, Ind$POS, Ind$ALT, sep = "")
Ind$VAL <- Ind$AO/Ind$DP
Ind$VAL[Ind$VAL == "NaN"] <- 0.0000
Ind$VAL <- round(Ind$VAL, digits = 4)
Ind_l <- dcast(Ind[,c("SAMPLE", "MUT", "VAL")], SAMPLE~MUT, value.var = "VAL")

data_2 <- data.frame(matrix(ncol = length(A_col), nrow = nrow(Ind_l)))
colnames(data_2) <- A_col

for(i in colnames(data_2)){
  if(i %in% colnames(Ind_l)){
    data_2[1,i] <- Ind_l[1,i]
  }else{
    data_2[1,i] <- 0.00
  }}


pred_prob_svm_R_ind <- predict(svm_R, data_2, type="prob")
pred_svm_R_ind <- predict(svm_R, data_2)
prob <- unname(unlist(pred_prob_svm_R_ind))
diff <- max(prob)-max(prob[prob!=max(prob)]) #highest - second highest
RI <- round(diff*10, digits = 0)

cat("\n")
cat("\tSAMPLE =", Ind_l$SAMPLE, "\n")
cat("\tClass =", suppressWarnings(names(sort(pred_prob_svm_R_ind, decreasing=TRUE)))[1], "\n")
cat("\tProbability =", max(prob), "\n")
cat("\tRelability Index (RI) =", RI, "\n")
cat("\tSAMPLE", "\t", "Class", "\t", "Probability", "\t", "Relability Index","\n\n")
