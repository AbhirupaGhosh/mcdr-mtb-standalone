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

Input <- read.table(input_file, sep = "\t", comment.char = "", header = TRUE)
Ind_batch <- Input[,c("POS", "REF", "ALT", "QUAL", "SAMPLE", "AO.1", "DP.1")]
Ind_batch$AO.1[Ind_batch$AO.1=="."] <- 0
Ind_batch$DP.1[Ind_batch$DP.1=="."] <- 0
Ind_batch <- transform(Ind_batch, DP.1 = as.numeric(DP.1))
Ind_batch <- transform(Ind_batch, AO.1 = as.numeric(AO.1))
Ind_batch$MUT <- paste(Ind_batch$REF, Ind_batch$POS, Ind_batch$ALT, sep = "")
Ind_batch$VAL <- Ind_batch$AO.1/Ind_batch$DP.1
Ind_batch$VAL[Ind_batch$VAL == "NaN"] <- 0.0000
Ind_batch$VAL <- round(Ind_batch$VAL, digits = 4)
Ind_batch_l <- dcast(Ind_batch[,c("SAMPLE", "MUT", "VAL")], SAMPLE~MUT, value.var = "VAL")

A_col <- unlist(read.table(paste0(script_path, "A_col.txt")))
svm_R <- readRDS(paste0(script_path, "svm_R_model.rds"))

data_2_bat <- data.frame(matrix(ncol = length(A_col), nrow = nrow(Ind_batch_l)))
colnames(data_2_bat) <- A_col

for(j in 1:nrow(data_2_bat)){
  for(i in colnames(data_2_bat)){
    if(i %in% colnames(Ind_batch_l)){
      data_2_bat[j,i] <- Ind_batch_l[j,i]
    }else{
      data_2_bat[j,i] <- 0.00
    }}}

pred_prob_svm_R_ind_bat <- predict(svm_R, data_2_bat, type="prob")
pred_svm_R_ind_bat <- predict(svm_R, data_2_bat)

pred_svm_R_ind_bat <- unlist(pred_prob_svm_R_ind_bat)
output <- data.frame()
for(i in 1:nrow(pred_prob_svm_R_ind_bat)){
 n <- which.max(pred_prob_svm_R_ind_bat[i,])
 output[i,1] <- Ind_batch_l$SAMPLE[i]
 output[i,2] <- colnames(pred_prob_svm_R_ind_bat)[n]
 output[i,3] <- max(pred_prob_svm_R_ind_bat[i,])- max(pred_prob_svm_R_ind_bat[i,-n])
 output[i,4] <- round(output[i,3]*10, digits = 0)
}
colnames(output) <- c("SAMPLE", "PREDICTION", "DIFF", "RI")

cat("\n")
print(output)
cat("\n")
