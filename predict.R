args=commandArgs(trailingOnly=TRUE)
# print(length(args))
if (length(args) != 3 )
  stop("Invalid number of arguments to Rscript.")

script_path <- args[1]
input_file <- args[2]
output_path <- args[3]
# print(script_path)
# print(input_file)
# print(output_path)

convertRI <- function(diff){
  if(diff < 0.1)
    return(1)

  return(ceiling(round(diff*10, digits = 0)/2))
}

library(reshape2, quietly=TRUE)
library(caret, quietly=TRUE)

A_col <- unlist(read.table(paste0(script_path, "A_col.txt")))
svm_L <- readRDS(paste0(script_path, "svm_L_model.rds"))

In <- read.table(input_file, comment.char = "", sep = "\t", header = TRUE)

if(length(unique(In$SAMPLE))>1)
  stop("More than 1 samples are found in the VCF file.")

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
  }
}


pred_prob <- predict(svm_L, data_2, type="prob")

prob <- unname(unlist(pred_prob))
diff <- max(prob)-max(prob[prob!=max(prob)]) #highest - second highest
RI <- convertRI(diff)

out <- pred_prob
colnames(out)[1:3] <- c("MDR", "Susceptible", "XDR")
out[,c("MDR", "Susceptible", "XDR")] <- round(out[,c("MDR", "Susceptible", "XDR")], digits = 4)
if(predict(svm_L, data_2)== "X") {
  out$Class <- "XDR"
} else if(predict(svm_L, data_2)== "M") {
  out$Class <- "MDR"
} else {
  out$Class <- "Susceptible"
}

out$Sample <- Ind_l$SAMPLE
out$RI <- RI
out <- out[c(5,1:4,6)]
cat("\n")
print(out)
cat("\n")

write.table(out, file = paste0(output_path, "prediction.tsv"), row.names = FALSE)

# cat("\n")
# cat("\tSAMPLE =", Ind_l$SAMPLE, "\n")
# cat("\tPredicted class =", suppressWarnings(names(sort(pred_prob_svm_R_ind, decreasing=TRUE)))[1], "\n")
# cat("\tProbability =", max(prob), "\n")
# cat("\tRelability Index (RI) =", RI, "\n\n")

