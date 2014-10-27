# auc
library("AUC")
pauc <- function(test.vector, lab){
  roc.out <- roc(1-test.vector, lab)
  roc.ind <- sum(roc.out$fpr<=.1)
  roc.min <- roc.out$cutoffs[roc.ind]
  pauc_out <- auc(roc.out, min =roc.min)
  return(pauc_out)
}
