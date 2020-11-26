precision_recall<- function (path, theta, verbose = TRUE, plot = FALSE, flip = T) {
  gcinfo(verbose = FALSE)
  ROC = list()
  theta = as.matrix(theta)
  d = ncol(theta)
  pos.total = sum(theta != 0)
  neg.total = d * (d - 1) - pos.total
  if (verbose)
    message("Computing F1 scores, false positive rates and true positive rates....", appendLF=FALSE)
  ROC$prec = rep(0, length(path))
  ROC$rec  = rep(0, length(path))
  ROC$F1 = rep(0, length(path))
  for (r in 1:length(path)) {
    tmp = as.matrix(path[[r]])
    tp.all = (theta != 0) * (tmp != 0)
    diag(tp.all) = 0
    ROC$tp[r] <- sum(tp.all != 0)/pos.total
    fp.all = (theta == 0) * (tmp != 0)
    diag(fp.all) = 0
    ROC$fp[r] <- sum(fp.all != 0)/neg.total
    fn = 1 - ROC$tp[r]
    precision = ROC$tp[r]/(ROC$tp[r] + ROC$fp[r])
    recall = ROC$tp[r]/(ROC$tp[r] + fn)
    
    # Correction
    tp = sum(tp.all)
    fp = sum(fp.all)
    p = sum(theta != 0)
    
    denominateur = tp + fp
    denominateur = max(denominateur, 1e-6)
    
    precision = tp / (denominateur)
    recall = tp / p
    
    ROC$prec[r] <- precision
    ROC$rec[r]  <- recall
    
    ROC$F1[r] = 2 * precision * recall/(precision + recall)
    if (is.na(ROC$F1[r]))
      ROC$F1[r] = 0
  }
  if (verbose)
    message("done.")
  rm(precision, recall, tp.all, fp.all, path, theta, fn)
  gc()
  if(flip) {
    ord.p = order(ROC$rec, na.last=NA)
    tmp2 = ROC$prec[ord.p]
    tmp1 = ROC$rec[ord.p]
    tmp2 = c(1, tmp2)
    tmp1 = c(0, tmp1)
  }
  else{
    ord.p = order(ROC$prec, na.last=NA)
    tmp1 = ROC$prec[ord.p]
    tmp2 = ROC$rec[ord.p]
  }
  
  if (plot) {
    par(mfrow = c(1, 1))
    plot(tmp1, tmp2, type = "b", main = "PR Curve", xlab = "Precision",
         ylab = "Recall", ylim = c(0, 1))
  }
  ROC$AUC = sum(diff(tmp1) * (tmp2[-1] + tmp2[-length(tmp2)]))/2
  rm(ord.p, tmp1, tmp2)
  gc()
  class(ROC) = "roc"
  return(ROC)
}


max_off_diagonal_value = function(S) {
  S_diag_off = S
  diag(S_diag_off) = 0
  max_value = max(S_diag_off)
  return(max_value)
}

lamda_path = function(lamda_max = 1, lamda_ratio = 0.1, nlambda=12){
  pen = seq(0,nlambda - 1, 1) / (nlambda - 1)
  lamda_min = lamda_max * lamda_ratio
  pen = rep(lamda_max, nlambda) * (lamda_min)^pen
  return(pen)
}

