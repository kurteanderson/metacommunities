mle.store <- rep(NA,time)

for (j in 1:time) {
  n <- results_both[,j]
  
  bs.test <- data.frame(mass,n)
  
  bs.test <- bs.test[bs.test$n > 0,]
  
  for (k in 1:length(bs.test$n)) {
    
    if (k == 1) {
      m2 <- rep(bs.test$mass[k],bs.test$n[k])
    }
    else {
      m2 <- c(m2,rep(bs.test$mass[k],bs.test$n[k]))
    }
  }
  
  x_binned <- binData(counts = bs.test, binWidth = "2k")
  
  x_binned$binVals <- x_binned$binVals %>% filter(binCount !=0)
  
  num_bins <- nrow(x_binned$binVals)
  
  # bin breaks are the minima plus the max of the final bin:
  bin_breaks <- c(dplyr::pull(x_binned$binVals, binMin), dplyr::pull(x_binned$binVals, binMax)[num_bins])
  
  bin_counts <- dplyr::pull(x_binned$binVals, binCount)
  
  mle_bin <-  calcLike(negLL.PLB.binned, p = -1.5, w = bin_breaks, d = bin_counts, J = length(bin_counts), vecDiff = 1, suppress.warnings = TRUE)             # increase this if hit a bound
  
  mle.store[j] <- mle_bin$MLE
}

ss.mle.median <- median(mle.store)
ss.mle.range  <- max(mle.store)-min(mle.store)