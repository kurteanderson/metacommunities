library(sizeSpectra)
library(tidyverse)


setwd("G:/My Drive/Stream Metacommunities/Local Models")

#### run small community
r.in <- 0.75
fec.in <- 25
temp <- 298.15

source("Open_model_body_size_params.R")

mle.store.small <- rep(NA,36)

ss_list_small <- list()

for (j in 13:48) {

n <- results_both[,j]

bs.test.small <- data.frame(mass,n)

bs.test.small <- bs.test.small[bs.test.small$n > 0,]

for (k in 1:length(bs.test.small$n)) {

  if (k == 1) {
    m2.small <- rep(bs.test.small$mass[k],bs.test.small$n[k])
  }
  else {
    m2.small <- c(m2.small,rep(bs.test.small$mass[k],bs.test.small$n[k]))
  }
}

x_binned <- binData(counts = bs.test.small, binWidth = "2k")

x_binned$binVals <- x_binned$binVals %>% filter(binCount !=0)

num_bins <- nrow(x_binned$binVals)

# bin breaks are the minima plus the max of the final bin:
bin_breaks <- c(dplyr::pull(x_binned$binVals, binMin), dplyr::pull(x_binned$binVals, binMax)[num_bins])

bin_counts <- dplyr::pull(x_binned$binVals, binCount)

mle_bin <-  calcLike(negLL.PLB.binned, p = -1.5, w = bin_breaks, d = bin_counts, J = length(bin_counts), vecDiff = 1, suppress.warnings = TRUE)             # increase this if hit a bound

mle.store.small[j-12] <- mle_bin$MLE

mle_bin_small <- mle_bin

binTibble <- dplyr::select(x_binned$binVals, wmin = binMin, wmax = binMax, Number = binCount)

x.PLB <- exp(seq(log(min(binTibble$wmin)), log(max(binTibble$wmax)), length = 10000))

B.PLB <- dPLB(x.PLB, b = mle_bin_small$MLE, xmin=min(x.PLB), xmax=max(x.PLB)) * sum(binTibble$Number) * x.PLB

ss_list_small[[j-12]] <- data.frame(x = x.PLB, y = B.PLB, temp = "25")

}

# LBN_bin_plot(
#   binValsTibble = x_binned$binVals,
#   b.MLE = mle_bin_small$MLE,
#   b.confMin = mle_bin_small$conf[1],
#   b.confMax = mle_bin_small$conf[2],
#   log.xy = "xy",
#   xLab = expression(paste("Body mass ", italic(x), "(ug)")),
#   yLab = "Normalised biomass",
#   leg.text = "25 °C")

#### run large community

r.in <- 0.75
fec.in <- 25
temp <- 288.15

source("Open_model_body_size_params.R")

mle.store.large <- rep(NA,36)
 
ss_list_large <- list()

for (j in 13:48) {

  n <- results_both[,j]

  bs.test.large <- data.frame(mass,n)

  bs.test.large <- bs.test.large[bs.test.large$n > 0,]


  for (k in 1:length(bs.test.large$n)) {

    if (k == 1) {
      m2.large <- rep(bs.test.large$mass[k],bs.test.large$n[k])
    }
    else {
      m2.large <- c(m2.large,rep(bs.test.large$mass[k],bs.test.large$n[k]))
    }
  }

  x_binned <- binData(counts = bs.test.large, binWidth = "2k")
  
  x_binned$binVals <- x_binned$binVals %>% filter(binCount !=0)
  
  num_bins <- nrow(x_binned$binVals)
  
  # bin breaks are the minima plus the max of the final bin:
  bin_breaks <- c(dplyr::pull(x_binned$binVals, binMin), dplyr::pull(x_binned$binVals, binMax)[num_bins])
  
  bin_counts <- dplyr::pull(x_binned$binVals, binCount)
  
  mle_bin <-  calcLike(negLL.PLB.binned, p = -1.5, w = bin_breaks, d = bin_counts, J = length(bin_counts), vecDiff = 1, suppress.warnings = TRUE)             # increase this if hit a bound
  
  mle.store.large[j-12] <- mle_bin$MLE
  
  mle_bin_large <- mle_bin
  
  binTibble <- dplyr::select(x_binned$binVals, wmin = binMin, wmax = binMax, Number = binCount)
  
  x.PLB <- exp(seq(log(min(binTibble$wmin)), log(max(binTibble$wmax)), length = 10000))
  
  B.PLB <- dPLB(x.PLB, b = mle_bin_large$MLE, xmin=min(x.PLB), xmax=max(x.PLB)) * sum(binTibble$Number) * x.PLB
  
  ss_list_large[[j-12]] <- data.frame(x = x.PLB, y = B.PLB, temp = "15")

}

l_small<-as.data.frame(do.call(rbind,ss_list_small))
l_large<-as.data.frame(do.call(rbind,ss_list_large))

dataset_ss<-data.frame(rbind(l_small,l_large))

ggplot(dataset_ss, aes(x, y, color = as.factor(temp))) +
  theme_classic() +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", alpha = 0.1, se = FALSE) +
  scale_x_log10() +
  scale_y_log10(limits = c(2,900)) +
  labs(title = "a)", x = "Body Mass (ug)", y = "Density", fill = "temp", color = "temp") +
  theme(legend.position = "inside", legend.position.inside = c(0.01, 0.01), legend.justification = c("left", "bottom"), legend.box.just = "left") +
  scale_color_manual(NULL, values = c("#0072B2", "#E69F00"), labels = c("15 °C", "25 °C"))

# LBN_bin_plot(
#   binValsTibble = x_binned$binVals,
#   b.MLE = mle_bin_large$MLE,
#   b.confMin = mle_bin_large$conf[1],
#   b.confMax = mle_bin_large$conf[2],
#   log.xy = "xy",
#   xLab = expression(paste("Body mass ", italic(x), "(ug)")),
#   yLab = "Normalised biomass",
#   leg.text = "15 °C")

df1 <- data.frame(
  comm.size = factor(rep("small",length(mle.store.small[1:36]))),
  lambda = mle.store.small[1:36]
)

df2 <- data.frame(
  comm.size = factor(rep("large",length(mle.store.large[1:36]))),
  lambda = mle.store.large[1:36]
)

df <- rbind(df1,df2)

sizes <- c("small", "large")

ggplot(df, aes(x=lambda, fill=comm.size, color=comm.size)) +
  theme_classic() +
  geom_histogram(position = "identity", alpha = 0.5, binwidth = 0.02) +
  geom_density(alpha = 0.3) +
  labs(x = expression(lambda), fill = "Community Size", color = "Community Size") +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text = element_text(size = 16), axis.title.x = element_text(size = 16)) +
  scale_color_manual(NULL, values = c("#E69F00", "#0072B2"), labels = c("25 °C", "15 °C")) +
  scale_fill_manual(NULL, values = c("#E69F00", "#0072B2"), labels = c("25 °C", "15 °C"))


