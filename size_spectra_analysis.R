library(sizeSpectra)
library(tidyverse)

setwd("G:/My Drive/Stream Metacommunities/Local Models")

#### run small community
r.in <- 0.25
fec.in <- 5
temp <- 298.15

source("Open_model_body_size_params.R")

mle.store.small <- rep(NA,36)

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

mle_lambda_2 <- calcLike(
  negLL.fn = negLL.PLB, # continuous estimates of all individuals
  x = m2.small, # the vector of data
  xmin = min(m2.small), # the minimum body size
  xmax = max(m2.small), # the maximum body size
  n = length(m2.small), # the number of observations
  sumlogx = sum(log(m2.small)), # sum of log-transformed data
  p = -1.5) # starting point, arbitrary number

mle.store.small[j] <- mle_lambda_2$MLE

mle_lambda_2_small <- mle_lambda_2



}

# MLE.plot(x = m2.small, # vector of simulated body sizes
#          b = mle_lambda_2$MLE, #lambda estimate
#          confVals = c(mle_lambda_2$conf[1],# confidence interval
#                       mle_lambda_2$conf[2]),
#          panel = "b", #This option includes the estimate and CI
#          log="xy") # you can change this to just x
#
#hist(mle.store.small, xlab=expression(lambda), main="Small community")


#### run large community

r.in <- 0.75
fec.in <- 25
temp <- 288.15

source("Open_model_body_size_params.R")

mle.store.large <- rep(NA,36)

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

  mle_lambda_2 <- calcLike(
    negLL.fn = negLL.PLB, # continuous estimates of all individuals
    x = m2.large, # the vector of data
    xmin = min(m2.large), # the minimum body size
    xmax = max(m2.large), # the maximum body size
    n = length(m2.large), # the number of observations
    sumlogx = sum(log(m2.large)), # sum of log-transformed data
    p = -1.5) # starting point, arbitrary number

  mle.store.large[j] <- mle_lambda_2$MLE

  mle_lambda_2_large <- mle_lambda_2


}

# MLE.plot(x = m2.large, # vector of simulated body sizes
#           b = mle_lambda_2$MLE, #lambda estimate
#           confVals = c(mle_lambda_2$conf[1],# confidence interval
#                        mle_lambda_2$conf[2]),
#           panel = "b", #This option includes the estimate and CI
#           log="xy") # you can change this to just x
#
#hist(mle.store.large, xlab=expression(lambda), main="Large community")


#### post processing

df1 <- data.frame(
  comm.size = factor(rep("small",length(mle.store.small[13:48]))),
  lambda = mle.store.small[13:48]
)

df2 <- data.frame(
  comm.size = factor(rep("large",length(mle.store.large[13:48]))),
  lambda = mle.store.large[13:48]
)

df <- rbind(df1,df2)

sizes <- c("small", "large")

ggplot(df, aes(x=lambda, fill=comm.size, color=comm.size)) +
  theme_classic() +
  geom_histogram(position="identity", alpha=0.5, binwidth = 0.02) +
  geom_density(alpha=0.3) +
  labs(title="b)", x=expression(lambda), y = "Density", fill = "Community Size", color = "Community Size") +
  theme(legend.position = "inside", legend.position.inside = c(.99, .99), legend.justification = c("right", "top"), legend.box.just = "left") +
  scale_color_manual(NULL, values = c("#E69F00", "#0072B2"), labels = c("25 °C", "15 °C")) +
  scale_fill_manual(NULL, values = c("#E69F00", "#0072B2"), labels = c("25 °C", "15 °C"))

mle.df1 <- data.frame("mass" = sort(m2.small, decreasing = TRUE), "num" = 1:length(m2.small), "size" = "small")

mle.df2 <- data.frame("mass" = sort(m2.large, decreasing = TRUE), "num" = 1:length(m2.large), "size" = "large")

mle.df <- rbind(mle.df1,mle.df2)

min(mle.df1$mass)
max(mle.df1$mass)
mle.store.small[48]
length(mle.df1$mass)

min(mle.df2$mass)
max(mle.df2$mass)
mle.store.large[48]
length(mle.df2$mass)

sizes <- c("small", "large")

ggplot(mle.df, aes(x = mass, y = num, fill = size, color = size)) +
  theme_classic() +
  geom_point() +
  geom_function(data = mle.df %>% filter(size == "small"), fun = function(x, b, xmin, xmax) ((1-pPLB(x,b,xmin,xmax))*length(mle.df1$mass)),
                args = list(b = mle_lambda_2_small$MLE, xmin = min(mle.df1$mass), xmax = max(mle.df1$mass))) +
  geom_function(data = mle.df %>% filter(size == "small"), fun = function(x, b, xmin, xmax) ((1-pPLB(x,b,xmin,xmax))*length(mle.df1$mass)),
                args = list(b = mle_lambda_2_small$conf[1], xmin = min(mle.df1$mass), xmax = max(mle.df1$mass)), linetype = "dashed") +
  geom_function(data = mle.df %>% filter(size == "small"), fun = function(x, b, xmin, xmax) ((1-pPLB(x,b,xmin,xmax))*length(mle.df1$mass)),
                args = list(b = mle_lambda_2_small$conf[2], xmin = min(mle.df1$mass), xmax = max(mle.df1$mass)), linetype = "dashed") +
  geom_function(data = mle.df %>% filter(size == "large"), fun = function(x, b, xmin, xmax) ((1-pPLB(x,b,xmin,xmax))*length(mle.df2$mass)),
                args = list(b = mle_lambda_2_large$MLE, xmin = min(mle.df2$mass), xmax = max(mle.df2$mass))) +
  geom_function(data = mle.df %>% filter(size == "large"), fun = function(x, b, xmin, xmax) ((1-pPLB(x,b,xmin,xmax))*length(mle.df2$mass)),
                args = list(b = mle_lambda_2_large$conf[1], xmin = min(mle.df2$mass), xmax = max(mle.df2$mass)), linetype = "dashed") +
  geom_function(data = mle.df %>% filter(size == "large"), fun = function(x, b, xmin, xmax) ((1-pPLB(x,b,xmin,xmax))*length(mle.df2$mass)),
                args = list(b = mle_lambda_2_large$conf[2], xmin = min(mle.df2$mass), xmax = max(mle.df2$mass)), linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Example Size Spectra", x = "Body Mass (ug)", y = "Number of Values", fill = "Community Size", color = "Community Size") +
  theme(legend.position = "inside", legend.position.inside = c(.99, .99), legend.justification = c("right", "top"), legend.box.just = "left") +
  scale_color_manual(NULL, values = c("#E69F00", "#0072B2"), labels = c("Large Community", "Small Community")) +
  scale_fill_manual(NULL, values = c("#E69F00", "#0072B2"), labels = c("Large Community", "Small Community"))

x_binned <- binData(counts = bs.test.large,
                    binWidth = "2k")

x_binned

x_binned$binVals <- x_binned$binVals %>%
  filter(binCount !=0)

x_binned$binVals

num_bins <- nrow(x_binned$binVals)

# bin breaks are the minima plus the max of the final bin:
bin_breaks <- c(dplyr::pull(x_binned$binVals, binMin),
                dplyr::pull(x_binned$binVals, binMax)[num_bins])

bin_counts <- dplyr::pull(x_binned$binVals, binCount)

mle_fish_bin <-  calcLike(negLL.PLB.binned,
                          p = -1.5,
                          w = bin_breaks,
                          d = bin_counts,
                          J = length(bin_counts),   # = num.bins
                          vecDiff = 1,
                          suppress.warnings = TRUE)             # increase this if hit a bound
mle_fish_bin

LBN_bin_plot(
  binValsTibble = x_binned$binVals,
  b.MLE = mle_fish_bin$MLE,
  b.confMin = mle_fish_bin$conf[1],
  b.confMax = mle_fish_bin$conf[2],
  leg.text = "(c)",
  log.xy = "xy")
