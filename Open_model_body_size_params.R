# some parameters to think about in quote from Victor's paper:
# 16,113 individuals, distributed among 83 genera in Brazil,
# and 86,048 individuals distributed among 77 genera in Finland.
# The mean number of genera per stream was 17.84 ± 7.46 (mean ± SD) in Brazil and 14.01 ± 5.07 in Finland
# while the mean number of individuals per stream was 181.50 ± 111.38 and 886.57 ± 700.73

library(vegan)
library(sizeSpectra)
library(tidyverse)

setwd("G:/My Drive/Stream Metacommunities/Local Models")

#################################################
# population growth function
population.growth <- function(S, fec, rmax, alphas, N, Nall, z, niche, env) {
  
  surv <- S/(1+sum(alphas*Nall)) # survival function of individuals present during last time step

  surv.out <- rbinom(1,as.integer(N),surv) #demographic stochasticity in survival as binomial random deviates

  r <- max(0,rmax*exp(-((z-env)/(2*niche))^2))
  
  if (runif(1) < r) {
    r.input <- rpois(1,fec/(1+sum(alphas*Nall)))
  } else {
    r.input <- 0
  }
  
  N_next <- surv.out + r.input
  
  return(N_next)
}

#################################################
# simulation input parameters
# time to run the model
time <- 48
# number of species
species <- 150
# number of runs
runs <- 100

#################################################
# make storage
community.size   <- rep(NA, runs)
species.richness <- rep(NA, runs)
alpha.diversity  <- rep(NA, runs)
beta.diversity   <- rep(NA, runs)
ss.slope         <- rep(NA, runs)
ss.range         <- rep(NA, runs)

for (index in 1:runs){

#################################################
# autocorrelation of environmental noise
# mu = 0 is uncorrelated (white) noise
mu <- 0.75
# scale beta such that sigma has the same variance across all alpha values
beta <- (1-mu^2)^0.5
# mean of environmental conditions
env_mean <- 0.5
# variance in environmental conditions
zeta <- 0.05

#################################################
# population dynamics parameters (constants)
# max recruitment probability - same for all species 
#r.in <- 0.75
#rmax <- rep(r.in,species)
rmax <- runif(species,r.in*0.9,r.in+1.1)
rmax[rmax > 1] <- 1
rmax[rmax < 0] <- 0

# body mass dependent survival
min.mass <- 0.0002*1000 # convert to ug
max.mass <- 1300*1000 # convert to ug

mass <- 10^runif(species,log10(min.mass),log10(max.mass))
#mass <- runif(species,min.mass,max.mass)

b1 <- 7.5e+10
eV <- 0.63
k  <- 8.62e-05

mort <- b1*exp(-(eV/(k*temp)))*mass^(-0.25)

S <- sapply(exp(-mort), function(x) runif(1,x*0.9,x*1.1))

# fecundity
#fec.in <- 10
#fec <- rep(fec.in,species)
fec <- runif(species,fec.in*0.9,fec.in*1.1)
#fec <- -450*S + 295
#fec <- 0.28*(s^-2.95)

# place species on environmental axis
# niche optima
z <- runif(species, min = 0, max = 1) # pre-defined 0,1 axis

# niche breadths
#niche.in <- 0.1
#niche <- rep(niche.in, species)
niche <- runif(species,0.01,0.1)

# body mass dependent density dependence 
b2 <- 5.5e+5

alpha.in <- b2*exp(-(eV/(k*temp)))*mass^0.75;
# alphas <- matrix(alpha.in, ncol=species, nrow=species)
alphas <- matrix(data = runif(species^2,0,alpha.in/2), nrow = species, ncol = species) # check this!
diag(alphas) <- runif(species,alpha.in*0.9,alpha.in*1.1)

#################################################
# run specific storage
results_both <- matrix(NA, nrow=species, ncol=time)
results_both[,1] <- 10#(median(5)-1)/(max(diag(alphas))*species)

# create timeseries of environmental conditions
sigma <- rep(NA, time)
sigma[1] <- 0
for (t in 1:(time-1)) {
  sigma[t+1] <- mu*sigma[t]+beta*rnorm(1,0,1)
}
env <- env_mean+zeta*sigma #Does it make sense to bound?

# run population dynamics
for (t in 1:(time-1)) {
  for(k in 1:species) {
    results_both[k,t+1] <- population.growth(S[k], fec = fec[k], rmax = rmax[k], alpha=alphas[k,], N=results_both[k,t],
                                          Nall=results_both[,t], z=z[k], niche=niche[k], env=env[t])
  }
}

#average community size
community.size[index] <- exp(mean(log(colSums(results_both))))

#average species richness
species.richness[index] <- median(colSums(results_both>=1))

#average species diversity 
alpha.diversity[index] <- median(diversity(results_both, index = "shannon", MARGIN = 2))

#average beta diversity
beta.store <- as.matrix(vegdist(t(results_both[,]), method = "bray", na.rm = TRUE))
beta.diversity[index] <- median(diag(beta.store[-nrow(beta.store),-1]))

#median size spectra slope
source("size_spectra_analysis_bin_out.R")
ss.slope[index] <- ss.mle.median
ss.range[index] <- ss.mle.range

}

#################################################
# Data Summary
# median.community.size <- median(community.size)
# upper.quantile.comm   <- quantile(community.size,0.975)
# lower.quantile.comm   <- quantile(community.size,0.025)
# 
# median.richness         <- median(species.richness)
# upper.quantile.richness <- quantile(species.richness,0.975)
# lower.quantile.richness <- quantile(species.richness,0.025)
# 
# median.alpha.diversity <- median(alpha.diversity)
# upper.quantile.alpha   <- quantile(alpha.diversity,0.975)
# lower.quantile.alpha   <- quantile(alpha.diversity,0.025)
# 
# median.beta.diversity  <- median(beta.diversity,na.rm = "TRUE")
# upper.quantile.beta    <- quantile(beta.diversity,0.975,na.rm = "TRUE")
# lower.quantile.beta    <- quantile(beta.diversity,0.025,na.rm = "TRUE")

# data.output <- data.frame(time, species, mu, zeta, r.in, S.in, fec.in, niche.in, alpha.in,
#                         median.community.size, upper.quantile.comm, lower.quantile.comm,
#                         median.richness, upper.quantile.richness, lower.quantile.richness,
#                         median.alpha.diversity, upper.quantile.alpha, lower.quantile.alpha,
#                         median.beta.diversity, upper.quantile.beta, lower.quantile.beta, row.names=NULL)
                 
#write data to files
#file.str <- paste("data",species,temp,mu,zeta,r.in,S,fec.in,niche.in,alpha.in,".csv", sep = "_")

#write.csv(data.output, file = file.str, row.names = FALSE)

data.output <- data.frame(r.in,temp,community.size,species.richness,alpha.diversity,beta.diversity,ss.slope,ss.range)

file.str <- paste("bs_temp_runs/data_",r.in,"_",temp,".csv", sep = "")

write.csv(data.output, file = file.str, row.names = FALSE)

#################################################
# plots

 #plot(community.size,beta.diversity)
# 
# plot niche axes
#  niche_seq <- seq(0,1,by=0.01)
#  color_v <- rainbow(species)
# #
# plot(1, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0, 1))
# 
# for (index in 1:species)
# {
#    r_range <- rmax[index]*exp(-((z[index]-niche_seq)/(2*niche[index]))^2)
#    lines(niche_seq, r_range, type="l", col = color_v[index])
# }

#plot(1:time, env, type = "l", lty = 1, ylim = c(0,1))
#
#plot(1:time,colSums(results_both), type = "l", lty = 1)
#
#matplot(1:time, cbind(colSums(results_both),t(results_both)), type="l", lwd = 2, col = 1:species, xlab="Time", ylab="Density", bty="n", xlim = c(1,time), log = "y")
#  
# hist(species.richness, xlab="Species Richness", main="summary across all runs")
#  
# hist(alpha.diversity, xlab="Shannon Diversity", main="summary across all runs")
#  
# hist(beta.diversity, xlab="Bray-Curtis Dissimilarity", main="summary across all runs")
#  
# #print(data.output)