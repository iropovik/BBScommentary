#' ---
#' title: "The lack of robust evidence for the cleansing effect"
#' author: "Ivan Ropovik"
#' date: "`r Sys.Date()`"
#' output:
#'    html_document:
#'       toc: true
#'       toc_float: true
#'       code_folding: show
#'       fig_retina: 2
#' always_allow_html: yes
#' ---

# Install required R libraries if not installed already
list.of.packages <- c("tidyverse", "moments")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load required libraries
#+ include = FALSE
lapply(list.of.packages, require, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)

# Source the script for the p-curve method
source("p-curve.R")

# Read in the data
data <- read_csv("data/data.csv")

# Median Ns ---------------------------------------------------------------

#'#### Median Ns
data %>% group_by(replicated) %>% summarise(meanN = mean(N, na.rm = TRUE),
                                            medianN = median(N, na.rm = TRUE)) %>% na.omit()

# Simulation of p-values --------------------------------------------------

#'# Data Simulation
#'
#'## If H0 is false
pMedianSample <- median(data$pReported[data$pReported <= .051], na.rm = TRUE)
pSkewSample <- skewness(data$pReported[data$pReported <= .051], na.rm = TRUE)
pVarianceSample <- var(data$pReported[data$pReported <= .051], na.rm = TRUE)
pMedians <- list()
pSkew <- list()
pVariance <- list()
nSim <- 10000
set.seed(1)

for(j in 1:nSim){
sigma = 1
designMatrix <- expand.grid(
  "n" = data[!is.na(data$N),]$N,
  "es" = seq(from = .1, to = .6, by = .1)
)

pvals <- list()
for(i in 1:nrow(designMatrix)){
  tstats = replicate(22, t.test(rnorm(designMatrix$n[i], sigma*designMatrix$es[i], sigma))$statistic)
  pvals[[i]] = 1 - pf(tstats^2, df1 = 1, df2 = designMatrix$n[i] - 1)
  }

ps <- pvals %>% map(.,~median(.[sum(. < .05) >= 7][. < .05])) %>% unlist()
skew <- pvals %>% map(.,~skewness(.[sum(. < .05) >= 7][. < .05])) %>% unlist()
variance <- pvals %>% map(.,~var(.[sum(. < .05) == 7][. < .05])) %>% unlist()
pMedians[[j]] <- table(ps >= pMedianSample)[2]
pSkew[[j]] <- table(skew <= pSkewSample)[2]
pVariance[[j]] <- table(variance <= pVarianceSample)[2]
}

#'#### Cumulative probability of observing seven or more significant effects for which the median p-value was the same or higher than the median of the observed p-value distribution
sum(unlist(pMedians), na.rm = TRUE)/(nSim*nrow(designMatrix))

#'#### Cumulative probability of observing seven or more significant effects for which the skewness was the same or or more negative than the skewness of the observed p-value distribution
sum(unlist(pSkew), na.rm = TRUE)/(nSim*nrow(designMatrix))

#'#### Cumulative probability of observing seven significant effects for which the variance of p-values was the same or smaller than the variance in the observed p-value distribution
sum(unlist(pVariance), na.rm = TRUE)/(nSim*nrow(designMatrix))

#'## If H0 is true
nSimH0 <- 1e7
ps <- list()
pMediansH0 <- list()
pSkewH0 <- list()
pVarianceH0 <- list()
set.seed(1)
for(i in 1:nSimH0){
  ps[[i]] <- runif(22)
}

pMediansH0 <- ps %>% map(.,~median(.[sum(. < .05) >= 7][. < .05])) %>% unlist()
pSkewH0 <- ps %>% map(.,~skewness(.[sum(. < .05) >= 7][. < .05])) %>% unlist()
pVarianceH0 <- ps %>% map(.,~var(.[sum(. < .05) == 7][. < .05])) %>% unlist()

#'#### Cumulative probability of observing seven or more significant effects for which the median p-value was the same or higher than the median of the observed p-value distribution
table(pMediansH0 >= pMedianSample)[2]/nSimH0

#'#### Cumulative probability of observing seven or more significant effects for which the skewness was the same or or more negative than the skewness of the observed p-value distribution
table(pSkewH0 <= pSkewSample)[2]/nSimH0

#'#### Cumulative probability of observing seven significant effects for which the variance of p-values was the same or smaller than the variance in the observed p-value distribution
table(pVarianceH0 <= pVarianceSample)[2]/nSimH0

# p-curve analysis --------------------------------------------------------

# Saves the output as a figure in the working directory
data$pReported[19] <- .049999
data$zvalue <- NA
data$zvalue <- qnorm(data$pReported/2)
zvals <- paste("Z=", data$zvalue[!is.na(data$zvalue)], sep = "")
pcurve_app(zvals)



