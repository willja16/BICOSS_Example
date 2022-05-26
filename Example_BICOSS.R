#########################################################################################################
# We present here a toy example with an analysis of a subset of the first dataset of the first simulation setting in Section 3.1 of 
# our manuscript.
# This example compares the performance of our proposed Bayesian Iterative Conditional Stochastic Search (BICOSS) method
# versus single marker analysis with Bonferroni multiplicity control (SMA-Bonf) and single marker analysis with Benjamini-Hochberg 
# multiplicity control (SMA-BH). In this toy example, BICOSS has better recall than the widely used SMA-Bonf method and 
# BICOSS shows much stronger FDR control vs the SMA methods.
# Results of our extensive simulation study can be found in Section 3.1 of our manuscript. Those results show that, 
# when compared to SMA-Bonf and SMA-BH, BICOSS has on average better recall and stronger FDR control.
#########################################################################################################

library(parallel)
library(caret)
library(doParallel)
library(GA)
library(memoise)
library(limma)
rm(list = ls())

## All files located at https://github.com/willja16/BICOSS_Example

set.seed(1330)

githubURL <- url("https://github.com/willja16/BICOSS_Example/raw/main/BICOSS_Example.RData")
load(githubURL)
close(githubURL)
rm(githubURL)

githubURL <- url("https://github.com/willja16/BICOSS_Example/raw/main/BICOSS.R")
source(githubURL)
close(githubURL)
rm(githubURL)

## Single Marker Association Test using SMA in GWAS.BAYES
# Returns a matrix where first column is 0,1 where 1 indicates significance according to threshold and multiple correction procedure
# Second column is p-values for each SNP
SMA_output <- SMA(Y = Y, SNPs = SNPs,kinship = kinship,number_cores = 4,threshold = 0.05,P3D = FALSE,selfing = TRUE,controlrate = "bonferroni")

## The BICOSS method is implemented with the BICOSS function.
## The parameters of the BICOSS function are as follows:
# Y:  continuous response vector
# SNPs:  matrix where each column is a single SNP coded as 0,1,2
# number_cores: number of cores to use for parallel computation
# threshold: type 1 error rate
# kinship: kinship matrix of size length(Y) by length(Y)
# maxiterations: maximum number of iterations of the genetic algorithm that performs model search
# runs_til_stop: number of consecutive iterations of the genetic algorithm at the best solution to declare convergence and stop
# P3D: logical parameter. If TRUE then BICOSS uses the concept of population parameters previously determined (P3D), that is, it estimates the parameter tau (Equation 1 of our manuscript) from the base model and uses that estimate to speed up computations in the screen step of the BICOSS algorithm. If P3D is FALSE then BICOSS estimates tau for each fitted model. 
# selfing: logical parameter. Use selfing=TRUE if the dataset is from a selfing species and use selfing=FALSE otherwise
# fixed_effects: matrix with length(Y) rows, with each column corresponding to a fixed effect that will be present in every fitted model. Currently BICOSS supports numeric columns of fixed effects. 
BICOSS_output <- BICOSS(Y = Y, SNPs = SNPs,number_cores = 4,threshold = 0.05,kinship = kinship,maxiterations = 400,runs_til_stop = 40,P3D = FALSE,selfing = TRUE)

True_causal_SNPs <- c(450,1350,2250,3150,4050,4950,5850,6750,7650,8550)

## Single marker analysis with Bonferroni multiplicity control (SMA-Bonf)
# SMA-Bonf: Total number of identified SNPs
sum(p.adjust(SMA_output$P_values,method = "bonferroni") < 0.05)
# SMA-Bonf: Recall
sum(which(p.adjust(SMA_output$P_values,method = "bonferroni") < 0.05) %in% True_causal_SNPs)/length(True_causal_SNPs)
# SMA-Bonf: FDR
sum(!(which(p.adjust(SMA_output$P_values,method = "bonferroni") < 0.05) %in% True_causal_SNPs))/sum(p.adjust(SMA_output$P_values,method = "bonferroni") < 0.05) 

## Single marker analysis with Benjamini-Hochberg multiplicity control (SMA - BH)
# SMA-BH: Total number of identified SNPs
sum(p.adjust(SMA_output$P_values,method = "BH") < 0.05)
# SMA-BH: Recall
sum(which(p.adjust(SMA_output$P_values,method = "BH") < 0.05) %in% True_causal_SNPs)/length(True_causal_SNPs)
# SMA-BH: FDR
sum(!(which(p.adjust(SMA_output$P_values,method = "BH") < 0.05) %in% True_causal_SNPs))/sum(p.adjust(SMA_output$P_values,method = "BH") < 0.05) 

## Bayesian Iterative Conditional Stochastic Search (BICOSS) 
# BICOSS: Total number of identified SNPs
length(BICOSS_output$BICOSS_SNPs)
# BICOSS: Recall
sum(BICOSS_output$BICOSS_SNPs %in% True_causal_SNPs)/length(True_causal_SNPs)
# BICOSS: FDR
sum(!(BICOSS_output$BICOSS_SNPs %in% True_causal_SNPs))/length(BICOSS_output$BICOSS_SNPs) 
