library(readr)
library(readxl)
library(readr)
library(dplyr)
library(lme4)
library(parallel)
library(qqman)
library(rrBLUP)
library(caret)
library(doParallel)
library(MASS)
library(limma)
rm(list = ls())

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

## The parameters of the BICOSS are as follows
# Y is the continous response vector 
# SNPs is a matrix where each column is a single SNP coded as 0,1,2
# number_cores sets the number of cores to use for this procedure
# threshold is the type 1 error rate
# kinship is the kinship matrix of size n x n 
# maxiterations is the maximum iterations to allow the genetic algorithm to iterate for
# runs_til_stop is the number of consecutive iterations at the best solution of the genetic algorithm upon which to stop
# P3D stands for population parameters previously determined, this estimates tau from the base model and uses that estimate to speed up computations
# selfing is a TRUE/FALSE value, where TRUE indicates SNPs from a selfing species and FALSE indicates SNPs from a non-selfing species
# fixed_effects is a matrix with n rows, each column should be a fixed effect that will be present in every model, currently BICOSS supports numeric columns of fixed effects
#  so if one has categorical fixed effects use model.matrix() and upload the design matrix without the intercept
BICOSS_output <- BICOSS(Y = Y, SNPs = SNPs,number_cores = 4,threshold = 0.05,kinship = kinship,maxiterations = 400,runs_til_stop = 40,P3D = FALSE,selfing = TRUE)

True_causal_SNPs <- c(450,1350,2250,3150,4050,4950,5850,6750,7650,8550)

## SMA - Bonf
# SMA - Bonf Total number of SNPs
sum(p.adjust(SMA_output$P_values,method = "bonferroni") < 0.05)
# SMA - Bonf Recall
sum(which(p.adjust(SMA_output$P_values,method = "bonferroni") < 0.05) %in% True_causal_SNPs)/length(True_causal_SNPs)
# SMA - Bonf FDR
sum(!(which(p.adjust(SMA_output$P_values,method = "bonferroni") < 0.05) %in% True_causal_SNPs))/sum(p.adjust(SMA_output$P_values,method = "bonferroni") < 0.05)

## SMA - BH
# SMA - BH Total number of SNPs
sum(p.adjust(SMA_output$P_values,method = "BH") < 0.05)
# SMA - BH Recall
sum(which(p.adjust(SMA_output$P_values,method = "BH") < 0.05) %in% True_causal_SNPs)/length(True_causal_SNPs)
# SMA - BH FDR
sum(!(which(p.adjust(SMA_output$P_values,method = "BH") < 0.05) %in% True_causal_SNPs))/sum(p.adjust(SMA_output$P_values,method = "BH") < 0.05)

## Best BICOSS Model
# BICOSS Total number of SNPs
length(BICOSS_output$BICOSS_SNPs)
# BICOSS Recall
sum(BICOSS_output$BICOSS_SNPs %in% True_causal_SNPs)/length(True_causal_SNPs)
# BICOSS FDR
sum(!(BICOSS_output$BICOSS_SNPs %in% True_causal_SNPs))/length(BICOSS_output$BICOSS_SNPs)

##################
# This toy example is a compressed verions of the first dataset of the first simulation setting.
# BICOSS showed strong FDR control vs the SMA methods.
# The recall was slightly lower compared to SMA-BH but upon further investigation BICOSS identified SNP 3148 instead of SNP 3150. A non-negigible difference in practice.
