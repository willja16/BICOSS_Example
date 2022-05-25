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

## All files located at ...

githubURL <- "https://github.com/willja16/BICOSS_Example/raw/main/BICOSS_Example.RData"
load(url(githubURL))
rm(githubURL)

githubURL <- "https://github.com/willja16/BICOSS_Example/raw/main/BICOSS.R"
source(url(githubURL))
rm(githubURL)

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
example_output <- BICOSS(Y = Y, SNPs = SNPs,number_cores = 24,threshold = 0.05,kinship = kinship,maxiterations = 400,runs_til_stop = 40,P3D = FALSE,selfing = TRUE)

### Significant SNPs: 450,1350,2250,3150,4050,4950,5850,6750,7650,8550

## SMA - Bonf
which(p.adjust(example_output$p_values[[1]]$P_values,method = "bonferroni") < 0.05)
## SMA - BH
which(p.adjust(example_output$p_values[[1]]$P_values,method = "BH") < 0.05)

## Best BICOSS Model
example_output$modelselection[[1]]$Models
