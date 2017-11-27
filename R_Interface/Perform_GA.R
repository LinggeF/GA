#!/usr/bin/R
# Perform module identification, result is outputed to R data object: **.RData
# Three file required: 	Score matrix: Score_matrix_even_inter.txt 
#						PCC matrix of mRNA-mRNA pairs: mR_mR_corr.txt
#						Category of gene names: gene_category.txt
#
# Created by Wenbo Mu, May 15,2012
#
# R CMD BATCH '--args [popN] [Pco] [Pmu] [K] [modN] [Pnew] [miR_per] [TF_per] [nTF_per] [Output_Filename]' --no-save --no-restore Perform_GA.R
# 
# Suggested parameter setup: popN (population size): 100;
#							 Pco (crossover probability): 0.7;
#							 Pmu (mutation probability): 0.001;
#							 modN (# of modules generated): 10;
#							 K (GA iteration times): 5000;
#							 LK (Local search iteration times): 1000;
#							 Pnew (probability of random immigrant): 0.01
#							 miR_per (Option of prefered miRNA percentage): 0.02
#							 TF_per (Option of prefered TF-gene percentage): 0.04
#							 nTF_per (Option of prefered nTF-gene percentage): 0.02
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# initialize parameter
#-----------------------------------------------------------------------------

args <- commandArgs(TRUE);
# population size
popN <- as.numeric(args[1])
# crossover probability
Pco <- as.numeric(args[2])
# putation probability
Pmu <- as.numeric(args[3])
# iteration times
K <- as.numeric(args[4])
# number of modules
modN <- as.numeric(args[5])
#probability of introducing fresh chromesomes
Pnew <- as.numeric(args[6])
# iteration times for Local search
LK <- as.numeric(args[7])
#score_matrix
Score_matrix_even_inter <- as.character(args[8])
#mR_mR_corr
mR_mR_corr <- as.character(args[9])
#gene_category
gene_category <- as.character(args[10])
#####
# Options for control the size of regulators
#miR_per <- as.numeric(args[7])*100
#TF_per <- as.numeric(args[8])*100
# Initial population setup for GA
#nTF_per <- as.numeric(args[9])*100
#####




# File  for output
output <- as.character(args[11])






# loading GA functions
source("GA_1.2.R")

# set parameter
#popN <- 100
#Pco <- 0.7
#Pmu <- 0.001
#K <- 500
#modN <- 10
#Pnew <- 0.01
#LK <-1000







# Options for control the size of regulators
TF_per <- 0.04*100
miR_per <- 0.02*100
# Initial population setup for GA
nTF_per <- 0.02*100

#-----------------------------Data Import-----------------------------------------------
# Loading Score Matrix
# Loading PCC matrix
# Calculate # of miRNAs/TFs/nTFs

Score_matrix <- read.table(Score_matrix_even_inter,sep="\t",header=T)
colnames(Score_matrix) <- rownames(Score_matrix)
mR_corr <- read.table(mR_mR_corr,sep="\t",header=T)
colnames(mR_corr) <- rownames(mR_corr)

geneCa <- read.table(gene_category,sep="\t",header=T)
miR_list <- as.character(geneCa[which(geneCa[,"Category"]=="miR"),"Name"])
TF_list <- as.character(geneCa[which(geneCa[,"Category"]=="TF"),"Name"])
nTF_list <- as.character(geneCa[which(geneCa[,"Category"]=="nTF"),"Name"])
miR_length <- length(miR_list)
TF_length <- length(TF_list)
nTF_length <- length(nTF_list)

N <- miR_length + TF_length + nTF_length
TF_s_index <- miR_length+1
nTF_s_index <- miR_length+TF_length+1
mR_length <- TF_length+nTF_length
rg_length <- miR_length+TF_length

#-----------------------------------------------------------------------------
# Formal Run of Algorithm
#-----------------------------------------------------------------------------

# Multiple core setup
registerDoMC(cores=6)
final <- Module_Identification(Pmu,Pco,popN,modN,K,LK,Pnew,verb=TRUE)


os = file(output,"w")


lapply(final,write,os, append =F)
#write(final,os,append = F)

save(list=ls(),file=paste(output,".RData",sep=""))


