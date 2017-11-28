#!/usr/bin/R
#2017/11/27
# R CMD BATCH '--args [popN] [Pco] [Pmu] [K] [modN] [Pnew] [mR_mR_corr] [miR_miR_corr] [miR_mR_corr] [mR_mR_corr_pvalue] [miR_miR_corr_pvalue]
# [miR_mR_pvalue] [Output_Filename]' --no-save --no-restore Perform_GA.R
#
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
###read in correlation matrix
#miR_miR_corr
miR_miR_corr <- as.character(args[8])
#mR_mR_corr
mR_mR_corr <- as.character(args[9])
#miR_mR_corr
miR_mR_corr <-as.character(args[10])

###Read in p value matrix
#miR_miR_corr_pvalue
miR_miR_corr_pvalue <- as.character(args[11])
#mR_mR_corr_pvalue
mR_mR_corr_pvalue <- as.character(args[12])
#miR_mR_corr_pvalue
miR_mR_corr_pvalue <- as.character(args[13])
# File  for output
output <- as.character(args[14])

# Script for generating score matrix
#-----------------------------------------------------------------
# Functions
#-----------------------------------------------------------------

# Split strings by delimited to matrix
strsplit2 <- function (x, split, ...)
{
    x <- as.character(x)
    n <- length(x)
    s <- strsplit(x, split = split, ...)
    nc <- unlist(lapply(s, length))
    out <- matrix("", n, max(nc))
    for (i in 1:n) {
        if (nc[i])
            out[i, 1:nc[i]] <- s[[i]]
    }
    out
}

# Transform matrix with three columns to vector, the third column contains the value
# The first two column are paste by : to be names of the vector
unmatrix <- function (x)
{
    rnames <- rownames(x)
    cnames <- colnames(x)
    if (is.null(rnames))
        rnames <- paste("r", 1:nrow(x), sep = "")
    if (is.null(cnames))
        cnames <- paste("c", 1:ncol(x), sep = "")
    nmat <- outer(rnames, cnames, paste, sep = ":")
    vlist <- c(t(x))
    names(vlist) <- c(t(nmat))
    return(vlist)
}

# Transform list to a matrix with three columns
list2matrix <- function(list,names1,names2)
{
	list <- list[(list[,1] %in% names1) & (list[,2] %in% names2),]
	out <- matrix(0,nrow=length(names1),ncol=length(names2),dimnames=list(names1,names2))
	for(i in 1:dim(list)[1]){
		out[as.character(list[i,1]),as.character(list[i,2])] <- list[i,3]
	}
	out
}


#-----------------------------------------------------------------------------
# data process
#-----------------------------------------------------------------------------

# Read in corelation matrices
miR_miR_corr_raw <- as.matrix(read.table(miR_miR_corr,sep="\t",check.names=F))
mR_mR_corr_raw <- as.matrix(read.table(mR_mR_corr,sep="\t"))
miR_mR_corr_raw <- as.matrix(read.table(miR_mR_corr,sep="\t"))

# correct bugs which may existed in colnames
colnames(miR_miR_corr_raw) <- rownames(miR_miR_corr_raw)
colnames(mR_mR_corr_raw) <- rownames(mR_mR_corr_raw)
colnames(miR_mR_corr_raw) <- rownames(mR_mR_corr_raw)

# Read in p value matrices
miR_miR_corr_p <- as.matrix(read.table(miR_miR_corr_pvalue,sep="\t",check.names=F))
miR_mR_corr_p <- as.matrix(read.table(miR_mR_corr_pvalue,sep="\t"))
mR_mR_corr_p <- as.matrix(read.table(mR_mR_corr_pvalue,sep="\t"))

# Correct the scenario that different p values for the same pairs of PCC
miR_miR_corr_p[lower.tri(miR_miR_corr_p)] <- t(miR_miR_corr_p)[lower.tri(miR_miR_corr_p)]
mR_mR_corr_p[lower.tri(mR_mR_corr_p)] <- t(mR_mR_corr_p)[lower.tri(mR_mR_corr_p)]

# Set the diagonal of p value matrices to 1
diag(miR_miR_corr_p) <- 1
diag(mR_mR_corr_p) <- 1

# Read in target prediction of TFs & miRNAs
TF_gene_raw <- read.table("TF_gene.txt",sep=" ")
miR_gene_raw <- read.table("miR_gene.txt",sep=" ")

# Find the overlap genes among target prediction of miRNAs, target prediction of TFs and differential expressed genes
gene_list <- sort(intersect(intersect(unique(miR_gene_raw$V3),unique(rownames(mR_mR_corr_raw))),unique(TF_gene_raw$V3)))
TF_list <- sort(intersect(unique(TF_gene_raw$V2),gene_list))
miR_list <- sort(intersect(unique(miR_gene_raw$V2),unique(rownames(miR_mR_corr_raw))))
nTF_list <- sort(gene_list[-match(TF_list,gene_list)])
gene_list <- c(TF_list,nTF_list)

# Transform target prediction of TFs to symmetric matrix
TF_gene_inter <- list2matrix(TF_gene_raw[,c(2,3,1)],TF_list,gene_list)
# Apply similarity score cutoff for target prediction of TFs
TF_gene_inter[TF_gene_inter<0.99] <- 0

# Transform target prediction of miRNAs to symmetric matrix
miR_gene_raw[,1] <- as.numeric(miR_gene_raw[,1])/max(as.numeric(miR_gene_raw[,1]))
miR_gene_inter <- list2matrix(miR_gene_raw[,c(2,3,1)],miR_list,gene_list)

# Calculate the length of miRNAs, TFs and nTFs
miR_length <- length(miR_list)
TF_length <- length(TF_list)
nTF_length <- length(nTF_list)

N <- miR_length + TF_length + nTF_length
TF_s_index <- miR_length+1
nTF_s_index <- miR_length+TF_length+1
mR_length <- TF_length+nTF_length
rg_length <- miR_length+TF_length

# save the orignal correlation matrices
miR_miR_corr_b <- miR_miR_corr_raw
miR_mR_corr_b <- miR_mR_corr_raw
mR_mR_corr_b <- mR_mR_corr_raw

# Apply a cutoff to correlation coefficient
miR_miR_corr_raw[miR_miR_corr_p > 1e-03] <- 0
miR_mR_corr_raw[miR_mR_corr_p > 1e-03] <- 0
mR_mR_corr_raw[mR_mR_corr_p > 1e-03] <- 0

# Keep only the miRNAs/genes with predicted information
miR_miR_corr_new <- miR_miR_corr_raw[miR_list,miR_list]
miR_mR_corr_new <- miR_mR_corr_raw[miR_list,gene_list]
mR_mR_corr_new <- mR_mR_corr_raw[gene_list,gene_list]

#diag(miR_miR_corr_new) <- 0
#diag(mR_mR_corr_new) <- 0

# Any PCC between regulator and its target genes are kept regardless of the p values
miR_mR_corr_new[which(miR_gene_inter>0,arr.ind=T)] <- miR_mR_corr_b[which(miR_gene_inter>0,arr.ind=T)]
mR_mR_corr_new[which(TF_gene_inter>0,arr.ind=T)] <- mR_mR_corr_b[which(TF_gene_inter>0,arr.ind=T)]
mR_mR_corr_new[which(t(TF_gene_inter>0),arr.ind=T)] <- mR_mR_corr_b[which(t(TF_gene_inter>0),arr.ind=T)]

#miR_miR_corr <- as.matrix(miR_miR_corr_new[miR_list,miR_list])
miR_miR_corr <- as.matrix(miR_miR_corr_new[miR_list,miR_list])
miR_mR_corr <- as.matrix(miR_mR_corr_new[miR_list,gene_list])
miR_mR_corr[miR_list,TF_list] <- miR_mR_corr[miR_list,TF_list]

#######mR_mR_corr
mR_mR_corr <- as.matrix(mR_mR_corr_new[gene_list,gene_list])


miR_corrs <- cbind(miR_miR_corr,miR_mR_corr[miR_list,gene_list]/2)
# Use negative correlation only
#
# miR_corrs <- cbind(miR_miR_corr,miR_mR_corr[miR_list,gene_list]/2)
# miR_corrs[miR_corrs>0] <- 0
#

# Parameters for control the contibuition of sequence-based prediction(k1) and correaation(k2)
k1=1
k2=1

# Construct Socre matrix as below
#  -----------------------
# | miR-miR miR-TF miR-nTF|
# | TF-miR  TF-TF  TF-nTF |
# |	nTF-miR nTF-TF nTF-nTF|
#  -----------------------
Score_matrix <- as.matrix(rbind(k1*cbind(matrix(0,nrow=miR_length,ncol=miR_length), miR_gene_inter/2) + k2*abs(miR_corrs), cbind(k2*t(abs(miR_mR_corr[miR_list,gene_list]))/2 + k1*t(miR_gene_inter)/2, rbind(k1*cbind(TF_gene_inter[TF_list,TF_list],TF_gene_inter[TF_list,nTF_list]/2),k1*cbind(t(TF_gene_inter[TF_list,nTF_list])/2,matrix(0,nrow=nTF_length,ncol=nTF_length))) + k2*abs(mR_mR_corr))))

# Scale each part of score matrix to [0.5,1]
tmp <- Score_matrix[1:miR_length,1:miR_length]
Score_matrix[1:miR_length,1:miR_length] <- (1-0.5)*(tmp-min(tmp[tmp>0]))/(max(tmp)-min(tmp[tmp>0]))+0.5
tmp <- Score_matrix[TF_s_index:rg_length,1:miR_length]
Score_matrix[TF_s_index:TF_s_index:rg_length,1:miR_length] <- (1-0.5)*(tmp-min(tmp[tmp>0]))/(max(tmp)-min(tmp[tmp>0]))+0.5
tmp <- Score_matrix[1:miR_length,TF_s_index:rg_length]
Score_matrix[1:miR_length,TF_s_index:rg_length] <- (1-0.5)*(tmp-min(tmp[tmp>0]))/(max(tmp)-min(tmp[tmp>0]))+0.5
tmp <- Score_matrix[1:miR_length,nTF_s_index:N]
Score_matrix[1:miR_length,nTF_s_index:N] <- (1-0.5)*(tmp-min(tmp[tmp>0]))/(max(tmp)-min(tmp[tmp>0]))+0.5
tmp <- Score_matrix[nTF_s_index:N,1:miR_length]
Score_matrix[nTF_s_index:N,1:miR_length] <- (1-0.5)*(tmp-min(tmp[tmp>0]))/(max(tmp)-min(tmp[tmp>0]))+0.5

tmp <- Score_matrix[TF_s_index:rg_length,TF_s_index:rg_length]
Score_matrix[TF_s_index:rg_length,TF_s_index:rg_length] <- (1-0.5)*(tmp-min(tmp[tmp>0]))/(max(tmp)-min(tmp[tmp>0]))+0.5
tmp <- Score_matrix[TF_s_index:rg_length,nTF_s_index:N]
Score_matrix[TF_s_index:rg_length,nTF_s_index:N] <- (1-0.5)*(tmp-min(tmp[tmp>0]))/(max(tmp)-min(tmp[tmp>0]))+0.5
tmp <- Score_matrix[nTF_s_index:N,TF_s_index:rg_length]
Score_matrix[nTF_s_index:N,TF_s_index:rg_length] <- (1-0.5)*(tmp-min(tmp[tmp>0]))/(max(tmp)-min(tmp[tmp>0]))+0.5

tmp <- Score_matrix[nTF_s_index:N,nTF_s_index:N]
Score_matrix[nTF_s_index:N,nTF_s_index:N] <- (1-0.5)*(tmp-min(tmp[tmp>0]))/(max(tmp)-min(tmp[tmp>0]))+0.5

#####scorematrix
Score_matrix[Score_matrix < 0.5] <- 0
colnames(Score_matrix) <- rownames(Score_matrix)

# Generate information of category of names
# Example:
# Name 			Category
# hsa-miR-15a	miR
# MSIS1			TF
# HAND1			nTF
#


#####gene_category
geneCa <- rbind(cbind(miR_list,rep("miR",length(miR_list))),cbind(TF_list,rep("TF",length(TF_list))),cbind(nTF_list,rep("nTF",length(nTF_list))))
colnames(geneCa) <- c("Name","Category")


# Write score matrix, mRNA-mRNA PCC and category of names to files
#write.table(Score_matrix,"Score_matrix_even_inter.txt",sep="\t",row.names=T,col.names=T)
#write.table(geneCa,"gene_category.txt",sep="\t",row.names=T,col.names=T)
#write.table(mR_mR_corr,"mR_mR_corr.txt",sep="\t",row.names=T,col.names=T)

##########
#gene_category <- geneCa
#Score_matrix_even_inter <- Score_matrix
##########


# Perform module identification, result is outputed to R data object: **.RData
# Three file required: 	Score matrix
#						PCC matrix of mRNA-mRNA pairs
#						Category of gene names
#
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

####
#score_matrix
#Score_matrix_even_inter <- as.character(args[8])
#mR_mR_corr
#mR_mR_corr <- as.character(args[9])
#gene_category
#gene_category <- as.character(args[10])
####
########
# Options for control the size of regulators
#miR_per <- as.numeric(args[7])*100
#TF_per <- as.numeric(args[8])*100
# Initial population setup for GA
#nTF_per <- as.numeric(args[9])*100
##########











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
##########
#Score_matrix <- Score_matrix_even_inter
colnames(Score_matrix) <- rownames(Score_matrix)
mR_corr <- mR_mR_corr
colnames(mR_corr) <- rownames(mR_corr)
#geneCa <- gene_category
##########
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

save(list=ls(),file=paste(output,".RData",sep=""))


