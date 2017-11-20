#!/usr/bin/R
#
# Functions for module identification
#
# Created by Wenbo Mu
# Aug 15, 2012
# 
#
#-----------------------------------------------------------------
# Functions
#-----------------------------------------------------------------

# Import library for parallele computing
library("foreach")
library("doMC")

# Fitness function for GA
cal_score_GA <- function(mat,score_sub)
{
#t0 <- Sys.time()
	nNode <- length(which(mat))
	if(nNode == 0 || nNode == 1)
	{
		score <- -Inf	
	}
	else
	{
		score_m <- score_sub[mat,mat]
		score <- sum(score_m)/length(which(score_m>0))
	}
	score
}

# Fitness measurement of local search

cal_score_LS <- function(mat,score_sub)
{
	nNode <- length(which(mat[1:rg_length]))
	# In case no node is included in solution
	if(nNode == 0 || nNode == 1)
	{
		score <- -Inf	
	}
	else
	{	
		score_1 <- score_sub[which(mat[1:rg_length]),mat]
		score_2 <- score_1[,which(mat)>rg_length]
		score <- (sum(score_1)+sum(score_2))/(2*(length(which(score_1>0))+length(which(score_2>0))))
	}
	score
}


# Obtain indexes for parents selection
obtain_index <- function(x)
{
	cum <- cumsum(x)
	randn <- runif(2)
	cum_wr <- append(cum,randn)
	rrank <- rank(cum_wr)
	index1 <- rrank[length(x)+1]
	index2 <- rrank[length(x)+2]
	r <- list(index1=index1,index2=index2)
}

# Obtain index of crossover location
obtain_cvid <- function(x,randn)
{
	cum_wr <- c(x,randn)
	rrank <- ceiling(rank(cum_wr)[length(x)+1])
}

# Shuffule matrix to generate random dataset for permutation test
matrix_shuffle <- function(data)
{
    o_col <- colnames(data)
    o_row <- rownames(data)
    new_data <-	data
    rownames(new_data) <- sample(o_row)
    colnames(new_data) <- sample(o_col)
    new_data <- new_data[o_row,o_col]
}

# Count # of miRNAs/TFs/nTFs in solution
count_category <- function(x)
{
	index <- which(x)
	miR_count <- length(which(index<=miR_length))
	nTF_count <- length(which(index > rg_length))
	TF_count <- length(index)-miR_count-nTF_count
	c(miR_count,TF_count,nTF_count)
}

# Repair GA results by removing isolated genes
module_cleanup <- function(x,score_m)
{
	mod_m <- score_m[x,x]
	indexes <- which(x)
	for(i in 1:length(indexes))
	{
		if(length(which(mod_m[i,]!=0)) < 1)
		{
			x[indexes[i]] <- FALSE
		}
	}
	x
}

# Genetic Algorithm for co-expression gene set identification
GA <- function(Pmu,Pco,popN,modN,K,Pnew,score_m,verb=TRUE,out_detail=FALSE)
{
	# Generate random population based on pre-decided percentage of TF/nTF 
	pop <- foreach(i = 1:mR_length, .combine="rbind", .multicombine=TRUE, .inorder=F) %do% c(sample(c(rep(TRUE,TF_per),rep(FALSE,100-TF_per)),TF_length,replace=T),sample(c(rep(TRUE,nTF_per),rep(FALSE,100-nTF_per)),nTF_length,replace=T))
	score <- foreach(i = 1:popN, .combine="c", .multicombine=TRUE, .inorder=T) %dopar% cal_score_GA(pop[i,],score_m)
	pop_o <- pop[order(score,decreasing=T),]
	score_weight <- (score[order(score,decreasing=T)]-min(score))/sum(score-min(score))
	# record scores for each iteration
	idv_record <- vector()
	pop_new <- matrix(logical(),nrow=popN,ncol=mR_length)
	print("Run Genetic Algorithm for mRNAs")	
	for (iter in 1:K)
		{		
		# iteration report
			if(iter%%100==1 & verb==TRUE)
			{
				print(paste("iteration",iter))
			}
				
			pop_new[c(1,2),] <- pop_o[1:2,]
			for(Lin in seq(3,popN,by=2))
				{
					# parent choosen
					indexes <- obtain_index(score_weight)
					parent1_id <- indexes$index1
					parent2_id <- indexes$index2

					chrom1 <- pop_o[parent1_id,]
					chrom2 <- pop_o[parent2_id,]
					# Crossover
					if(runif(1) < Pco)
					{
						crossover_id <- c(sample(1:(TF_length-1),1),sample(1:(nTF_length-1),1))

						TF_chrom1 <- chrom1[1:TF_length]
						nTF_chrom1 <- chrom1[(TF_length+1):mR_length]
						TF_chrom2 <- chrom2[1:TF_length]
						nTF_chrom2 <- chrom2[(TF_length+1):mR_length]
					
						offspring1 <- c(TF_chrom1[1:crossover_id[1]],TF_chrom2[(crossover_id[1]+1):TF_length],nTF_chrom1[1:crossover_id[2]],nTF_chrom2[(crossover_id[2]+1):nTF_length])
						offspring2 <- c(TF_chrom2[1:crossover_id[1]],TF_chrom1[(crossover_id[1]+1):TF_length],nTF_chrom2[1:crossover_id[2]],nTF_chrom1[(crossover_id[2]+1):nTF_length])					
					}				
					else
					{
						offspring1 <- chrom1
						offspring2 <- chrom2
					}
					# Mutation
					Pmu_select1 <- ifelse(runif(mR_length,0,1) < Pmu, TRUE, FALSE)
					Pmu_select2 <- ifelse(runif(mR_length,0,1) < Pmu, TRUE, FALSE)
					offspring1[Pmu_select1] <- !offspring1[Pmu_select1]
					offspring2[Pmu_select2] <- !offspring2[Pmu_select2]	

					pop_new[Lin,] <- offspring1
					pop_new[Lin+1,] <- offspring2
				}
				
			score <- foreach(i = 1:popN, .combine="c", .multicombine=TRUE, .inorder=T) %dopar% cal_score_GA(pop_new[i,],score_m)
			score_order <- order(score,decreasing=T)

			# Random immigrant
			if(runif(1,0,1) < Pnew)
			{
				pop_o[which(score==min(score)),] <- c(sample(c(rep(TRUE,TF_per),rep(FALSE,100-TF_per)),TF_length,replace=T),sample(c(rep(TRUE,nTF_per),rep(FALSE,100-nTF_per)),nTF_length,replace=T))
				score[which(score==min(score))] <- cal_score_GA(pop_o[popN,],score_m)
				score_order <- order(score,decreasing=T)
			}
			pop_o <- pop_new[score_order,]
			score_weight <- (score[score_order]-min(score))/sum(score-min(score))	
			idv_record[iter] <- score[score_order[1]]
			# Termination condition, Score not improved for 200 times
			if( iter>200 && idv_record[iter] == idv_record[iter-200])
			{
				break
			}
		
		} # iteration end
	# If detail information are needed, output the original population and score track
	# Otherwise only output solutions	
	if(out_detail==TRUE)
	{
		out <- list(pop=pop_o,record=idv_record)
	}
	else
	{
		pop_o
	}
}

# Local search for detecting best regulators
Local_Search <- function(x,LK,score_m,verb=TRUE,out_detail=FALSE)
{
	print("Local search to improve solutions")
	# Construct searching space
	# Any direct regulator are included
	tmp <- score_m[which(x)+miR_length,1:rg_length]
	cand_rg_du <- which(x[1:TF_length])+miR_length
	for(i in 1:dim(tmp)[1])
	{
		cand_rg_du <- c(cand_rg_du,which(tmp[i,] > 0))
	}
	cand_rg <- sort(unique(cand_rg_du))
	miR_rg <- cand_rg[cand_rg<= miR_length] 
	TF_rg <- cand_rg[cand_rg > miR_length]
	# Construct initial solution
	best <- c(rep(FALSE,miR_length),x)
	best[sample(cand_rg[cand_rg<=miR_length],min((miR_per*miR_length)/100),length(cand_rg))] <- TRUE
	best_score <- cal_score_LS(best,score_m)
	start <- which(best[1:rg_length]) 
	# Record scores
	idv_record <- vector()
	new_chr <- best
	for(iter in 1:LK)
	{
			# Report # of iterations
			if(iter%%100==1 & verb==TRUE)
			{
				print(paste("iteration",iter))
			}
			new_chr <- best
			cand_exist <- which(best)
			# miRNA and TF are chosen to change based on a probability according to their percentage in solution
			if(runif(1) <=miR_length/(miR_length+TF_length))
			{
				# Find miRNAs already included
				miR_in <- cand_exist[which(cand_exist<=miR_length)]
				# If # of miRNA < threshold, the bit at chosen location switch to the opposite
				if(length(miR_in) < miR_length*miR_per/100)
				{
					loc <- sample(miR_rg,1)
					new_chr[loc] <- !new_chr[loc]
				}
				# if # of miRNA > threshold, remove one node by chosing a bit from current solution and set false
				else if (length(miR_in) > miR_length*miR_per/100)
				{
					loc <- sample(miR_rg[!(miR_rg %in% miR_in)],1)
					new_chr[loc] <- FALSE
				}
				# if # of miRNA == threshold, substitute one node by another
				else
				{
					loc <- sample(miR_rg[!(miR_rg %in% miR_in)],1)
					loc_rm <- sample(miR_in,1)
					new_chr[loc_rm] <- FALSE
					new_chr[loc] <- TRUE
				}			
			}
			else
			{
				# Find TFs already included
				TF_in <- cand_exist[which(cand_exist >miR_length & cand_exist <= rg_length)]
				# If # of TF < threshold, the bit at chosen location switch to the opposite
				if(length(TF_in) < TF_length*TF_per/100)
				{
					#loc <- sample(TF_rg[!(TF_rg %in% TF_in)],1)
					loc <- sample(TF_rg,1)
					new_chr[loc] <- !new_chr[loc]
					#new_chr[loc] <- TRUE
				}
				# if # of TF > threshold, remove one node by chosing a bit from current solution and set false
				else if (length(TF_in) > TF_length*TF_per/100)
				{
					loc <- sample(TF_rg[!(TF_rg %in% TF_in)],1)
					new_chr[loc] <- FALSE
				}
				# if # of TF == threshold, substitute one node by another			
				else
				{
					loc <- sample(TF_rg[!(TF_rg %in% TF_in)],1)
					loc_rm <- sample(TF_in,1)
					new_chr[loc_rm] <- FALSE
					new_chr[loc] <- TRUE
				}			
			}
			new_score <- cal_score_LS(new_chr,score_m) 
			if( new_score > best_score)
			{
				best <- new_chr
				best_score <- new_score
			}
			idv_record[iter] <- best_score
			# Terminate if score not improved for 100 times
			if( iter>100 && idv_record[iter] <= idv_record[iter-100])
			{
				break
			}

	} 
	# If detail information needed, the best solution, score track, start solution and end solution are outputed
	if(out_detail==TRUE)
	{
		out <- list(best=best,record=idv_record,start=start,end=which(best[1:rg_length]))
	}
	else
	{
		best
	}
}


#-----------------------------------------------------------------
# Primary Algorithm
# Intermedia result of GA is also outputed
#-----------------------------------------------------------------
Module_Identification <- function(Pmu,Pco,popN,modN,K,LK,Pnew,verb=TRUE)
{
	print("Algorithm Start")
	Modules <- matrix(logical(),nrow=modN,ncol=N)
	# Matrix for store GA solutions 
	ga_result_out <- matrix(logical(),nrow=modN,ncol=mR_length)
	m <- 1
	score_sub <- as.matrix(Score_matrix)
	corr <- abs(as.matrix(mR_corr))
	while(m <= modN)
	{
		print(paste("Module",m))
		# Perform GA to identify co-expressed gene set
		ga_result <- GA(Pmu,Pco,popN,modN,K,Pnew,corr,verb,FALSE)		
		# Repair GA solution
		best_chr <- module_cleanup(ga_result[1,],corr)
		# In case there no nTF genes exists, regenerate new solutions
		if(length(which(best_chr[(TF_length+1):mR_length])) < 3)
		{
			print("No nTF genes exists. Regenerate.")
		}
		else
		{
			ga_result_out[m,] <- best_chr
			# Perform LS to detect best regulators
			solution <- Local_Search(best_chr,LK,score_sub,verb,FALSE)	
			Modules[m,] <- solution
			# Update correlation matrix in case duplicated gene sets are identified by GA
			corr[best_chr,best_chr] <- matrix(0,nrow=length(which(best_chr)),ncol=length(which(best_chr)))	
			m <- m+1		
		}		
	} # module end
	result <- list(module=Modules,ga=ga_result_out)
} # GA end

		
	