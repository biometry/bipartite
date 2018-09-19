## node overlap and segretation (Strona & Veech 2015) 
## LaTeX symbol: \m{J} 

NOS <- function(web, keep.Nij=FALSE, keep.diag=FALSE){
	# N_ij must be computed once for each level: Nbar_in and Nbar_out
	# keep.diag: I think this should be excluded, but the paper is not clear about it
	# returns a list with various objects

	# S_ij : number of neighbours (interacting partners) shared by i and j
	# P_ij : expected number of shared neighbours
	# d_.  : node degrees of i and j
	# omega_ij: maximum possible value of N_ij (for standardisation)

	# Nbar : mean of all N_ij
	# sd of N_ij : measure of modularity
	# n : number of nodes that can be potentially shared by i and j

	# p_ij is the expected number of shared for nodes i and j, for each value of k=1 to min(d_i, d_j)
	
	out <- list()
	
	# for the higher trophic level, Nbar_in:
	n_in <- ncol(web) # = number of higher-trophic level species!
	# for the lower trophic level, Nbar_out:
	n_out <- nrow(web)
	
	degrees_lower  <- rowSums(web > 0) # degree of species i
	degrees_higher <- colSums(web > 0)
	
	# NOS is computed for each possible species pair. 
	## For the higher trophic level:
	N_ij_in <- matrix(0, n_in, n_in)
	for (i in 1:n_in){
		for (j in i:n_in){ # compute only upper triangle of the matrix!
			# compute actually shared species, S_ij:
			S_ij <- sum( rowSums((web[,c(i,j )]>0)) == 2) 
			# compute expected number of shared species, P_ij:
			p_ij <- 0
			d_i <- degrees_higher[i]
			d_j <- degrees_higher[j]
			n <- n_out # number of possible interactors, i.e. number of species in the other trophic level
			for (k in 1:min(d_i, d_j)){
				p_ij[k] <- (choose(n, k) *  choose(n - k, d_j - k) * choose(n - d_j, d_i - k)) / (choose(n, d_j) * choose(n, d_i)) * k
			}
			P_ij <- sum(p_ij)
			#Compute N_ij from (Strona & Veech, Eq. 6).  
			#Instead of computing Omega_ij, simplify Eq. 6 and compute N_ij directly, to avoid division by 0 when min(d_i,d_j)==0
			if (S_ij == P_ij){
			  N_ij_in[i,j] <- 0
			}
			else if (S_ij > P_ij){ 
			  N_ij_in[i,j] <- (S_ij - P_ij)/(min(d_i,d_j) - P_ij)
			}
			else{
			  if ((d_i + d_j - n) < 0){
			    N_ij_in[i,j] <- (S_ij - P_ij)/P_ij
			  }
			  else{
			    N_ij_in[i,j] <- (S_ij - P_ij)/(P_ij - (d_i + d_j - n))
			  }
			}
		}
	}
	#if (keep.Nij) out$"N_ij_in" <- N_ij_in
	#out$Nbar_in <- mean(N_ij_in)
	
	## same for lower trophic level:
	N_ij_out <- matrix(0, n_out, n_out)
	for (i in 1:n_out){
		for (j in i:n_out){ # compute only upper triangle of the matrix!
			# compute actually shared species, S_ij:
			S_ij <- sum( colSums((web[c(i,j ), ]>0)) == 2) 
			# compute expected number of shared species, P_ij:
			p_ij <- 0
			d_i <- degrees_lower[i]
			d_j <- degrees_lower[j]
			n <- n_in
			for (k in 1:min(d_i, d_j)){
				p_ij[k] <- (choose(n, k) *  choose(n - k, d_j - k) * choose(n - d_j, d_i - k)) / (choose(n, d_j) * choose(n, d_i)) * k
			}
			P_ij <- sum(p_ij)
			if (S_ij == P_ij){
			  N_ij_out[i,j] <- 0
			}
			else if (S_ij > P_ij){ 
			  N_ij_out[i,j] <- (S_ij - P_ij)/(min(d_i,d_j) - P_ij)
			}
			else{
			  if ((d_i + d_j - n) < 0){
			    N_ij_out[i,j] <- (S_ij - P_ij)/P_ij
			  }
			  else{
			    N_ij_out[i,j] <- (S_ij - P_ij)/(P_ij - (d_i + d_j - n))
			  }
			}
			rm(S_ij,P_ij,d_i,d_j,n)
		}
	}
	#N_ij_out
	#Nbar_out <- mean(N_ij_out)
	
	N_ij_in_vector <- N_ij_in[upper.tri(N_ij_in, diag=keep.diag)]
	N_ij_out_vector <- N_ij_out[upper.tri(N_ij_out, diag=keep.diag)]
	
	out$Nbar <- mean(c(N_ij_in_vector, N_ij_out_vector))
	out$mod  <- mean(c(sd(N_ij_in_vector), sd(N_ij_out_vector)))
	out$Nbar_higher <- mean(N_ij_in_vector)
	out$Nbar_lower <- mean(N_ij_out_vector)
 	out$mod_lower <- sd(N_ij_out_vector)
	out$mod_higher <- sd(N_ij_in_vector)
	if (keep.Nij) out$"N_ij_higher" <- N_ij_in + t(N_ij_in) - diag(1, nrow(N_ij_in)) - {if (keep.diag==FALSE) diag(1, nrow(N_ij_in)) else diag(0, nrow(N_ij_in))}
	if (keep.Nij) out$"N_ij_lower" <- N_ij_out + t(N_ij_out) - diag(1, nrow(N_ij_out)) - {if (keep.diag==FALSE) diag(1, nrow(N_ij_out)) else diag(0, nrow(N_ij_out))}

	return(out)
}	
