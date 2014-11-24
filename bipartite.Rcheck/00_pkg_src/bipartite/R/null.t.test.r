null.t.test <- function(web, N=30, ...){
     # S L O W !! More repetitions will take longer ...
     #
     # A little null-model function to check, if the observed values actually are
     # much different to what one would expect under random numbers given the
     # observed row and column totals (i.e. information in the structure of the
     # web, not only in its species' abundances)
     #
     # web  a bipartite matrix with interactions
     # N    number of random webs to be created
     # ...  options passed on to other functions called internally:
     #             can go to t.test (for setting the conf.level to other than 0.95)
     #             or to index of networklevel or specieslevel (see there)
     #
     # returns a table with selected indices in rows and columns giving observed
     # values, mean null model values with CI-values, and t.test statistics (two-tailed)
     #
     # Does only work for indices directly returned by the function called (e.g. NOT for
     # degree distributions).

	 call.args <- list(...)
	 #return(call.args)
     obs <- networklevel(web, ...)
   	 null.list <- r2dtable(n=N, r=rowSums(web), c=colSums(web))
   	 if (length(call.args$index) == 1){ 
   	 		res <- matrix(sapply(null.list, networklevel, ...), nrow=1)
     } else {
     		res <- sapply(null.list, networklevel, ...)
     }
     # t.test the difference between observed and expected from null model:
     t.test.mat <- matrix(nrow=NROW(res), ncol=6)
     rownames(t.test.mat) <- if (NROW(res) == 1) call.args$index  else rownames(res)
     colnames(t.test.mat) <- c("obs", "null mean", "lower CI", "upper CI", "t", "P")
     for (i in 1:NROW(res)){
         t.res <- try(t.test(unlist(res[i,]), mu=obs[[i]], na.action=na.omit))#, ...
         t.test.mat[i,] <- c(obs[[i]], mean(unlist(res[i,]), na.rm=TRUE), t.res$conf.int,
                            t.res$statistic, t.res$p.value)
     }
     t.test.mat
}
#example:
#null.t.test(web, index=c("number of species", "number of links", "linkage density", "web asymmetry",
#          "number of compartments", "generality", "vulnerability", "interaction evenness",
#          "compartment diversity", "cluster coefficient", "H2", "ISA", "SA",
#          "extinction slope"), nrep=2, N=10)
