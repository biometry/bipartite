null.distr <- function(N, web, distr="lognormal"){
    
    null.lnorm.one <- function(web, distr){
        lmeansd <- function(x){ 
          N <- length(x)
          mu <- mean(log(x))    
          sigma2 <- sqrt( (N-1)/N) * sd(log(x)) # correction for sample to population
          c(meanlog=mu, sdlog=sigma2)
        }
                    
        if (distr == "log-normal" | distr == "lognormal") {
            rpar <- lmeansd(rowSums(web))
            cpar <- lmeansd(colSums(web))
            rnew <- rlnorm(NROW(web), rpar[1], rpar[2])
            cnew <- rlnorm(NCOL(web), cpar[1], cpar[2])    
        }

        if (distr =="negative binomial" | distr == "negbin"){
            rnbpars <- fitdistr(rowSums(web), "negative binomial")$estimate
            cnbpars <- fitdistr(colSums(web), "negative binomial")$estimate
            rnew <- rnbinom(NROW(web), size=rnbpars[1], mu=rnbpars[2]) + 1
            cnew <- rnbinom(NCOL(web), size=cnbpars[1], mu=cnbpars[2]) + 1           
        }
        
        Pmat <- tcrossprod(rnew/sum(rnew), cnew/sum(cnew))
        
        newweb <- unname(web)
        newweb[] <- 0
        
        simulated.interactions <- sample(1:prod(dim(web)), sum(web), prob=as.vector(Pmat), replace=TRUE)
        sim.ints.tabled <- as.data.frame(table(simulated.interactions))  # table makes first column a factor!
        web.index <- as.numeric(levels(sim.ints.tabled[,1]))
        newweb[web.index] <- sim.ints.tabled[, 2]
        
        empty(newweb)
    }  
    #
    replicate(n=N, expr=null.lnorm.one(web, distr=distr), simplify=FALSE)
}
