`genweb` <-
function(N1=10, N2=30, dens=2) {
  # generates a random web; N1=nrows bzw. plant species (or basal); ncols (number of "top" species); dens is also called SI (mean frequency per cell)
  if(length(N1)==2)  {N2 = N1[2]; N1 = N1[1]}
  mm = round(N1*N2*dens)  # total number of interactions in the web
  for (i in 1:2) {        # this loops through (row,col)     # N=10      # for tests
    N = N1 ; if (i==2) N = N2
    M = mm/N                     # mean number of interactions per species on level i on "real" scale, i.e. untransformed
    sigma = 1.5                  # standard deviation on log-scale, i.e. the "normal distribution scale"; this value is the median of the NCEAS-webs (both for rows and columns!)
    mu = log(M) - 0.5 *sigma^2   # mean on the log-scale ; thanks to Thomas Hovestadt, who insisted on this formula  # E(X) = exp(mu + 1/2 s^2) 
    ms = 1                           # ms = marginal sums ; 1 is only that ms is "known"
    if (mm==N) {ms=rep(1,N)} else {
      while (sum(ms) != mm){         # das dauert dann auch gerne mal ein paar Minuten!
      ms = sample(1:mm, size=N, replace=TRUE, prob=dlnorm(1:mm, meanlog=mu, sdlog=sigma))
      }
    }
    if (i==1) rs=ms else cs=ms
  }
  r2dtable(1, rs, cs)[[1]]
}

